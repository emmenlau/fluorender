#include "brkxml_reader.h"
#include "FLIVR/TextureRenderer.h"
#include "FLIVR/ShaderProgram.h"
#include "FLIVR/Utils.h"
//#include "FLIVR/Point.h"
#include "../compatibility.h"
#include <fstream>
#include <iostream>
#include <locale>

template <typename _T> void clear2DVector(std::vector<std::vector<_T>> &vec2d)
{
	if(vec2d.empty())return;
	for(int i = 0; i < vec2d.size(); i++)std::vector<_T>().swap(vec2d[i]);
	std::vector<std::vector<_T>>().swap(vec2d);
}

template <typename _T> inline void SafeDelete(_T* &p)
{
	if(p != NULL){
		delete (p);
		(p) = NULL;
	}
}

BRKXMLReader::BRKXMLReader()
{
   m_resize_type = 0;
   m_resample_type = 0;
   m_alignment = 0;

   m_level_num = 0;
   m_cur_level = -1;

   m_time_num = 0;
   m_cur_time = -1;
   m_chan_num = 0;
   m_slice_num = 0;
   m_x_size = 0;
   m_y_size = 0;
   
   m_valid_spc = false;
   m_xspc = 0.0;
   m_yspc = 0.0;
   m_zspc = 0.0;

   m_max_value = 0.0;
   m_scalar_scale = 1.0;

   m_batch = false;
   m_cur_batch = -1;
   m_file_type = BRICK_FILE_TYPE_NONE;
   
}

BRKXMLReader::~BRKXMLReader()
{
	Clear();
}

void BRKXMLReader::Clear()
{
	if(m_pyramid.empty()) return;
	for(int i = 0; i < m_pyramid.size(); i++){
		if(!m_pyramid[i].bricks.empty()){
			for(int j = 0; j < m_pyramid[i].bricks.size(); j++) SafeDelete(m_pyramid[i].bricks[j]);
			vector<BrickInfo *>().swap(m_pyramid[i].bricks);
		}
		if(!m_pyramid[i].filename.empty()){
			for(int j = 0; j < m_pyramid[i].filename.size(); j++){
				if(!m_pyramid[i].filename[j].empty()){
					for(int k = 0; k < m_pyramid[i].filename[j].size(); k++){
						if(!m_pyramid[i].filename[j][k].empty()){
							for(int m = 0; m < m_pyramid[i].filename[j][k].size(); m++)
								SafeDelete(m_pyramid[i].filename[j][k][m]);
							vector<wstring *>().swap(m_pyramid[i].filename[j][k]);
						}
					}
					vector<vector<wstring *>>().swap(m_pyramid[i].filename[j]);
				}
			}
			vector<vector<vector<wstring *>>>().swap(m_pyramid[i].filename);
		}
	}
	vector<LevelInfo>().swap(m_pyramid);
}

void BRKXMLReader::SetFile(string &file)
{
   if (!file.empty())
   {
      if (!m_path_name.empty())
         m_path_name.clear();
      m_path_name.assign(file.length(), L' ');
      copy(file.begin(), file.end(), m_path_name.begin());
#ifdef _WIN32
   wchar_t slash = L'\\';
#else
   wchar_t slash = L'/';
#endif
      m_data_name = m_path_name.substr(m_path_name.find_last_of(slash)+1);
	  m_dir_name = m_path_name.substr(0, m_path_name.find_last_of(slash)+1);
   }
   m_id_string = m_path_name;
}

void BRKXMLReader::SetFile(wstring &file)
{
   m_path_name = file;
#ifdef _WIN32
   wchar_t slash = L'\\';
#else
   wchar_t slash = L'/';
#endif
   m_data_name = m_path_name.substr(m_path_name.find_last_of(slash)+1);
   m_dir_name = m_path_name.substr(0, m_path_name.find_last_of(slash)+1);
   m_id_string = m_path_name;
}

void BRKXMLReader::Preprocess()
{
	long ival;
	double dval;
	
	wxString strValue;

	Clear();
	m_slice_num = 0;
	m_chan_num = 0;
	m_max_value = 0.0;

#ifdef _WIN32
	wchar_t slash = L'\\';
#else
	wchar_t slash = L'/';
#endif
	//separate path and name
	size_t pos = m_path_name.find_last_of(slash);
	if (pos == -1)
		return;
	wstring path = m_path_name.substr(0, pos+1);
	wstring name = m_path_name.substr(pos+1);

	wxString fname(m_path_name);
	if (m_doc.LoadFile(fname.ToStdString().c_str()) != 0){
		return;
	}
		
	tinyxml2::XMLElement *root = m_doc.RootElement();
	if (!root || strcmp(root->Name(), "BRK"))
		return;
	m_imageinfo = ReadImageInfo(root);

	ReadPyramid(root, m_pyramid);
	
	m_time_num = m_imageinfo.nFrame;
	m_chan_num = m_imageinfo.nChannel;

	m_cur_time = 0;

	if(m_pyramid.empty()) return;

	m_xspc = m_pyramid[m_pyramid.size()-1].xspc;
	m_yspc = m_pyramid[m_pyramid.size()-1].yspc;
	m_zspc = m_pyramid[m_pyramid.size()-1].zspc;

	m_x_size = m_pyramid[m_pyramid.size()-1].imageW;
	m_y_size = m_pyramid[m_pyramid.size()-1].imageH;
	m_slice_num = m_pyramid[m_pyramid.size()-1].imageD;

	m_file_type = m_pyramid[m_pyramid.size()-1].file_type;

	m_level_num = m_pyramid.size();
	m_cur_level = m_pyramid.size()-1;

	OutputInfo();
}

BRKXMLReader::ImageInfo BRKXMLReader::ReadImageInfo(tinyxml2::XMLElement *infoNode)
{
	ImageInfo iinfo;
	long ival;
	double dval;
	
	wxString strValue;
	
	strValue = infoNode->Attribute("nChannel");
	strValue.ToLong(&ival);
	iinfo.nChannel = ival;

	strValue = infoNode->Attribute("nFrame");
	strValue.ToLong(&ival);
	iinfo.nFrame = ival;

	strValue = infoNode->Attribute("nLevel");
	strValue.ToLong(&ival);
	iinfo.nLevel = ival;

	return iinfo;
}

void BRKXMLReader::ReadPyramid(tinyxml2::XMLElement *lvRootNode, vector<LevelInfo> &pylamid)
{
	wxString str;
	long ival;
	int level;

	tinyxml2::XMLElement *child = lvRootNode->FirstChildElement();
	while (child)
	{
		if (strcmp(child->Name(), "Level") == 0)
		{
			str = child->Attribute("lv");
			str.ToLong(&ival);
			level = ival;
			
			if(level + 1 > pylamid.size()) pylamid.resize(level + 1);
			ReadLevel(child, pylamid[level]);
		}
		child = child->NextSiblingElement();
	}
}

void BRKXMLReader::ReadLevel(tinyxml2::XMLElement* lvNode, LevelInfo &lvinfo)
{
	long ival;
	double dval;
	
	wxString strValue;

	strValue = lvNode->Attribute("imageW");
	strValue.ToLong(&ival);
	lvinfo.imageW = ival;

	strValue = lvNode->Attribute("imageH");
	strValue.ToLong(&ival);
	lvinfo.imageH = ival;

	strValue = lvNode->Attribute("imageD");
	strValue.ToLong(&ival);
	lvinfo.imageD = ival;

	strValue = lvNode->Attribute("xspc");
	strValue.ToDouble(&dval);
	lvinfo.xspc = dval;

	strValue = lvNode->Attribute("yspc");
	strValue.ToDouble(&dval);
	lvinfo.yspc = dval;

	strValue = lvNode->Attribute("zspc");
	strValue.ToDouble(&dval);
	lvinfo.zspc = dval;

	strValue = lvNode->Attribute("bitDepth");
	strValue.ToLong(&ival);
	lvinfo.bit_depth = ival;

	strValue = lvNode->Attribute("FileType");
	if(strValue == "RAW") lvinfo.file_type = BRICK_FILE_TYPE_RAW;
	else if (strValue == "JPEG") lvinfo.file_type = BRICK_FILE_TYPE_JPEG;
	else lvinfo.file_type = BRICK_FILE_TYPE_NONE;

	tinyxml2::XMLElement *child = lvNode->FirstChildElement();
	while (child)
	{
		if (strcmp(child->Name(), "Bricks") == 0){
			strValue = child->Attribute("brick_baseW");
			strValue.ToLong(&ival);
			lvinfo.brick_baseW = ival;

			strValue = child->Attribute("brick_baseH");
			strValue.ToLong(&ival);
			lvinfo.brick_baseH = ival;

			strValue = child->Attribute("brick_baseD");
			strValue.ToLong(&ival);
			lvinfo.brick_baseD = ival;

			ReadPackedBricks(child, lvinfo.bricks);
		}
		if (strcmp(child->Name(), "Files") == 0)  ReadFilenames(child, lvinfo.filename);
		child = child->NextSiblingElement();
	}
}

void BRKXMLReader::ReadPackedBricks(tinyxml2::XMLElement* packNode, vector<BrickInfo *> &brks)
{
	wxString str;
	long ival;
	int id;
	
	tinyxml2::XMLElement *child = packNode->FirstChildElement();
	while (child)
	{
		if (strcmp(child->Name(), "Brick") == 0)
		{
			str = child->Attribute("id");
			str.ToLong(&ival);
			id = ival;

			if(id + 1 > brks.size())
				brks.resize(id + 1, NULL);

			if (!brks[id]) brks[id] = new BrickInfo();
			ReadBrick(child, *brks[id]);
		}
		child = child->NextSiblingElement();
	}
}

void BRKXMLReader::ReadBrick(tinyxml2::XMLElement* brickNode, BrickInfo &binfo)
{
	long ival;
	double dval;
	
	wxString strValue;
		
	strValue = brickNode->Attribute("id");
	strValue.ToLong(&ival);
	binfo.id = ival;

	strValue = brickNode->Attribute("width");
	strValue.ToLong(&ival);
	binfo.x_size = ival;

	strValue = brickNode->Attribute("height");
	strValue.ToLong(&ival);
	binfo.y_size = ival;

	strValue = brickNode->Attribute("depth");
	strValue.ToLong(&ival);
	binfo.z_size = ival;

	strValue = brickNode->Attribute("st_x");
	strValue.ToLong(&ival);
	binfo.x_start = ival;

	strValue = brickNode->Attribute("st_y");
	strValue.ToLong(&ival);
	binfo.y_start = ival;

	strValue = brickNode->Attribute("st_z");
	strValue.ToLong(&ival);
	binfo.z_start = ival;

	strValue = brickNode->Attribute("offset");
	strValue.ToLong(&ival);
	binfo.offset = ival;

	strValue = brickNode->Attribute("size");
	strValue.ToLong(&ival);
	binfo.fsize = ival;

	tinyxml2::XMLElement *child = brickNode->FirstChildElement();
	while (child)
	{
		if (strcmp(child->Name(), "tbox") == 0)
		{
			Readbox(child, binfo.tx0, binfo.ty0, binfo.tz0, binfo.tx1, binfo.ty1, binfo.tz1);
		}
		else if (strcmp(child->Name(), "bbox") == 0)
		{
			Readbox(child, binfo.bx0, binfo.by0, binfo.bz0, binfo.bx1, binfo.by1, binfo.bz1);
		}
		child = child->NextSiblingElement();
	}
}

void BRKXMLReader::Readbox(tinyxml2::XMLElement* boxNode, double &x0, double &y0, double &z0, double &x1, double &y1, double &z1)
{
	double dval;
	
	wxString strValue;
	BrickInfo binfo;
	
	strValue = boxNode->Attribute("x0");
	strValue.ToDouble(&dval);
	x0 = dval;
	
	strValue = boxNode->Attribute("y0");
	strValue.ToDouble(&dval);
	y0 = dval;

	strValue = boxNode->Attribute("z0");
	strValue.ToDouble(&dval);
	z0 = dval;

	strValue = boxNode->Attribute("x1");
	strValue.ToDouble(&dval);
	x1 = dval;

	strValue = boxNode->Attribute("y1");
	strValue.ToDouble(&dval);
	y1 = dval;

	strValue = boxNode->Attribute("z1");
	strValue.ToDouble(&dval);
	z1 = dval;

}

void BRKXMLReader::ReadFilenames(tinyxml2::XMLElement* fileRootNode, vector<vector<vector<wstring *>>> &filename)
{
	wxString str;
	long ival;
	int frame, channel, id;

	tinyxml2::XMLElement *child = fileRootNode->FirstChildElement();
	while (child)
	{
		if (strcmp(child->Name(), "File") == 0)
		{
			str = child->Attribute("frame");
			str.ToLong(&ival);
			frame = ival;

			str = child->Attribute("channel");
			str.ToLong(&ival);
			channel = ival;

			str = child->Attribute("brickID");
			str.ToLong(&ival);
			id = ival;

			if(frame + 1 > filename.size())
				filename.resize(frame + 1);
			if(channel + 1 > filename[frame].size())
				filename[frame].resize(channel + 1);
			if(id + 1 > filename[frame][channel].size())
				filename[frame][channel].resize(id + 1, NULL);

			str = child->Attribute("filename");
			if (!filename[frame][channel][id]) filename[frame][channel][id] = new wstring(m_dir_name + str.ToStdWstring());
		}
		child = child->NextSiblingElement();
	}
}

void BRKXMLReader::SetSliceSeq(bool ss)
{
   //do nothing
}

bool BRKXMLReader::GetSliceSeq()
{
   return false;
}

void BRKXMLReader::SetTimeId(wstring &id)
{
   m_time_id = id;
}

wstring BRKXMLReader::GetTimeId()
{
   return m_time_id;
}

void BRKXMLReader::SetCurTime(int t)
{
	if (t < 0) m_cur_time = 0;
	else if (t >= m_imageinfo.nFrame) m_cur_time = m_imageinfo.nFrame - 1;
	else m_cur_time = t;
}

void BRKXMLReader::SetLevel(int lv)
{
	if(m_pyramid.empty()) return;
	
	if(lv < 0 || lv > m_level_num-1) return;
	m_cur_level = lv;

	m_xspc = m_pyramid[lv].xspc;
	m_yspc = m_pyramid[lv].yspc;
	m_zspc = m_pyramid[lv].zspc;

	m_x_size = m_pyramid[lv].imageW;
	m_y_size = m_pyramid[lv].imageH;
	m_slice_num = m_pyramid[lv].imageD;

	m_file_type = m_pyramid[lv].file_type;
}

void BRKXMLReader::SetBatch(bool batch)
{
#ifdef _WIN32
   wchar_t slash = L'\\';
#else
   wchar_t slash = L'/';
#endif
   if (batch)
   {
      //read the directory info
      wstring search_path = m_path_name.substr(0, m_path_name.find_last_of(slash)) + slash;
      FIND_FILES(search_path,L".oib",m_batch_list,m_cur_batch);
      m_batch = true;
   }
   else
      m_batch = false;
}

int BRKXMLReader::LoadBatch(int index)
{
   int result = -1;
   if (index>=0 && index<(int)m_batch_list.size())
   {
      m_path_name = m_batch_list[index];
      Preprocess();
      result = index;
      m_cur_batch = result;
   }
   else
      result = -1;

   return result;
}

int BRKXMLReader::LoadOffset(int offset)
{
   int result = m_cur_batch + offset;

   if (offset > 0)
   {
      if (result<(int)m_batch_list.size())
      {
         m_path_name = m_batch_list[result];
         Preprocess();
         m_cur_batch = result;
      }
      else if (m_cur_batch<(int)m_batch_list.size()-1)
      {
         result = (int)m_batch_list.size()-1;
         m_path_name = m_batch_list[result];
         Preprocess();
         m_cur_batch = result;
      }
      else
         result = -1;
   }
   else if (offset < 0)
   {
      if (result >= 0)
      {
         m_path_name = m_batch_list[result];
         Preprocess();
         m_cur_batch = result;
      }
      else if (m_cur_batch > 0)
      {
         result = 0;
         m_path_name = m_batch_list[result];
         Preprocess();
         m_cur_batch = result;
      }
      else
         result = -1;
   }
   else
      result = -1;

   return result;
}

double BRKXMLReader::GetExcitationWavelength(int chan)
{
   return 0.0;
}

//This function does not load image data into Nrrd.
Nrrd *BRKXMLReader::Convert(int t, int c, bool get_max)
{
	Nrrd *data = 0;

	if (t>=0 && t<m_time_num &&
		c>=0 && c<m_chan_num &&
		m_slice_num>0 &&
		m_x_size>0 &&
		m_y_size>0)
	{
		char dummy = 0;
		data = nrrdNew();
		if(m_pyramid[m_cur_level].bit_depth == 8)nrrdWrap(data, &dummy, nrrdTypeUChar, 3, (size_t)m_x_size, (size_t)m_y_size, (size_t)m_slice_num);
		else if(m_pyramid[m_cur_level].bit_depth == 16)nrrdWrap(data, &dummy, nrrdTypeUShort, 3, (size_t)m_x_size, (size_t)m_y_size, (size_t)m_slice_num);
		else if(m_pyramid[m_cur_level].bit_depth == 32)nrrdWrap(data, &dummy, nrrdTypeFloat, 3, (size_t)m_x_size, (size_t)m_y_size, (size_t)m_slice_num);
		nrrdAxisInfoSet(data, nrrdAxisInfoSpacing, m_xspc, m_yspc, m_zspc);
		nrrdAxisInfoSet(data, nrrdAxisInfoMax, m_xspc*m_x_size, m_yspc*m_y_size, m_zspc*m_slice_num);
		nrrdAxisInfoSet(data, nrrdAxisInfoMin, 0.0, 0.0, 0.0);
		nrrdAxisInfoSet(data, nrrdAxisInfoSize, (size_t)m_x_size, (size_t)m_y_size, (size_t)m_slice_num);
		data->data = NULL;//dangerous//
	}

	m_cur_time = t;
	//if (m_max_value > 0.0)
	//	m_scalar_scale = 65535.0 / m_max_value;
	m_scalar_scale = 1.0;

	if (m_xspc>0.0 && m_xspc<100.0 &&
		m_yspc>0.0 && m_yspc<100.0)
	{
		m_valid_spc = true;
		if (m_zspc<=0.0 || m_zspc>100.0)
			m_zspc = max(m_xspc, m_yspc);
	}
	else
	{
		m_valid_spc = false;
		m_xspc = 1.0;
		m_yspc = 1.0;
		m_zspc = 1.0;
	}

	return data;
}

wstring BRKXMLReader::GetCurName(int t, int c)
{
   return wstring(L"");
}

wstring BRKXMLReader::GetBrickFilePath(int fr, int ch, int id, int lv)
{
	int level = lv;
	int frame = fr;
	int channel = ch;
	int brickID = id;
	
	if(lv < 0 || lv >= m_level_num) level = m_cur_level;
	if(fr < 0 || fr >= m_time_num)  frame = m_cur_time;
	if(ch < 0 || ch >= m_chan_num)	channel = 0;
	if(id < 0 || id >= m_pyramid[level].bricks.size()) brickID = 0;
	
	return *m_pyramid[level].filename[frame][channel][brickID];
}

wstring BRKXMLReader::GetBrickFileName(int fr, int ch, int id, int lv)
{
	int level = lv;
	int frame = fr;
	int channel = ch;
	int brickID = id;
	
	if(lv < 0 || lv >= m_level_num) level = m_cur_level;
	if(fr < 0 || fr >= m_time_num)  frame = m_cur_time;
	if(ch < 0 || ch >= m_chan_num)	channel = 0;
	if(id < 0 || id >= m_pyramid[level].bricks.size()) brickID = 0;

	#ifdef _WIN32
	wchar_t slash = L'\\';
#else
	wchar_t slash = L'/';
#endif
	//separate path and name
	size_t pos = m_pyramid[level].filename[frame][channel][brickID]->find_last_of(slash);
	if (pos == -1)
		return wstring(L"");
	wstring name = m_pyramid[level].filename[frame][channel][brickID]->substr(pos+1);
	
	return name;
}

int BRKXMLReader::GetFileType(int lv)
{
	if(lv < 0 || lv > m_level_num-1) return m_file_type;

	return m_pyramid[lv].file_type;
}

void BRKXMLReader::OutputInfo()
{
	ofstream ofs;
	ofs.open("PyramidInfo.txt");

	ofs << "nChannel: " << m_imageinfo.nChannel << "\n";
	ofs << "nFrame: " << m_imageinfo.nFrame << "\n";
	ofs << "nLevel: " << m_imageinfo.nLevel << "\n\n";

	for(int i = 0; i < m_pyramid.size(); i++){
		ofs << "<Level: " << i << ">\n";
		ofs << "\timageW: " << m_pyramid[i].imageW << "\n";
		ofs << "\timageH: " << m_pyramid[i].imageH << "\n";
		ofs << "\timageD: " << m_pyramid[i].imageD << "\n";
		ofs << "\txspc: " << m_pyramid[i].xspc << "\n";
		ofs << "\tyspc: " << m_pyramid[i].yspc << "\n";
		ofs << "\tzspc: " << m_pyramid[i].zspc << "\n";
		ofs << "\tbrick_baseW: " << m_pyramid[i].brick_baseW << "\n";
		ofs << "\tbrick_baseH: " << m_pyramid[i].brick_baseH << "\n";
		ofs << "\tbrick_baseD: " << m_pyramid[i].brick_baseD << "\n";
		ofs << "\tbit_depth: " << m_pyramid[i].bit_depth << "\n";
		ofs << "\tfile_type: " << m_pyramid[i].file_type << "\n\n";
		for(int j = 0; j < m_pyramid[i].bricks.size(); j++){
			ofs << "\tBrick: " << " id = " <<  m_pyramid[i].bricks[j]->id
				<< " w = " << m_pyramid[i].bricks[j]->x_size
				<< " h = " << m_pyramid[i].bricks[j]->y_size
				<< " d = " << m_pyramid[i].bricks[j]->z_size
				<< " st_x = " << m_pyramid[i].bricks[j]->x_start
				<< " st_y = " << m_pyramid[i].bricks[j]->y_start
				<< " st_z = " << m_pyramid[i].bricks[j]->z_start
				<< " offset = " << m_pyramid[i].bricks[j]->offset
				<< " fsize = " << m_pyramid[i].bricks[j]->fsize << "\n";

			ofs << "\t\ttbox: "
				<< " x0 = " << m_pyramid[i].bricks[j]->tx0
				<< " y0 = " << m_pyramid[i].bricks[j]->ty0
				<< " z0 = " << m_pyramid[i].bricks[j]->tz0
				<< " x1 = " << m_pyramid[i].bricks[j]->tx1
				<< " y1 = " << m_pyramid[i].bricks[j]->ty1
				<< " z1 = " << m_pyramid[i].bricks[j]->tz1 << "\n";
			ofs << "\t\tbbox: "
				<< " x0 = " << m_pyramid[i].bricks[j]->bx0
				<< " y0 = " << m_pyramid[i].bricks[j]->by0
				<< " z0 = " << m_pyramid[i].bricks[j]->bz0
				<< " x1 = " << m_pyramid[i].bricks[j]->bx1
				<< " y1 = " << m_pyramid[i].bricks[j]->by1
				<< " z1 = " << m_pyramid[i].bricks[j]->bz1 << "\n";
		}
		ofs << "\n";
		for(int j = 0; j < m_pyramid[i].filename.size(); j++){
			for(int k = 0; k < m_pyramid[i].filename[j].size(); k++){
				for(int n = 0; n < m_pyramid[i].filename[j][k].size(); n++)
					ofs << "\t<Frame = " << j << " Channel = " << k << " ID = " << n << " Filepath = " << *m_pyramid[i].filename[j][k][n] << ">\n";
			}
		}
		ofs << "\n";
	}

	ofs.close();
}

void BRKXMLReader::build_bricks(vector<FLIVR::TextureBrick*> &tbrks, int lv)
{
	int lev;

	if(lv < 0 || lv > m_level_num-1) lev = m_cur_level;
	else lev = lv;

	// Initial brick size
	int bsize[3];

	bsize[0] = m_pyramid[lev].brick_baseW;
	bsize[1] = m_pyramid[lev].brick_baseH;
	bsize[2] = m_pyramid[lev].brick_baseD;

	bool force_pow2 = false;
	if (FLIVR::ShaderProgram::init())
		force_pow2 = !FLIVR::ShaderProgram::texture_non_power_of_two();

	int max_texture_size = 2048;
	if (FLIVR::ShaderProgram::init())
		max_texture_size = FLIVR::ShaderProgram::max_texture_size();

	int numb[1];
	if (m_pyramid[lev].bit_depth == 8 || m_pyramid[lev].bit_depth == 16 || m_pyramid[lev].bit_depth == 32)
		numb[0] = m_pyramid[lev].bit_depth / 8;
	else
		numb[0] = 0;
	
	//further determine the max texture size
//	if (FLIVR::TextureRenderer::get_mem_swap())
//	{
//		double data_size = double(m_pyramid[lev].imageW)*double(m_pyramid[lev].imageH)*double(m_pyramid[lev].imageD)*double(numb[0])/1.04e6;
//		if (data_size > FLIVR::TextureRenderer::get_mem_limit() ||
//			data_size > FLIVR::TextureRenderer::get_large_data_size())
//			max_texture_size = FLIVR::TextureRenderer::get_force_brick_size();
//	}
	
	if(bsize[0] > max_texture_size || bsize[1] > max_texture_size || bsize[2] > max_texture_size) return;
	if(force_pow2 && (FLIVR::Pow2(bsize[0]) > bsize[0] || FLIVR::Pow2(bsize[1]) > bsize[1] || FLIVR::Pow2(bsize[2]) > bsize[2])) return;
	
	if(!tbrks.empty())
	{
		for(int i = 0; i < tbrks.size(); i++)
		{
			tbrks[i]->freeBrkData();
			delete tbrks[i];
		}
		tbrks.clear();
	}
	vector<BrickInfo *>::iterator bite = m_pyramid[lev].bricks.begin();
	while (bite != m_pyramid[lev].bricks.end())
	{
		FLIVR::BBox tbox(FLIVR::Point((*bite)->tx0, (*bite)->ty0, (*bite)->tz0), FLIVR::Point((*bite)->tx1, (*bite)->ty1, (*bite)->tz1));
		FLIVR::BBox bbox(FLIVR::Point((*bite)->bx0, (*bite)->by0, (*bite)->bz0), FLIVR::Point((*bite)->bx1, (*bite)->by1, (*bite)->bz1));

		//numc? gm_nrrd?
		FLIVR::TextureBrick *b = new FLIVR::TextureBrick(0, 0, (*bite)->x_size, (*bite)->y_size, (*bite)->z_size, 1, numb, 
														 (*bite)->x_start, (*bite)->y_start, (*bite)->z_start,
														 (*bite)->x_size, (*bite)->y_size, (*bite)->z_size, bbox, tbox, (*bite)->id, (*bite)->offset, (*bite)->fsize);
		tbrks.push_back(b);
		
		bite++;
	}

	return;
}

void BRKXMLReader::build_pyramid(vector<FLIVR::Pyramid_Level> &pyramid, vector<vector<vector<vector<wstring *>>>> &filenames, int t, int c)
{
	if (!pyramid.empty())
	{
		for (int i = 0; i < pyramid.size(); i++)
		{
			if (pyramid[i].data) nrrdNix(pyramid[i].data);
			for (int j = 0; j < pyramid[i].bricks.size(); j++)
				if (pyramid[i].bricks[j]) delete pyramid[i].bricks[j];
		}
		vector<FLIVR::Pyramid_Level>().swap(pyramid);
	}

	if(!filenames.empty())
	{
		for (int i = 0; i < filenames.size(); i++)
			for (int j = 0; j < filenames[i].size(); j++)
				for (int k = 0; k < filenames[i][j].size(); k++)
					for (int n = 0; n < filenames[i][j][k].size(); n++)
					if (filenames[i][j][k][n]) delete filenames[i][j][k][n];
		vector<vector<vector<vector<wstring *>>>>().swap(filenames);
	}

	pyramid.resize(m_pyramid.size());

	for (int i = 0; i < m_pyramid.size(); i++)
	{
		SetLevel(i);
		pyramid[i].data = Convert(t, c, false);
		build_bricks(pyramid[i].bricks);
		pyramid[i].filenames = &m_pyramid[i].filename[t][c];
		pyramid[i].filetype = GetFileType();
	}

	filenames.resize(m_pyramid.size());
	for (int i = 0; i < filenames.size(); i++)
	{
		filenames[i].resize(m_pyramid[i].filename.size());
		for (int j = 0; j < filenames[i].size(); j++)
		{
			filenames[i][j].resize(m_pyramid[i].filename[j].size());
			for (int k = 0; k < filenames[i][j].size(); k++)
			{
				filenames[i][j][k].resize(m_pyramid[i].filename[j][k].size());
				for (int n = 0; n < filenames[i][j][k].size(); n++)
				{
					filenames[i][j][k][n] = new wstring(*m_pyramid[i].filename[j][k][n]);
				}
			}
		}
	}

}

