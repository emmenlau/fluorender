/**
 * This file is used for compatibility across windows and mac/linux platforms.
 * @author Brig Bagley
 * @version 4 March 2014
 */
#ifdef _WIN32
#include "WacUtils/WacUtils.h"
#include <windows.h>
#include <ole2.h>
#include <wintab.h>
#include <pktdef.h>
#define WSTOD(s)                                _wtof(s)
#define WSTOI(s)                                _wtoi(s)
#define CREATE_DIR(f)                           CreateDirectory(f, NULL)
#define SPRINTF                                 sprintf_s
#define SSCANF                                  scanf_s
#define STRDUP                                  _strdup
#define STRCPY(a,b,c)                           strcpy_s(a,b,c)
#define STRCAT(a,b,c)                           strcat_s(a,b,c)
#define STRNCPY(a,b,c,d)                        strncpy_s(a,b,c,d)
#define TIME                                    _time32

#ifndef WFOPEN__
#define WFOPEN__
#define WFOPEN                                  _wfopen_s
#endif

#ifndef FOPEN__
#define FOPEN__
#define FOPEN                                   fopen_s
#endif

#ifndef GET_TICK_COUNT__
#define GET_TICK_COUNT__
uint32_t GET_TICK_COUNT() { return GetTickCount(); }
#endif

#ifndef FIND_FILES__
#define FIND_FILES__
#include <vector>
void FIND_FILES(std::wstring m_path_name,
      std::wstring search_ext,
      std::vector<std::wstring> &m_batch_list,
      int &m_cur_batch) {
   std::wstring search_path = m_path_name.substr(0,
         m_path_name.find_last_of(L'\\')) + L'\\';
   std::wstring search_str = L"*" + search_ext;
   WIN32_FIND_DATAW FindFileData;
   HANDLE hFind;
   hFind = FindFirstFileW(search_str.c_str(), &FindFileData);
   if (hFind != INVALID_HANDLE_VALUE)
   {
      int cnt = 0;
      m_batch_list.clear();
      std::wstring name = search_path + FindFileData.cFileName;
      m_batch_list.push_back(name);
      if (name == m_path_name)
         m_cur_batch = cnt;
      cnt++;

      while (FindNextFileW(hFind, &FindFileData) != 0)
      {
         name = search_path + FindFileData.cFileName;
         m_batch_list.push_back(name);
         if (name == m_path_name)
            m_cur_batch = cnt;
         cnt++;
      }
   }
   FindClose(hFind);
}
#endif

//#define GET_WM_ACTIVATE_STATE(wp, lp)           LOWORD(wp)
//#define GET_WM_COMMAND_ID(wp, lp)               LOWORD(wp)
//#define GET_WM_COMMAND_HWND(wp, lp)             (HWND)(lp)
//#define GET_WM_COMMAND_CMD(wp, lp)              HIWORD(wp)
//#define FORWARD_WM_COMMAND(hwnd, id, hwndCtl, codeNotify, fn) \
//   (void)(fn)((hwnd), WM_COMMAND, MAKEWPARAM((UINT)(id),(UINT)(codeNotify)), (LPARAM)(HWND)(hwndCtl))
/* -------------------------------------------------------------------------- */
#else

#include <cstdlib>
#include <cstdio>
#define WSTOD(s)                                   atof(s)
#define WSTOI(s)                                   atoi(s)
#define DWORD                                      unsigned int
#define CREATE_DIR(f)                              mkdir(f, S_IRWXU | S_IRGRP | S_IXGRP)
#define SetDoubleBuffered(t)
#define SPRINTF                                    sprintf
#define SSCANF                                     scanf
#define STRDUP                                     strdup
#define STRCPY(a,b,c)                              strcpy(a,c)
#define STRCAT(a,b,c)                              strcat(a,c)
#define STRNCPY(a,b,c,d)                           strncpy(a,c,d)
#define TIME                                       time

#ifndef FIND_FILES__
#define FIND_FILES__
#include <dirent.h>
#include <string>
#include <vector>
#include <codecvt>
#include <iostream>

std::wstring s2ws(std::string& str) {
   typedef std::codecvt_utf8<wchar_t> convert_typeX;
   std::wstring_convert<convert_typeX, wchar_t> converterX;
   return converterX.from_bytes(str);
}

std::string ws2s(std::wstring& str) {
   typedef std::codecvt_utf8<wchar_t> convert_typeX;
   std::wstring_convert<convert_typeX, wchar_t> converterX;
   return converterX.to_bytes(str);
}

void FIND_FILES(std::wstring m_path_name,
      std::wstring search_ext,
      std::vector<std::wstring> &m_batch_list,
      int &m_cur_batch) {
   std::wstring search_path = m_path_name.substr(0,m_path_name.find_last_of(L'/')) + L'/';
   std::wstring search_str(L".txt");
   DIR* dir;
   struct dirent *ent;
   if ((dir = opendir(ws2s(search_path).c_str())) != NULL) {
      int cnt = 0;
      m_batch_list.clear();
      while((ent = readdir(dir)) != NULL) {
         //check if it contains the string.
         if (strstr(ent->d_name,ws2s(search_str).c_str())) {
            std::string ss = ent->d_name;
            std::wstring f = s2ws(ss);
            std::wstring name = search_path + f;
            m_batch_list.push_back(name);
            if (name == m_path_name)
               m_cur_batch = cnt;
            cnt++;
         }
      }
   }
}
#endif

#ifndef WFOPEN__
#define WFOPEN__
FILE* WFOPEN(FILE ** fp, const wchar_t* filename, const wchar_t* mode) {
   *fp = fopen((const char*)filename,(const char*)mode);
   return *fp;
}
#endif

#ifndef FOPEN__
#define FOPEN__
FILE* FOPEN(FILE ** fp, const char* filename, const char* mode) {
   *fp = fopen(filename,mode);
   return *fp;
}
#endif

#ifndef GET_TICK_COUNT__
#define GET_TICK_COUNT__
#include <sys/time.h>
uint32_t GET_TICK_COUNT() {
   struct timeval ts;
   gettimeofday(&ts, NULL);
   return ts.tv_sec * 1000 + ts.tv_usec / 1000;
}
#endif

#ifdef _LINUX
#endif

#ifdef _DARWIN
#endif

//#define GET_WM_ACTIVATE_STATE(wp, lp)               (wp)
//#define GET_WM_COMMAND_ID(wp, lp)                   (wp)
//#define GET_WM_COMMAND_HWND(wp, lp)                 (HWND)LOWORD(lp)
//#define GET_WM_COMMAND_CMD(wp, lp)                  HIWORD(lp)
//#define FORWARD_WM_COMMAND(hwnd, id, hwndCtl, codeNotify, fn) \
//   (void)(fn)((hwnd), WM_COMMAND, (WPARAM)(int)(id), MAKELPARAM((UINT)(hwndCtl), (codeNotify)))
/* -------------------------------------------------------------------------- */
#endif