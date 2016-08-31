/*
For more information, please see: http://software.sci.utah.edu

The MIT License

Copyright (c) 2014 Scientific Computing and Imaging Institute,
University of Utah.


Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
*/
#ifndef FL_VolCache_h
#define FL_VolCache_h

#include <deque>
#include <boost/signals2.hpp>

namespace FL
{
	struct VolCache
	{
		VolCache() :
			nrrd_data(0),
			nrrd_label(0),
			data(0),
			label(0),
			frame(0),
			valid(false),
			modified(false) {}

		//manage memory externally
		//don't care releasing here
		//handles to release
		void* nrrd_data;
		void* nrrd_label;
		//actual data
		void* data;
		void* label;
		size_t frame;
		bool valid;
		bool modified;
	};

	//queue with a max size
	class CacheQueue
	{
	public:
		CacheQueue():
		m_max_size(1) {};
		~CacheQueue();

		inline void set_max_size(size_t size);
		inline size_t get_max_size();
		inline size_t size();
		inline void push_back(const VolCache& val);
		inline VolCache get(size_t frame);
		inline void set_modified(size_t frame, bool value = true);

		//external calls
		boost::signals2::signal<void(VolCache&)> m_new_cache;
		boost::signals2::signal<void(VolCache&)> m_del_cache;

	private:
		size_t m_max_size;
		std::deque<VolCache> m_queue;
	};

	inline CacheQueue::~CacheQueue()
	{
		for (size_t i = 0; i < m_queue.size(); ++i)
			m_del_cache(m_queue[i]);
	}

	inline void CacheQueue::set_max_size(size_t size)
	{
		if (size == 0)
			return;
		m_max_size = size;
		while (m_queue.size() > m_max_size)
		{
			m_del_cache(m_queue.front());
			m_queue.pop_front();
		}
	}

	inline size_t CacheQueue::get_max_size()
	{
		return m_max_size;
	}

	inline size_t CacheQueue::size()
	{
		return m_queue.size();
	}

	inline void CacheQueue::push_back(const VolCache& val)
	{
		while (m_queue.size() > m_max_size - 1)
		{
			m_del_cache(m_queue.front());
			m_queue.pop_front();
		}
		m_queue.push_back(val);
	}

	inline VolCache CacheQueue::get(size_t frame)
	{
		bool found = false;
		int index = -1;
		for (size_t i = 0; i < m_queue.size(); ++i)
		{
			if (m_queue[i].frame == frame)
			{
				index = (int)i;
				found = true;
				break;
			}
		}
		if (found)
		{
			if (!m_queue[index].valid)
				m_new_cache(m_queue[index]);
			return m_queue[index];
		}
		else
		{
			VolCache vol_cache;
			vol_cache.frame = frame;
			m_new_cache(vol_cache);
			push_back(vol_cache);
			return vol_cache;
		}
	}

	inline void CacheQueue::set_modified(size_t frame, bool value)
	{
		for (size_t i = 0; i < m_queue.size(); ++i)
		{
			if (m_queue[i].frame == frame)
			{
				if (m_queue[i].valid)
					m_queue[i].modified = value;
				return;
			}
		}
	}

}//namespace FL

#endif//FL_VolCache_h