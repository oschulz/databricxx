// Copyright (C) 2014 Oliver Schulz <oschulz@mpp.mpg.de>

// This is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or
// (at your option) any later version.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.


#ifndef DBRX_NAMETABLE_H
#define DBRX_NAMETABLE_H

#include <atomic>
#include <memory>

#include "Name.h"


namespace dbrx {


class GlobalObject final {
public:
	template<typename T> using Ptr = std::atomic<T*>;

	template<typename T> class Deleter final {
	protected:
		Ptr<T>& m_ptr;

	public:
		Deleter(Ptr<T> &ptr) : m_ptr(ptr) {}

		~Deleter() {
			T* value = m_ptr.load();
			if ( (value != nullptr) && m_ptr.compare_exchange_strong(value, nullptr) )
				delete value;
		}

	};

	template<typename T> static T& get(std::atomic<T*> &s_globalValue) {
		std::unique_ptr<NameTable> emptyPtr;
		T* value = s_globalValue.load();
		if (value != nullptr) {
			return *value;
		} else {
			s_globalValue.compare_exchange_strong(value, new T);
			value = s_globalValue.load();
			return *value;
		}
	}
};


class NameTable {
protected:
	struct Internals;
	Internals *m_internals;
public:
	static GlobalObject::Ptr<NameTable> s_globalNameTable;

	static NameTable& global() { return GlobalObject::get(s_globalNameTable); }

	Name resolve(const std::string &s);

	using StringContent = std::string;

	NameTable();
	virtual ~NameTable();
};


} // namespace dbrx

#endif // DBRX_NAMETABLE_H
