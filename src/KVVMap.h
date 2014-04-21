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


#ifndef DBRX_KVVMAP_H
#define DBRX_KVVMAP_H

#include <utility>
#include <algorithm>
#include <stdexcept>


namespace dbrx {



template<typename It1, typename It2, bool wrKey, bool wrVal> class PairedIterators {
protected:
	std::pair<It1, It2> m_iterators;

	bool matchComp(bool cmpA, bool cmpB) const {
		if (cmpA == cmpB) return cmpA;
		else throw std::invalid_argument("Paired iterator comparison yields inhomogeneous result");		
	}

	void staticTypeCheck() {
		static_assert(std::is_same<typename It1::iterator_category, std::random_access_iterator_tag>::value, "Requires random access iterator");
		static_assert(std::is_same<typename It2::iterator_category, std::random_access_iterator_tag>::value, "Requires random access iterator");
		static_assert(std::is_same<typename It1::difference_type, typename It2::difference_type>::value, "Iterator difference type mismatch");
	}

public:
	// A custom pair class is required to support operations like std::sort
	// (given that both keys and values are writeable): Because std::sort calls
	// std::iter_swap with explicit namespace, we can't provide an ADL overload
	// for iter_swap. So we have to provide a custom swap for pairs of
	// references - as this can't be ADL'ed for the standard std::pair, we need
	// our own pair class.

	template<typename T1, typename T2> struct Pair {
		T1 first;
		T2 second;

		template<typename U1, typename U2> Pair& operator=(const Pair<U1, U2> &other) {
			first = other.first;
			second = other.second;
			return *this;
		}

		template<typename U1, typename U2> Pair& operator=(Pair<U1, U2> &&other) {
			first = std::move(other.first);
			second = std::move(other.second);
			return *this;
		}

		Pair() {}

		template<typename U1, typename U2> Pair(const Pair<U1, U2> &other)
			: first(other.first), second(other.second) {}

		template<typename U1, typename U2> Pair(Pair<U1, U2> &&other)
			: first(std::move(other.first)), second(std::move(other.second)) {}

		template<typename U1, typename U2> Pair(U1&& x, U2&& y)
			: first(std::forward<U1>(x)), second(std::forward<U2>(y)) {}
	};

	template<typename T1, typename T2> friend void swap(Pair<T1&, T2&> a, Pair<T1&, T2&> b) {
		using namespace std;
		swap(a.first, b.first);
		swap(a.second, b.second);
	}


	template<typename T1, typename T2>
		using PairT = typename std::conditional<wrKey && wrVal, Pair<T1, T2>, std::pair<T1, T2> >::type;

	using ValT1 = typename std::remove_reference<decltype(*m_iterators.first)>::type;
	using ValT2 = typename std::remove_reference<decltype(*m_iterators.second)>::type;

	using RefPair = Pair<ValT1&, ValT2&>;
	using ValPair = Pair<ValT1, ValT2>;

	using T1Ref = typename std::conditional<wrKey, ValT1&, const ValT1&>::type;
	using T1Val = typename std::conditional<wrKey, ValT1, const ValT1>::type;
	using T2Ref = typename std::conditional<wrVal, ValT2&, const ValT2&>::type;
	using T2Val = typename std::conditional<wrVal, ValT2, const ValT2>::type;

	using reference = PairT<T1Ref, T2Ref>;
	using value_type = PairT<T1Val, T2Val>;
	using pointer = PairT<typename It1::pointer, typename It2::pointer>;
	using difference_type = typename It1::difference_type;
	using iterator_category = typename It1::iterator_category;


	bool operator==(const PairedIterators &other) const {
		return matchComp(
			m_iterators.first == other.m_iterators.first,
			m_iterators.second == other.m_iterators.second
		);
	}

	bool operator!=(const PairedIterators &other) const {
		return matchComp(
			m_iterators.first != other.m_iterators.first,
			m_iterators.second != other.m_iterators.second
		);
	}

	bool operator<(const PairedIterators &other) const {
		return matchComp(
			m_iterators.first < other.m_iterators.first,
			m_iterators.second < other.m_iterators.second
		);
	}

	bool operator<=(const PairedIterators &other) const {
		return matchComp(
			m_iterators.first <= other.m_iterators.first,
			m_iterators.second <= other.m_iterators.second
		);
	}

	bool operator>(const PairedIterators &other) const {
		return matchComp(
			m_iterators.first > other.m_iterators.first,
			m_iterators.second > other.m_iterators.second
		);
	}

	bool operator>=(const PairedIterators &other) const {
		return matchComp(
			m_iterators.first >= other.m_iterators.first,
			m_iterators.second >= other.m_iterators.second
		);
	}

	difference_type operator-(const PairedIterators &other) const {
		typename It1::difference_type diffA = (m_iterators.first - other.m_iterators.first);
		typename It2::difference_type diffB = (m_iterators.second - other.m_iterators.second);
		if (diffA == diffB) return diffA;
		else throw std::invalid_argument("Paired iterator comparison yields inhomogeneous result");
	}


	PairedIterators& operator++() { ++m_iterators.first; ++m_iterators.second; return *this; }
	PairedIterators operator++(int) { decltype(*this) tmp; tmp = *this; operator++(); return tmp; }
	PairedIterators& operator--() { --m_iterators.first; --m_iterators.second; return *this; }
	PairedIterators operator--(int) { decltype(*this) tmp; tmp = *this; operator--(); return tmp; }
	PairedIterators& operator+=(difference_type n) { m_iterators.first += n;	m_iterators.second += n; return *this; }
	PairedIterators& operator-=(difference_type n) {	m_iterators.first -= n;	m_iterators.second -= n; return *this; }

	reference operator*() const { return reference(*m_iterators.first, *m_iterators.second); }

	reference operator[](difference_type n) const { return reference(m_iterators.first[n], m_iterators.second[n]); }

	PairedIterators() { staticTypeCheck(); }

	PairedIterators(It1 a, const It2 b) : m_iterators({std::move(a), std::move(b)}) { staticTypeCheck(); }

	~PairedIterators() {}

	friend PairedIterators operator+(PairedIterators a, difference_type n) { return a+=(n); }
	friend PairedIterators operator-(PairedIterators a, difference_type n) { return a-=(n); }
};



template<typename VecM> class ConstSizeVector {
public:
	using VecT = typename std::remove_reference<VecM>::type;

	using value_type = typename VecT::value_type;
	using reference = typename VecT::reference;
	using const_reference = typename VecT::const_reference;
	using pointer = typename VecT::pointer;
	using const_pointer = typename VecT::const_pointer;
	using iterator = typename VecT::iterator;
	using const_iterator = typename VecT::const_iterator;
	using reverse_iterator = typename VecT::reverse_iterator;
	using const_reverse_iterator = typename VecT::const_reverse_iterator;
	using difference_type = typename VecT::difference_type;
	using size_type = typename VecT::size_type;

protected:
	VecM m_vec;

public:
	bool empty() const noexcept { return m_vec.empty(); }
	size_type size() const noexcept { return m_vec.size(); }

	value_type data() noexcept { return m_vec.data(); }
	const value_type* data() const noexcept { return m_vec.data(); }

	iterator begin() noexcept { return m_vec.begin(); }
	const_iterator begin() const noexcept  { return m_vec.begin(); }
	const_iterator cbegin() const noexcept { return m_vec.cbegin(); }

	iterator end() noexcept { return m_vec.end(); }
	const_iterator end() const noexcept  { return m_vec.end(); }
	const_iterator cend() const noexcept { return m_vec.cbegin(); }

	reverse_iterator rbegin() { return m_vec.rbegin(); }
	const_reverse_iterator rbegin() const { return m_vec.rbegin(); }
	const_reverse_iterator crbegin() const { return m_vec.crbegin(); }

	reverse_iterator rend() { return m_vec.rend(); }
	const_reverse_iterator rend() const { return m_vec.rend(); }
	const_reverse_iterator crend() const { return m_vec.crend(); }

	reference front() { return m_vec.front(); }
	const_reference front() const { return m_vec.front(); }

	reference back() { return m_vec.back(); }
	const_reference back() const { return m_vec.back(); }

	reference operator[](size_t i) { return m_vec.operator[](i); }
	const_reference operator[](size_t i) const { return m_vec.operator[](i); }

	reference at(size_t i) { return m_vec.at(i); }
	const_reference at(size_t i) const { return m_vec.at(i); }

	template<typename T> ConstSizeVector& operator=(T&& other) {
		using namespace std;
		if (size() != other.size()) throw invalid_argument("Unequal size on assignment to ConstSizeVector");
		m_vec = std::forward<T>(other);
		return *this;
	}

	void swap(ConstSizeVector& other) noexcept {
		using namespace std;
		if (size() != other.size()) throw invalid_argument("Unequal size on swap of ConstSizeVector");
		swap(m_vec, other.m_vec);
	}

	ConstSizeVector() : m_vec() {}
	ConstSizeVector(VecT &vec) : m_vec(vec) {}

	friend void swap(ConstSizeVector& a, ConstSizeVector& b) noexcept { a.swap(b); }

	friend bool operator==(const ConstSizeVector& a, const ConstSizeVector& b)
		{ return (a.m_vec == b.m_vec); }

	friend bool operator!=(const ConstSizeVector& a, const ConstSizeVector& b)
		{ return (a.m_vec != b.m_vec); }
};



template<
	typename KeysM, typename ValsM,
	typename Compare = std::less<typename std::remove_reference<decltype(*std::declval<KeysM>().begin())>::type>
> class KVVMap {
public:
	using KeysT = typename std::remove_reference<KeysM>::type;
	using ValsT = typename std::remove_reference<ValsM>::type;

	using Keys = ConstSizeVector<KeysT&>;
	using Values = ConstSizeVector<ValsT&>;

protected:
	KeysM m_keys;
	ValsM m_vals;
	Compare m_comp;

	Keys m_constSizeKeys;
	Values m_constSizeVals;

public:
	using key_type = typename std::remove_reference<decltype(*m_keys.begin())>::type;
	using mapped_type = typename std::remove_reference<decltype(*m_vals.begin())>::type;
	using value_type = std::pair<const key_type, mapped_type>;
	using size_type = typename std::remove_reference<decltype(m_keys.size())>::type;
	using difference_type = typename std::remove_reference<decltype(m_keys.end() - m_keys.begin())>::type;
	using key_compare = Compare;
	// using allocator_type = ...
	using reference	= std::pair<const key_type&, mapped_type&>;
	using const_reference	= const std::pair<const key_type&, const mapped_type&>;
	// using pointer = std::pair<key_type*, mapped_type*>;
	// using const_pointer = const std::pair<const key_type*, const mapped_type*>;

	using const_iterator = PairedIterators<typename KeysT::const_iterator, typename ValsT::const_iterator, false, false>;
	using iterator = PairedIterators<typename KeysT::const_iterator, typename ValsT::iterator, false, true>;
	// using reverse_iterator = ...
	// const_reverse_iterator = ...
	using sort_iterator = PairedIterators<typename KeysT::iterator, typename ValsT::iterator, true, true>;

protected:
	struct PairComp {
		const Compare &comp;

		template<typename T1, typename T2> using Pair = typename sort_iterator::template Pair<T1, T2>;

		template<typename T1, typename T2, typename U1, typename U2>
			bool operator()(const Pair<T1, T2> &a, const Pair<U1, U2> &b)
				{ return comp(a.first, b.first); };

		PairComp(const Compare& c): comp(c) {}
	};

public:
	const Keys& keys() const { return m_constSizeKeys; }

	const Values& values() const { return m_constSizeVals; }
	Values& values() { return m_constSizeVals; }

	iterator begin() { return iterator(m_keys.begin(), m_vals.begin()); }
	const_iterator begin() const { return const_iterator(m_keys.begin(), m_vals.begin()); }
	iterator end() { return iterator(m_keys.end(), m_vals.end()); }
	const_iterator end() const { return const_iterator(m_keys.end(), m_vals.end()); }

	const_iterator cbegin() const { return const_iterator(m_keys.cbegin(), m_vals.cbegin()); }
	const_iterator cend() const { return const_iterator(m_keys.cend(), m_vals.cend()); }

	// iterator rbegin() { return iterator(m_keys.rbegin(), m_vals.rbegin()); }
	// const_iterator rbegin() const { return const_iterator(m_keys.rbegin(), m_vals.rbegin()); }
	// iterator rend() { return iterator(m_keys.rend(), m_vals.rend()); }
	// const_iterator rend() const { return const_iterator(m_keys.rend(), m_vals.rend()); }

	// const_iterator crbegin() const { return const_iterator(m_keys.crbegin(), m_vals.crbegin()); }
	// const_iterator crend() const { return const_iterator(m_keys.crend(), m_vals.crend()); }

	sort_iterator sort_begin() { return sort_iterator(m_keys.begin(), m_vals.begin()); }
	sort_iterator sort_end() { return sort_iterator(m_keys.end(), m_vals.end()); }


protected:
	template<typename KT> mapped_type& getValue(KT&& key) {
		auto r = std::lower_bound(m_keys.begin(), m_keys.end(), std::forward<KT>(key), m_comp);
		if (r == m_keys.end()) {
			m_keys.push_back(std::forward<KT>(key)); m_vals.push_back(mapped_type());
			return m_vals.back();
		} else {
			typename ValsT::size_type idx = r - m_keys.begin();
			if (*r != std::forward<KT>(key)) {
				m_keys.insert(r, std::forward<KT>(key));
				m_vals.insert(m_vals.begin() + idx, mapped_type());
			}
			return m_vals[idx];
		}
	}


public:
	mapped_type& at(const key_type& key) {
		auto r = std::lower_bound(m_keys.begin(), m_keys.end(), key, m_comp);
		if ( (r == m_keys.end()) || (*r != key) ) throw std::out_of_range("");
		else return m_vals[r - m_keys.begin()];
	}

	const mapped_type& at(const key_type& key) const {
		auto r = std::lower_bound(m_keys.begin(), m_keys.end(), key, m_comp);
		if ( (r == m_keys.end()) || (*r != key) ) throw std::out_of_range("");
		else return m_vals[r - m_keys.begin()];
	}


	iterator find(const key_type& key) {
		iterator r = std::lower_bound(begin(), end(), key,
			[&](reference p, const key_type& k) { return m_comp(p.first, k); });
		if ((*r).first != key) throw std::out_of_range("");
		else return r;
	}

	const_iterator find(const key_type& key) const {
		const_iterator r = std::lower_bound(begin(), end(), key,
			[&](const_reference p, const key_type& k) { return m_comp(p.first, k); });
		if ((*r).first != key) throw std::out_of_range("");
		else return r;
	}


	mapped_type& operator[](const key_type& key) { return getValue(key); }

	mapped_type& operator[](key_type&& key) { return getValue(std::move(key)); }


	void normalize() { std::sort(sort_begin(), sort_end(), PairComp(m_comp)); }

	bool empty() const { return m_keys.empty(); }
	size_type size() const { return m_keys.size(); }
	size_type max_size() const { return m_keys.max_size(); }

	void clear() { m_keys.clear(); m_vals.clear(); }


	template<typename P> std::pair<iterator,bool> insert(P&& value) {
		if (empty() || (std::forward<P>(value).first > m_keys.back())) {
			m_keys.push_back(std::forward<P>(value).first);
			m_vals.push_back(std::forward<P>(value).second);
			return {end() - 1, true};
		} else {
			auto keyIt = std::lower_bound(m_keys.begin(), m_keys.end(), std::forward<P>(value).first, m_comp);
			size_type idx = keyIt - m_keys.begin();
			auto resIt = begin() + idx;
			if (keyIt == m_keys.end()) {
				m_keys.insert(keyIt, std::forward<P>(value).first);
				m_vals.insert(m_vals.begin() + idx, std::forward<P>(value).second);
				return {resIt, true};
			}
			else return {resIt, false};
		}
	}

	std::pair<iterator, bool> insert(const value_type& value)
		{ return insert<value_type&>(value); }

	// iterator insert(const_iterator hint, const value_type& value);

	// template<typename P> iterator insert( const_iterator hint, P&& value);

	// template<typename InputIt> void insert(InputIt first, InputIt last);

	// void insert(std::initializer_list<value_type> ilist);


	// template<typename... Args> std::pair<iterator, bool> emplace(Args&&... args);

	// template <typename... Args> iterator emplace_hint(const_iterator hint, Args&&... args);


	// iterator erase( const_iterator pos );

	// iterator erase( const_iterator first, const_iterator last );

	// size_type erase( const key_type& key );


	void swap(KVVMap& other)
		{ std::swap(m_keys, other.m_keys); std::swap(m_vals, other.m_vals); std::swap(m_comp, other.m_comp); }


	void operator=(const KVVMap& other)
		{ m_keys = other.m_keys; m_vals = other.m_vals; m_comp = other.m_comp; }

	void operator=(KVVMap&& other)
		{ m_keys = std::move(other.m_keys); m_vals = std::move(other.m_vals); m_comp = std::move(other.m_comp); }


	KVVMap(const KVVMap& other)
		: m_keys(other.m_keys), m_vals(other.m_vals), m_comp(other.m_comp),
		  m_constSizeKeys(m_keys), m_constSizeVals(m_vals) {}

	KVVMap(KVVMap&& other)
		: m_keys(std::move(other.m_keys)), m_vals(std::move(other.m_vals)), m_comp(std::move(other.m_comp)),
		  m_constSizeKeys(m_keys), m_constSizeVals(m_vals) {}

	KVVMap(const Compare& comp = Compare())
		: m_comp(comp),	m_constSizeKeys(m_keys), m_constSizeVals(m_vals) {}

	template<typename KeysU, typename ValsU> KVVMap(KeysU&& k, ValsU&& v, const Compare& comp = Compare())
		: m_keys(std::forward<KeysU>(k)), m_vals(std::forward<ValsU>(v)), m_comp(comp),
		m_constSizeKeys(m_keys), m_constSizeVals(m_vals) {}


	friend void swap(KVVMap& a, KVVMap& b) noexcept { a.swap(b); }

	friend bool operator==(const KVVMap& a, const KVVMap& b)
		{ return (a.m_keys == b.m_keys) && (a.m_vals == b.m_vals); }

	friend bool operator!=(const KVVMap& a, const KVVMap& b)
		{ return (a.m_keys != b.m_keys) || (a.m_vals != b.m_vals); }
};


} // namespace dbrx


#endif // DBRX_KVVMAP_H
