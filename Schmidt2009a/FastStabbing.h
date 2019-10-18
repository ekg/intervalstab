/************************************************************
Copyright (C) 2009 Jens M. Schmidt

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
or 3 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
************************************************************/

#ifdef _MSC_VER
#pragma once
#endif

#ifndef FAST_STABBING_H
#define FAST_STABBING_H

#include <ogdf/basic/Array.h>
#include <vector>
#include <ogdf/basic/List.h>
#include <ogdf/basic/SList.h>
#include <ogdf/basic/Stack.h>
#include <fstream>


namespace ogdf {

// an interval
struct interval {
	int l,r;
	interval* leftsibling;
	interval* rightchild;
	interval* parent;
	interval* smaller;
	ListIterator<interval*> pIt;
	#ifdef OGDF_DEBUG
	bool stabbed;
	#endif
};

// lexicographic order
inline bool operator <(const interval& x,const interval& y) {
	return (x.l < y.l || x.l == y.l && x.r < y.r);
}

// lexicographic order
inline bool operator >(const interval& x,const interval& y) {
	return y < x;
}

// output stream for intervals
inline ostream& operator<<(ostream& os, const ogdf::interval& a) {
	os << &a << "\t" << a.l << "\t" << a.r << "\tP " << a.parent << " L " << a.leftsibling
		<< " C " << a.rightchild << "  Sm " << a.smaller << "\n";
	return os;
}
inline ostream& operator<<(ostream& os, const std::vector<ogdf::interval*>& a) {
	for (unsigned int i = 0; i < a.size(); ++i) {
		os << *a[i];
	}
	return os;
}


// compare function for quicksort (lexicographic order)
class IntervalComparer: public Comparer<interval>
{
public:
	IntervalComparer() { }

	int compare(const interval& x, const interval& y) {
		if (x.l < y.l) return -1;
		if (x.l > y.l) return 1;
		if (x.r < y.r) return -1; // same starting point: ascending order
		if (x.r > y.r) return 1;
		return 0;
	}
};

// compare function for quicksort (lexicographic order with inverse end point order)
class IntervalComparerInv: public Comparer<interval>
{
public:
	IntervalComparerInv() { }

	int compare(const interval& x, const interval& y) {
		if (x.l < y.l) return -1;
		if (x.l > y.l) return 1;
		if (x.r < y.r) return 1; // same starting point: descending order
		if (x.r > y.r) return -1;
		return 0;
	}
};

// fast stabbing
class FastStabbing
{
private:
	Array<interval>& a; // array of intervals [0,n-1]
	int n,bigN;
	Array<ListPure<interval*> > eventlist; // sweepline
	Array<interval*> stop;
	interval dummy;

	void preprocessing();
	void preprocessingChazelle();
	bool verify(std::vector<interval*> output, const int& q);

public:
	FastStabbing(Array<interval>& intervals, int numberIntervals, int numberDomain)
			: a(intervals), eventlist(numberDomain+1), stop(0,numberDomain,NULL)  {
		n = numberIntervals;
		bigN = numberDomain;
		dummy.parent = NULL;
		dummy.leftsibling = NULL;
		dummy.rightchild = NULL;
		preprocessing();
	};
	void query(const int& q, std::vector<interval*>& output, const bool onlySearch, long& numComparisons);
};

// Chazelle stabbing
class ChazelleStabbing
{
private:
	Array<interval>& a; // array of intervals [0,n-1]
	int n,bigN;
	double delta;
	Array<ListPure<interval*> > eventlist; // sweepline
	struct window {
		int l;
		ListPure<interval*> intervals;
	};
	SListPure<window> windows;
	Array<SListIterator<window> > pWindow;

	void preprocessing();
	bool verify(std::vector<interval*> output, const int& q);

public:
	ChazelleStabbing(Array<interval>& intervals, int numberIntervals, int numberDomain, double d)
			: a(intervals), eventlist(numberDomain+1), pWindow(numberDomain+1)  {
		delta = d; // delta>0 indicates that chazelle is used
		n = numberIntervals;
		bigN = numberDomain;
		preprocessing();
	};
	void query(const int& q, std::vector<interval*>& output, const bool onlySearch, long& numComparisons);
};

}

#endif