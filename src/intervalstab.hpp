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

#pragma once

#include <vector>
#include <list>
#include <stack>
#include <fstream>
#include <iostream>
#include <cassert>

namespace intervalstab {

// an interval
struct interval {
	uint64_t l = 0;
    uint64_t r = 0;
	interval* leftsibling = nullptr;
	interval* rightchild = nullptr;
	interval* parent = nullptr;
	interval* smaller = nullptr;
    std::list<interval*>::iterator pIt = (std::list<interval*>::iterator)nullptr;
    bool stabbed = false;
    interval(void) { }
    interval(const uint64_t& a,
             const uint64_t& b)
        : l(a),
          r(b) { }
    interval(const interval& other)
        : l(other.l)
        , r(other.r)
        , leftsibling(other.leftsibling)
        , rightchild(other.rightchild)
        , parent(other.parent)
        , smaller(other.smaller)
        , pIt(other.pIt)
        , stabbed(other.stabbed)
        { }
    interval(interval&& other) {
        l = other.l;
        r = other.r;
        leftsibling = other.leftsibling;
        rightchild = other.rightchild;
        parent = other.parent;
        smaller = other.smaller;
        pIt = other.pIt;
        stabbed = other.stabbed;
    }
    interval& operator=(const interval& other) {
        l = other.l;
        r = other.r;
        leftsibling = other.leftsibling;
        rightchild = other.rightchild;
        parent = other.parent;
        smaller = other.smaller;
        pIt = other.pIt;
        stabbed = other.stabbed;
        return *this;
    }
    interval& operator=(interval&& other) {
        l = other.l;
        r = other.r;
        leftsibling = other.leftsibling;
        rightchild = other.rightchild;
        parent = other.parent;
        smaller = other.smaller;
        pIt = other.pIt;
        stabbed = other.stabbed;
        return *this;
    }
    ~interval(void) { }
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
inline std::ostream& operator<<(std::ostream& os, const interval& a) {
	os << &a << "\t" << a.l << "\t" << a.r << "\tP " << a.parent << " L " << a.leftsibling
		<< " C " << a.rightchild << "  Sm " << a.smaller << "\n";
	return os;
}
inline std::ostream& operator<<(std::ostream& os, const std::vector<interval*>& a) {
	for (unsigned int i = 0; i < a.size(); ++i) {
		os << *a[i];
	}
	return os;
}

/*
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
*/

// fast stabbing
//template <typename interval>
class faststabbing
{
private:
    std::vector<interval>& a; // array of intervals [0,n-1]
	uint64_t n,bigN;
    std::vector<std::vector<interval*> > eventlist; // sweepline
    std::vector<interval*> stop;
	interval dummy;

    void preprocessing(void) {
        // create smaller lists and event lists
        uint64_t i,l,starting=-1;
        for (i=0; i<n; ++i) {
            l = a[i].l;
            if (l != starting) {
                // sorted event lists for sweepline
                eventlist[a[i].r].push_back(&a[i]);
                eventlist[l].push_back(&a[i]);
            } else {
                assert(a[i-1].l == l && a[i-1].r > a[i].r);
                a[i-1].smaller = &a[i];
            }
            starting = l;
        }

        // sweep line
        std::list<interval*> L; // status list
        interval* temp;
        interval* last;
        for (i=1; i<=bigN; ++i) {
            // interval with starting point i
            if (!eventlist[i].empty()) {
                temp = eventlist[i].back();
                if (temp->l == i) {
                    L.push_back(temp);
                    temp->pIt = std::prev(L.end());
                    eventlist[i].pop_back();
                }
            }
            //log << "\n" << i << ": " << eventlist[i] << ", \n\tL=" << L << "\n";
            assert(!L.empty() || eventlist[i].empty());
            if (!L.empty()) {
                // compute stop[i]
                stop[i] = L.back();
                // intervals with end points i
                for (auto it = eventlist[i].rbegin(); it != eventlist[i].rend(); --it) {
                    temp = *it;
                    if (temp->pIt != L.begin()) {
                        last = *std::prev(temp->pIt);
                    } else last = &dummy;
                    //log << "\n\t\t" << last << "\t\t" << temp;
                    temp->parent = last;
                    temp->leftsibling = last->rightchild;
                    last->rightchild = temp;
                    L.erase(temp->pIt);
                    last = temp;
                }
            }
        }
//#ifdef INTERVALSTAB_DEBUG
//	cout << "\nDummy\t\t" << &dummy << "\n" << a << "\n";
//#endi
    }

//#ifdef INTERVALSTAB_DEBUG
    bool verify(std::vector<interval*> output, const uint64_t& q) {
//	cout << "\nQuery q=" << q << ":\n" << output;
        interval* temp;
        interval* last = nullptr;
        while (!output.empty()) {
            temp = output.back();
            output.pop_back();
            if (last && (*temp < *last)) {
                std::cerr << "\nerror: interval " << temp << " not in order (not after " << last << ")\n";
                return 1;
            }
            temp->stabbed = true;
            last = temp;
        }
        bool stabs;
        for (uint64_t i=0; i<n; ++i) {
            stabs = a[i].l <= q && q <= a[i].r;
            if (a[i].stabbed != stabs) {
                std::cerr << "\nerror: interval " << i << " (" << &a[i] << ") should be" << (stabs ? " stabbed\n" : " not stabbed\n");
                return 1;
            }
            a[i].stabbed = false;
        }
        return 0;
    }
//#endif

public:
	faststabbing(std::vector<interval>& intervals,
                 const uint64_t& numberIntervals,
                 const uint64_t& numberDomain)
        : a(intervals), eventlist(numberDomain+1), stop(numberDomain)  {
		n = numberIntervals;
		bigN = numberDomain;
		dummy.parent = nullptr;
		dummy.leftsibling = nullptr;
		dummy.rightchild = nullptr;
		preprocessing();
	};

	std::vector<interval*> query(const uint64_t& q) { //, uint64_t& numComparisons) {
        assert(q >= 1 && q <= bigN+1);
        std::vector<interval*> output;
        if (stop[q] == nullptr) return output; // no stabbed intervals
        interval* i;
        interval* temp;
        std::deque<interval*> process;
        for (temp = stop[q]; temp->parent != nullptr; temp = temp->parent) {
            process.push_front(temp);
        }

        // traverse
        while (!process.empty()) {
            i = process.back();
            process.pop_back();
            //process.pop_back();
            output.push_back(i);
		
            temp = i->smaller;
            while (temp != nullptr) {
                //++numComparisons;
                if (q > temp->r) break;
                output.push_back(temp);
//#ifdef INTERVALSTAB_DEBUG
//			cout << "\tSmaller " << (*temp);
//#endif
                temp = temp->smaller;
            }

            // go along rightmost path of pa
            temp = i->leftsibling;
            while (temp) {
                //++numComparisons;
                if (temp->r < q) break;
                process.push_back(temp);
                temp = temp->rightchild;
            }
        }
        assert(verify(output,q) == 0);
        return output;
    }
};


}
