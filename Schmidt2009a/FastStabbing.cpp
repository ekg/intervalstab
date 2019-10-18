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

#include "FastStabbing.h"

namespace ogdf {

// verifies correctness of output, 0=correct, 1=not correct
#ifdef OGDF_DEBUG
bool FastStabbing::verify(std::vector<interval*> output, const int& q) {
//	cout << "\nQuery q=" << q << ":\n" << output;
	interval* temp;
	interval* last = NULL;
	while (!output.empty()) {
		temp = output.back();
		output.pop_back();
		if (last && (*temp < *last)) {
			cerr << "\nerror: interval " << temp << " not in order (not after " << last << ")\n";
			return 1;
		}
		temp->stabbed = true;
		last = temp;
	}
	bool stabs;
	for (int i=0; i<n; ++i) {
		stabs = a[i].l <= q && q <= a[i].r;
		if (a[i].stabbed != stabs) {
			cerr << "\nerror: interval " << i << " (" << &a[i] << ") should be" << (stabs ? " stabbed\n" : " not stabbed\n");
			return 1;
		}
		a[i].stabbed = false;
	}
	return 0;
}
bool ChazelleStabbing::verify(std::vector<interval*> output, const int& q) { // some different output information here...
//	cout << "\nQuery q=" << q << " (window start=" << (*pWindow[q]).l << "):\n" << output;
	OGDF_ASSERT(delta*output.size() >= (*pWindow[q]).intervals.size());
	interval* temp;
	interval* last = NULL;
	while (!output.empty()) {
		temp = output.back();
		output.pop_back();
		if (last && (*temp > *last)) {
			cerr << "\n" << a;
			cerr << "\nQuery q=" << q << " (window start=" << (*pWindow[q]).l << "):\n" << (*pWindow[q]).intervals << "\n" << output;
			cerr << "\nerror: interval " << temp << " not in order (not after " << last << ")\n";
			return 1;
		}
		temp->stabbed = true;
		last = temp;
	}
	bool stabs;
	for (int i=0; i<n; ++i) {
		stabs = a[i].l <= q && q <= a[i].r;
		if (a[i].stabbed != stabs) {
			cerr << "\n" << a;
			cerr << "\nQuery q=" << q << " (window start=" << (*pWindow[q]).l << "):\n" << (*pWindow[q]).intervals << "\n" << output;
			cerr << "\nerror: interval " << i << " (" << &a[i] << ") should be" << (stabs ? " stabbed\n" : " not stabbed\n");
			return 1;
		}
		a[i].stabbed = false;
	}
	return 0;
}
#endif

void FastStabbing::preprocessing() {
	// create smaller lists and event lists
	int64_t i,l,starting=-1;
	for (i=0; i<n; ++i) {
		l = a[i].l;
		if (l != starting) {
			// sorted event lists for sweepline
			eventlist[a[i].r].pushBack(&a[i]);
			eventlist[l].pushBack(&a[i]);
		} else {
			OGDF_ASSERT(a[i-1].l == l && a[i-1].r > a[i].r);
			a[i-1].smaller = &a[i];
		}
		starting = l;
	}

	// sweep line
	ListPure<interval*> L; // status list
	ListIterator<interval*> it;
	interval* temp;
	interval* last;
	for (i=1; i<=bigN; ++i) {
		// interval with starting point i
		if (!eventlist[i].empty()) {
			temp = eventlist[i].back();
			if (temp->l == i) {
				temp->pIt = L.pushBack(temp);
				eventlist[i].popBack();
			}
		}
		//log << "\n" << i << ": " << eventlist[i] << ", \n\tL=" << L << "\n";
		OGDF_ASSERT(!L.empty() || eventlist[i].empty());
		if (!L.empty()) {
			// compute stop[i]
			stop[i] = L.back();
			// intervals with end points i
			for (it = eventlist[i].rbegin(); it.valid(); --it) {
				temp = *it;
				if (temp->pIt.pred().valid()) {
					last = *temp->pIt.pred();
				} else last = &dummy;
				//log << "\n\t\t" << last << "\t\t" << temp;
				temp->parent = last;
				temp->leftsibling = last->rightchild;
				last->rightchild = temp;
				L.del(temp->pIt);
				last = temp;
			}
		}
	}
	#ifdef OGDF_DEBUG
//	cout << "\nDummy\t\t" << &dummy << "\n" << a << "\n";
	#endif
}

// stabbing query
void FastStabbing::query(const int& q, std::vector<interval*>& output, const bool onlySearch, long& numComparisons) {
	OGDF_ASSERT(q >= 1 && q <= bigN+1);
	output.clear();
	if (stop[q] == NULL) return; // no stabbed intervals
	if (onlySearch) {
		output.push_back(stop[q]);
		return;
	}
	interval* i;
	interval* temp;
	ListPure<interval*> process;
	for (temp = stop[q]; temp->parent != NULL; temp = temp->parent) {
		process.pushFront(temp);
	}

	// traverse
	while (!process.empty()) {
		i = process.popBackRet();
		//process.pop_back();
		output.push_back(i);
		
		temp = i->smaller;
		while (temp != NULL) {
			++numComparisons;
			if (q > temp->r) break;
			output.push_back(temp);
			#ifdef OGDF_DEBUG
//			cout << "\tSmaller " << (*temp);
			#endif
			temp = temp->smaller;
		}

		// go along rightmost path of pa
		temp = i->leftsibling;
		while (temp) {
			++numComparisons;
			if (temp->r < q) break;
			process.pushBack(temp);
			temp = temp->rightchild;
		}
	}
	OGDF_ASSERT(verify(output,q) == 0);
}

// make windows
void ChazelleStabbing::preprocessing() {
	// create smaller lists and event lists
	int i;
	for (i=0; i<n; ++i) {
		// sorted event lists for sweepline
		if (a[i].l == a[i].r) {
			eventlist[a[i].l].pushBack(&a[i]);
		} else {
			eventlist[a[i].r].pushBack(&a[i]);
			eventlist[a[i].l].pushBack(&a[i]);
		}
	}

	// sweep line
	interval* temp;
	int cur = 0, low = 0, T = 0;
	window dummy; // current window
	dummy.l = -1; // avoid null-aperture when first interval starts at value 1
	windows.pushBack(dummy);
	window* w = &windows.back();
	ListIterator<interval*> events,lastevents,it,del;
	SListIterator<window> lastpWindow = windows.end(); // marker for coordinates with no intervals stabbed
	for (i=1; i<=bigN; ++i) {
		//pWindow[i] = windows.rbegin();
		pWindow[i] = lastpWindow;
//		cout << "\n" << i;
		// interval with starting point i
		for (events=eventlist[i].begin(); events.valid();) {
			temp = *events;
			lastevents = events;
			++events;
			if (temp->l == i) {
				if (temp->r == i) {
					continue; // only degenerated intervals remain in eventlist and are not processed
				} else eventlist[i].del(lastevents);
				// starting point
				++cur; // number of active intervals
				++T; // number of intervals in current window
				//cout << "\tT=" << T << ",delta*low=" << delta*low;
				if (T > delta*low) { // full window
//					cout << "\nl-window (start=" << w->l << ",end=" << i << "): " << w->intervals;
					OGDF_ASSERT(w->intervals.size() <= T);
					// aperture of null: modify null-window to a bigger one
					// aperture non-null (and not the first interval after an empty window): create new window
					if (w->l < i && T > 1) {
						windows.pushBack(*w);
						w = &windows.back();
						// delete intervals that are not stabbed
						for (it = w->intervals.begin(); it.valid();) {
							if ((*it)->r <= i) {
								del = it;
								++it;
								w->intervals.del(del);
							} else ++it;
						}
					}
					// initialize new window
					w->l = i;
					w->intervals.pushBack(temp);
					pWindow[i] = windows.rbegin();
					lastpWindow = windows.rbegin();
//					cout << "\nadded interval: " << w->intervals;
					low = cur;
					T = cur;
				} else {
					w->intervals.pushBack(temp); // insert interval in current window
//					cout << "\nadded interval: " << w->intervals;
				}
			} else {
				eventlist[i].del(lastevents);
				// end point
				--cur;
				if (cur < low) low = cur;
				if (T > delta*low) { // full window
//					cout << "\nr-window (start=" << w->l << ",end=" << i << "): " << w->intervals;
					OGDF_ASSERT(w->intervals.size() <= T);
					// aperture of null: modify null-window to a bigger one
					// aperture non-null: create new window
					if (w->l < i) {
						windows.pushBack(*w);
						w = &windows.back();
						// delete intervals that are not stabbed
						for (it = w->intervals.begin(); it.valid();) {
							if ((*it)->r <= i) {
								del = it;
								++it;
								w->intervals.del(del);
								--T;
							} else ++it;
						}
//						cout << "\nafter deletion: " << w->intervals << " (start=" << i << ")";
					} else T = cur;
					// initialize new window
					w->l = i;
					//T = cur;
					low = T;
					if (T==0) {
						lastpWindow = windows.end();
					} else {
						pWindow[i] = windows.rbegin();
						lastpWindow = windows.rbegin();
					}
				} else {
//					cout << "\nremoved interval";
				}
			}
		}
	}
	#ifdef OGDF_DEBUG
//	cout << "\n" << a << "\npWindow(" << &(*windows.end()) << "):\n";
	ListIterator<window*> itW;
//	for (i=1; i<=bigN; ++i) {
//		cout << i << ": ";
//		if (!pWindow[i].valid()) {
//			cout << "NULL\n";
//		} else cout << (*pWindow[i]).l << "\n";
//	}
//	cout << "\n";
//	cout << "\neventlist:\n" << eventlist;
	int count = 0;
	for (SListIterator<window> itW = windows.begin(); itW.valid(); ++itW) {
		OGDF_ASSERT(itW == windows.rbegin() || !(*itW).intervals.empty());
		count = count + (*itW).intervals.size();
	}
	OGDF_ASSERT(count < 2*delta*n/(delta-1));
	#endif
}

// stabbing query chazelle
void ChazelleStabbing::query(const int& q, std::vector<interval*>& output, const bool onlySearch, long& numComparisons) {
	OGDF_ASSERT(q >= 1 && q <= bigN+1);
	output.clear();
	// compute degenerated interval, if exists
	interval* temp;

	if (onlySearch) {
		// ..., feature not used.
	}

	if (!eventlist[q].empty()) {
		OGDF_ASSERT(eventlist[q].size()==1);
		temp = eventlist[q].front();
//		cout << "\ndegenerated interval " << temp << " on " << q;
	} else temp = NULL;
	// no stabbed interval
	if (!pWindow[q].valid()) {
		if (temp != NULL) output.push_back(temp);
		return;
	}
	// search in preceding window as well, if q hits intersection
	if (q > 1 && (*pWindow[q]).l == q && pWindow[q-1].valid()) {
//		cout << "\n\tPreceding window:\n\t" << (*pWindow[q-1]).intervals;
		for (ListIterator<interval*> it = (*pWindow[q-1]).intervals.begin(); it.valid(); ++it) {
			++numComparisons;
			if ((*it)->l <= q) {
				++numComparisons;
				if (q <= (*it)->r) {
					if (temp != NULL && temp->l <= (*it)->l) {
						output.push_back(temp);
						temp = NULL;
						if (onlySearch) return;
					}
					output.push_back(*it);
					if (onlySearch) return;
					(*it)->parent = *it; // make parent nonnull as a marker
				}
			}
		}
	}
	// search in linked window
	for (ListIterator<interval*> it = (*pWindow[q]).intervals.begin(); it.valid(); ++it) {
			++numComparisons;
			if ((*it)->l <= q) {
				++numComparisons;
				if (q <= (*it)->r && (*it)->parent == NULL) {
					if (temp != NULL && temp->l <= (*it)->l) {
						output.push_back(temp);
						temp = NULL;
					}
					output.push_back(*it);
					if (onlySearch) return;
				}
			}
	}
	if (temp != NULL)
		output.push_back(temp);
	for (unsigned int i = 0; i < output.size(); ++i)
		output[i]->parent = NULL;
	OGDF_ASSERT(verify(output,q) == 0);
}

}
