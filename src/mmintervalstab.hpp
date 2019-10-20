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
#include <limits>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <fcntl.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <cassert>
#include "ips4o.hpp"
#include "mmappable_vector.h"

namespace intervalstab {

using namespace mmap_allocator_namespace;

// an interval
template <typename T>
struct interval {
	uint64_t l = 0;
    uint64_t r = 0;
	interval<T>* leftsibling = nullptr;
	interval<T>* rightchild = nullptr;
	interval<T>* parent = nullptr;
	interval<T>* smaller = nullptr;
    typename std::list<interval<T>*>::iterator pIt; // = (std::list<interval*>::iterator)nullptr;
    bool stabbed = false;
    T value;
    interval<T>(void) { }
    interval<T>(const uint64_t& a,
                const uint64_t& b,
                const T& d) : l(a), r(b), value(d) { }
    ~interval<T>(void) { }
};

// lexicographic order
template <typename T>
inline bool operator<(const interval<T>& x,const interval<T>& y) {
	return (x.l < y.l || x.l == y.l && x.r > y.r);
}

// equality
template <typename T>
inline bool operator==(const interval<T>& x,const interval<T>& y) {
	return (x.l == y.l && x.r == y.r);
}

// lexicographic order
template <typename T>
inline bool operator>(const interval<T>& x,const interval<T>& y) {
	return y < x;
}

// output stream for intervals
template <typename T>
inline std::ostream& operator<<(std::ostream& os, const interval<T>& a) {
	os << &a << "\t" << a.l << "\t" << a.r << "\tP " << a.parent << " L " << a.leftsibling
       << " C " << a.rightchild << "  Sm " << a.smaller;
	return os;
}
template <typename T>
inline std::ostream& operator<<(std::ostream& os, const std::vector<interval<T>*>& a) {
	for (unsigned int i = 0; i < a.size(); ++i) {
		os << *a[i];
	}
	return os;
}

template <typename T>
void fill_file(const char *fname, const uint64_t& count) {
    std::ofstream out(fname, std::ios_base::binary | std::ios::trunc);
    T x;
    for (uint64_t i=0; i<count; i++) {
        out.write((char*)&x, sizeof(T));
    }
    out.close();
}
    
// fast stabbing
//template <typename interval> // TODO
template <typename T>
class faststabbing
{
private:

    int get_thread_count(void) {
        int thread_count = 1;
#pragma omp parallel
        {
#pragma omp master
            thread_count = omp_get_num_threads();
        }
        return thread_count;
    }

    std::ofstream writer;
    std::vector<std::ofstream> writers;
    char* reader = nullptr;
    int reader_fd = 0;
    std::string filename;
    std::string index_filename;
    bool sorted = false;
    // key information
    uint64_t n_records = 0;
    bool indexed = false;
    uint32_t OUTPUT_VERSION = 1; // update as we change our format

public:
    //~iit_mm_builder_base(void) { close_writers(); }

    void set_base_filename(const std::string& f) {
        filename = f;
    }

    // close/open backing file
    void open_main_writer(void) {
        if (writer.is_open()) {
            writer.seekp(0, std::ios_base::end); // seek to the end for appending
            return;
        }
        assert(!filename.empty());
        writer.open(filename.c_str(), std::ios::binary | std::ios::trunc);
        if (writer.fail()) {
            throw std::ios_base::failure(std::strerror(errno));
        }
    }

    // per-thread writers
    void open_writers(const std::string& f) {
        set_base_filename(f);
        open_writers();
    }

    void open_writers(void) {
        assert(!filename.empty());
        writers.clear();
        writers.resize(get_thread_count());
        for (size_t i = 0; i < writers.size(); ++i) {
            auto& writer = writers[i];
            writer.open(writer_filename(i), std::ios::binary | std::ios::app);
            if (writer.fail()) {
                throw std::ios_base::failure(std::strerror(errno));
            }
        }
    }

    std::string writer_filename(size_t i) {
        std::stringstream wf;
        wf << filename << ".tmp_write" << "." << i;
        return wf.str();
    }

    std::string eventlist_filename(void) {
        return filename + ".eventlist";
    }

    std::string stop_filename(void) {
        return filename + ".stop";
    }

    std::ofstream& get_writer(void) {
        return writers[omp_get_thread_num()];
    }

    void sync_and_close_parallel_writers(void) {
        // close the temp writers and cat them onto the end of the main file
        open_main_writer();
        for (size_t i = 0; i < writers.size(); ++i) {
            writers[i].close();
            std::ifstream if_w(writer_filename(i), std::ios_base::binary);
            writer << if_w.rdbuf();
            if_w.close();
            std::remove(writer_filename(i).c_str());
        }
        writers.clear();
        writer.close();
        for (size_t i = 0; i < writers.size(); ++i) {
            std::remove(writer_filename(i).c_str());
        }
    }

    /// return the number of records, which will only work after indexing
    size_t size(void) const {
        return n_records;
    }

    /// return the backing buffer
    char* get_buffer(void) const {
        return reader;
    }

    /// get the record count
    size_t record_count(void) {
        int fd = open(filename.c_str(), O_RDWR);
        if (fd == -1) {
            assert(false);
        }
        struct stat stats;
        if (-1 == fstat(fd, &stats)) {
            assert(false);
        }
        assert(stats.st_size % sizeof(interval<T>) == 0); // must be even records
        size_t count = stats.st_size / sizeof(interval<T>);
        close(fd);
        return count;
    }

    mmappable_vector<interval<T>> a; // array of intervals [0,n-1]
	uint64_t n,bigN;
    mmappable_vector<std::vector<interval<T>*> > eventlist; // sweepline
    mmappable_vector<interval<T>*> stop;
	interval<T> dummy;

    void preprocessing(void) {
        // sort the array

        // if we want to work on unique data
        //a.erase(std::unique(a.begin(), a.end()), a.end());
        //((Interval*)buffer.data)+data_len,
        //IntervalLess());
        // create smaller lists and event lists
        uint64_t i,l,starting=-1;
        for (i=0; i<n; ++i) {
            l = a[i].l;
            //std::cerr << "processing " << a[i] << std::endl;
            if (l != starting) {
                // sorted event lists for sweepline
                eventlist[a[i].r].push_back(&a[i]);
                eventlist[l].push_back(&a[i]);
            } else {
                assert(a[i-1].l == l && a[i-1].r >= a[i].r);
                a[i-1].smaller = &a[i];
            }
            starting = l;
        }

        // sweep line
        std::list<interval<T>*> L; // status list
        interval<T>* temp;
        interval<T>* last;
        for (i=1; i<=bigN; ++i) {
            // interval with starting point i
            if (!eventlist[i].empty()) {
                //std::cerr << "eventlist size " << eventlist[i].size() << std::endl;
                temp = eventlist[i].back();
                if (temp->l == i) {
                    L.push_back(temp);
                    temp->pIt = std::prev(L.end());
                    eventlist[i].pop_back();
                }
            }
            /*
            std::cerr << "sweeep " << i << ": " << eventlist[i];
            for (auto& l : L) std::cerr << " " << l;
            std::cerr << std::endl;
            */
            //assert(!L.empty() || eventlist[i].empty());
            if (!L.empty()) {
                // compute stop[i]
                stop[i] = L.back();
                // intervals with end points i
                for (auto it = eventlist[i].rbegin(); it != eventlist[i].rend(); ++it) {
                    temp = *it;
                    //std::cerr << "Temp " << temp->l << " " << temp->r << std::endl;
                    if (temp->pIt != L.begin()) {
                        //std::cerr << "setting last " << *temp << std::endl;
                        last = *std::prev(temp->pIt);
                    } else last = &dummy;
                    //std::cerr << "\n\t\t" << last << "\t\t" << temp << std::endl;
                    temp->parent = last;
                    temp->leftsibling = last->rightchild;
                    last->rightchild = temp;
                    //temp->pIt =
                    //std::cerr << "L size " << L.size() << std::endl;
                    //if (temp->pIt != L.end())
                    L.erase(temp->pIt);
                    //temp->pIt = std::prev(L.end());
                    last = temp;
                }
            }
        }
//#ifdef INTERVALSTAB_DEBUG
        //std::cerr << "\nDummy\t\t" << &dummy << "\n" << a.size() << std::endl;
//#endi
    }

//#ifdef INTERVALSTAB_DEBUG
    bool verify(mmappable_vector<interval<T>*> output, const uint64_t& q) {
//	cout << "\nQuery q=" << q << ":\n" << output;
        interval<T>* temp;
        interval<T>* last = nullptr;
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
	faststabbing(const std::string& f)
        : filename(f) {
		dummy.parent = nullptr;
		dummy.leftsibling = nullptr;
		dummy.rightchild = nullptr;
        open_writers(f);
	};

    void add(const interval<T>& it) {
        auto& writer = get_writer();
        writer.write((char*)&it, sizeof(interval<T>));
    }

    void index(void) {
        // calculate numberDomain, numberIntervals, n, and bigN
        // sync the writers and mmap the file into our vector
        sync_and_close_parallel_writers();
        a.mmap_file(filename.c_str(), READ_WRITE_SHARED, 0, record_count());
        n = record_count(); // number of intervals
        ips4o::parallel::sort(a.begin(), a.end()); // sort the intervals
        uint64_t domain_count = 0; // find the domain of our integer space
        for (auto& i : a) { if (i.r > domain_count) domain_count = i.r; }
        bigN = domain_count; // number of domains
        // mmap our sweepline and stop
        fill_file<std::vector<interval<T>*>>(eventlist_filename().c_str(), bigN+1);
        fill_file<interval<T>*>(stop_filename().c_str(), bigN+1);
        eventlist.mmap_file(eventlist_filename().c_str(), READ_WRITE_SHARED, 0, bigN+1);
        stop.mmap_file(stop_filename().c_str(), READ_WRITE_SHARED, 0, bigN+1);
        for (auto& i : stop) { i = nullptr; }
        // build our data structures
        preprocessing();
    }

    std::vector<interval<T>*> query(const uint64_t& q) {
        assert(q >= 1 && q <= bigN+1);
        std::vector<interval<T>*> output;
        if (stop[q] == nullptr) return output; // no stabbed intervals
        interval<T>* i;
        interval<T>* temp;
        std::deque<interval<T>*> process;
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
                if (temp->r < q) break;
                process.push_back(temp);
                temp = temp->rightchild;
            }
        }
        //assert(verify(output,q) == 0);
        return output;
    }
};


}
