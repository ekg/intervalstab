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
    T value;
    interval<T>(void) { }
    interval<T>(const uint64_t& a,
                const uint64_t& b,
                const T& d) : l(a), r(b), value(d) { }
    ~interval<T>(void) { }
};

template <typename T>
struct interval_node : interval<T> {
	interval_node<T>* leftsibling = nullptr;
	interval_node<T>* rightchild = nullptr;
	interval_node<T>* parent = nullptr;
	interval_node<T>* smaller = nullptr;
    typename std::list<interval_node<T>*>::iterator pIt;
    bool stabbed = false;
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
inline std::ostream& operator<<(std::ostream& os, const interval_node<T>& a) {
	os << &a << "\t" << a.l << "\t" << a.r << "\tP " << a.parent << " L " << a.leftsibling
       << " C " << a.rightchild << "  Sm " << a.smaller;
	return os;
}
template <typename T>
inline std::ostream& operator<<(std::ostream& os, const std::vector<interval_node<T>*>& a) {
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
    bool sorted = false;
    // key information
    uint64_t n_records = 0;
    bool indexed = false;
    uint32_t OUTPUT_VERSION = 1; // update as we change our format

public:

    void set_base_filename(const std::string& f) {
        filename = f;
    }

    // close/open backing file
    void open_main_writer(void) {
        if (writer.is_open()) {
            writer.seekp(0, std::ios_base::end); // seek to the end for appending
            return;
        }
        assert(!intervals_filename().empty());
        writer.open(intervals_filename().c_str(), std::ios::binary | std::ios::trunc);
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

    std::string intervals_filename(void) {
        return filename + ".intervals";
    }

    std::string eventlist_filename(void) {
        return filename + ".eventlist";
    }

    std::string eventlist_layout_filename(void) {
        return filename + ".eventlist.layout";
    }

    std::string node_filename(void) {
        return filename + ".nodes";
    }

    std::string stop_filename(void) {
        return filename + ".stop";
    }

    std::ofstream& get_writer(void) {
        return writers[omp_get_thread_num()];
    }

    void sync_and_close_parallel_writers(void) {
        // check to see if we ran single-threaded
        uint64_t used_writers = 0;
        uint64_t writer_that_wrote = 0;
        for (size_t i = 0; i < writers.size(); ++i) {
            writers[i].close();
            if (filesize(writer_filename(i).c_str())) {
                ++used_writers;
                writer_that_wrote = i;
            }
        }
        bool single_threaded = used_writers == 1;
        // close the temp writers and cat them onto the end of the main file
        if (single_threaded) {
            std::rename(writer_filename(writer_that_wrote).c_str(), intervals_filename().c_str());
        } else {
            open_main_writer();
            for (size_t i = 0; i < writers.size(); ++i) {
                std::ifstream if_w(writer_filename(i), std::ios_base::binary);
                writer << if_w.rdbuf();
                if_w.close();
                std::remove(writer_filename(i).c_str());
            }
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
        int fd = open(intervals_filename().c_str(), O_RDWR);
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

    std::ifstream::pos_type filesize(const char* filename) {
        std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
        return in.tellg(); 
    }

    mmappable_vector<interval<T>> intervals; // array of intervals [0,n-1]
    mmappable_vector<interval_node<T>> a; // array of interval contexts [0,n-1]
	uint64_t n,bigN;
    mmappable_vector<interval_node<T>*> eventlist;
    mmappable_vector<uint64_t> eventlist_layout;
    //lciv_iv eventlist;
    //suc_bv eventlist_delim;

    mmappable_vector<interval_node<T>*> stop;
	interval_node<T> dummy;

    void preprocessing(void) {
        // calculate numberDomain, numberIntervals, n, and bigN
        // sync the writers and mmap the file into our vector
        sync_and_close_parallel_writers();
        intervals.mmap_file(intervals_filename().c_str(), READ_WRITE_SHARED, 0, record_count());
        n = record_count(); // number of intervals
        ips4o::parallel::sort(intervals.begin(), intervals.end()); // sort the intervals
        uint64_t domain_count = 0; // find the domain of our integer space
        for (auto& i : intervals) { if (i.r > domain_count) domain_count = i.r; }
        bigN = domain_count; // number of domains
        //std::cerr << "bigN = " << bigN << std::endl;
        fill_file<interval_node<T>>(node_filename().c_str(), n);
        a.mmap_file(node_filename().c_str(), READ_WRITE_SHARED, 0, n);
        // copy intervals into our stabbing tree
        for (uint64_t i = 0; i < n; ++i) {
            auto& o = intervals[i];
            a[i].l = o.l;
            a[i].r = o.r;
            a[i].value = o.value;
        }
        // clean up intervals file
        intervals.munmap_file();
        std::remove(intervals_filename().c_str());
        // mmap our sweepline and stop
        fill_file<interval_node<T>*>(stop_filename().c_str(), bigN+1);
        stop.mmap_file(stop_filename().c_str(), READ_WRITE_SHARED, 0, bigN+1);
        for (auto& i : stop) { i = nullptr; }
        // mmap our eventlist
        uint64_t eventlist_size = 0;
        uint64_t l=0,starting=-1;
        for (uint64_t i=0; i<n; ++i) {
            l = a[i].l;
            if (l != starting) {
                eventlist_size += 2;
            }
            starting = l;
        }

        fill_file<interval_node<T>*>(eventlist_filename().c_str(), eventlist_size);
        eventlist.mmap_file(eventlist_filename().c_str(), READ_WRITE_SHARED, 0, eventlist_size);
        for (uint64_t i=0; i<eventlist_size; ++i) {
            eventlist[i] = nullptr;
        }

        fill_file<uint64_t>(eventlist_layout_filename().c_str(), bigN+2);
        eventlist_layout.mmap_file(eventlist_layout_filename().c_str(), READ_WRITE_SHARED, 0, bigN+2);
        for (uint64_t i=0; i<bigN+2; ++i) {
            eventlist_layout[i] = 0;
        }

        // determine the layout, using our eventlist_layout to temporarily store the counts
        l=0; starting=-1;
        for (uint64_t i=0; i<n; ++i) {
            l = a[i].l;
            if (l != starting) {
                ++eventlist_layout[a[i].r];
                ++eventlist_layout[l];
            } else {
                assert(a[i-1].l == l && a[i-1].r >= a[i].r);
                a[i-1].smaller = &a[i];
            }
            starting = l;
        }
        // record the layout offsets in the eventlist_layout
        uint64_t offset = 0;
        for (uint64_t i=1; i<=bigN+1; ++i) {
            uint64_t count = eventlist_layout[i];
            //std::cerr << "count at " << i << " = " << count << std::endl;
            eventlist_layout[i] = offset;
            offset += count;
            //std::cerr << "eventlist_layout " << i << " -> " << offset << std::endl;
        }
        std::cerr << "eventlist size " << eventlist.size() << std::endl;
        // build our data structures
        l=0; starting=-1;
        for (uint64_t i=0; i<n; ++i) {
            if (i % 1000 == 0) {
                std::cerr << "eventlist " << i << "\r";
            }
            l = a[i].l;
            //std::cerr << "processing " << a[i] << std::endl;
            if (l != starting) {
                // sorted event lists for sweepline
                // find were to insert
                uint64_t write_at = eventlist_layout[a[i].r];
                //std::cerr << "write at for r " << write_at << std::endl;
                while (eventlist[write_at] != nullptr) ++write_at;
                eventlist[write_at] = &a[i];
                write_at = eventlist_layout[l];
                //std::cerr << "write at for l " << write_at << std::endl;
                while (eventlist[write_at] != nullptr) ++write_at;
                eventlist[write_at] = &a[i];
            } else {
                assert(a[i-1].l == l && a[i-1].r >= a[i].r);
                a[i-1].smaller = &a[i];
            }
            starting = l;
        }
        std::cerr << std::endl;

        // sweep line
        std::list<interval_node<T>*> L; // status list
        interval_node<T>* temp;
        interval_node<T>* last;
        for (uint64_t i=1; i<=bigN; ++i) {
            //for (uint64_t i=1; i<=bigN; ++i) {
            if (i % 1000 == 0) {
                std::cerr << "building " << i << "\r";
            }
            // interval with starting point i
            //uint64_t x = eventlist_delim.select1(i-1)+1;
            //uint64_t y = eventlist_delim.select1(i);
            uint64_t x = eventlist_layout[i];
            uint64_t y = eventlist_layout[i+1];
            if (y - x > 0) {
                //std::cerr << "eventlist size " << y - x << std::endl;
                //temp = eventlist[i].back();
                //temp = &a[eventlist.at(x)];
                uint64_t read_at = y-1;
                while (read_at != x-1 && eventlist[read_at] == nullptr) --read_at;
                if (read_at != x-1) {
                    temp = eventlist[read_at];
                    if (temp->l == i) {
                        L.push_back(temp);
                        temp->pIt = std::prev(L.end());
                        eventlist[read_at] = nullptr;
                    }
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
                uint64_t x = eventlist_layout[i];
                uint64_t y = eventlist_layout[i+1];
                if (y - x > 0) {
                    uint64_t read_at = y-1;
                    while (read_at != x-1 && eventlist[read_at] == nullptr) --read_at;
                    //if (read_at == x-1) last = &dummy;
                    //std::cerr << "read_at = " << read_at << std::endl;
                    for (uint64_t j = read_at; j != x-1; --j) {
                        //std::cerr << "looking at eventlist " << j << std::endl;
                        temp = eventlist[j];
                        //std::cerr << "temp " << temp << std::endl;
                        //std::cerr << "Temp " << temp->l << " " << temp->r << std::endl;
                        if (temp->pIt != L.begin()) {
                            //std::cerr << "setting last " << *temp << std::endl;
                            last = *std::prev(temp->pIt);
                        } else last = &dummy;
                        //std::cerr << "\n\t\t" << last << "\t\t" << temp << std::endl;
                        temp->parent = last;
                        temp->leftsibling = last->rightchild;
                        last->rightchild = temp;
                        //std::cerr << "L size " << L.size() << std::endl;
                        L.erase(temp->pIt);
                        last = temp;
                    }
                }
            }
        }
        std::cerr << std::endl;

        eventlist.munmap_file();
        std::remove(eventlist_filename().c_str());
        eventlist_layout.munmap_file();
        std::remove(eventlist_layout_filename().c_str());
//#ifdef INTERVALSTAB_DEBUG
        //std::cerr << "\nDummy\t\t" << &dummy << "\n" << a.size() << std::endl;
//#endi
    }

//#ifdef INTERVALSTAB_DEBUG
    bool verify(mmappable_vector<interval_node<T>*> output, const uint64_t& q) {
//	cout << "\nQuery q=" << q << ":\n" << output;
        interval_node<T>* temp;
        interval_node<T>* last = nullptr;
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

    ~faststabbing(void) {
        a.munmap_file();
        std::remove(node_filename().c_str());
        stop.munmap_file();
        std::remove(stop_filename().c_str());
    }

    void add(const interval<T>& it) {
        auto& writer = get_writer();
        writer.write((char*)&it, sizeof(interval<T>));
    }

    void index(void) {
        preprocessing();
    }

    std::vector<interval_node<T>*> query(const uint64_t& q) {
        assert(q >= 1 && q <= bigN+1);
        std::vector<interval_node<T>*> output;
        if (stop[q] == nullptr) return output; // no stabbed intervals
        interval_node<T>* i;
        interval_node<T>* temp;
        std::deque<interval_node<T>*> process;
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
