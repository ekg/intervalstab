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


#include <iostream>
#include <vector>
#include <random>
#include "intervalstab.hpp"
#include "args.hxx"

using namespace intervalstab;

int main(int argc, char** argv) {

    args::ArgumentParser parser("memmapped interpolated implicit interval tree");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> test_file(parser, "FILE", "test mmmultimap with random data in this file", {'T', "test-file"});
    args::ValueFlag<uint64_t> test_size(parser, "N", "test this many pairs", {'s', "test-size"});
    args::ValueFlag<uint64_t> max_val(parser, "N", "generate test data in the range [1,max_value]", {'M', "max-value"});
    args::ValueFlag<uint64_t> range_mean(parser, "N", "the mean length for intervals (under gaussian distribution)", {'m', "range-mean"});
    args::ValueFlag<uint64_t> range_stdev(parser, "N", "the standard deviation for intervals (under gaussian distribution)", {'D', "range-stdev"});
    args::ValueFlag<uint64_t> threads(parser, "N", "number of threads to use", {'t', "threads"});
    args::ValueFlag<uint64_t> domains(parser, "N", "number of domains for interpolation", {'d', "domains"});

    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    if (argc==1) {
        std::cout << parser;
        return 1;
    }
    
    assert(!args::get(test_file).empty());
    assert(args::get(test_size));
    assert(args::get(max_val));


    std::vector<interval> intervals;
    
    //std::remove(args::get(test_file).c_str());
    //p_iitii::builder bb = p_iitii::builder(args::get(test_file));
    //p_iitii::builder bb = p_iitii::builder(args::get(test_file));

    //bb.add(intpair(12,34));
    //bb.add(intpair(0,23));
    //bb.add(intpair(34,56));

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    uint64_t max_value = args::get(max_val);
    std::uniform_int_distribution<uint64_t> dis(0, max_value);
    std::normal_distribution<> dlen(args::get(range_mean),args::get(range_stdev));
    uint64_t x_len = args::get(test_size);
    uint64_t max_seen_value = 0;
#pragma omp parallel for
    for (int n=0; n<x_len; ++n) {
        uint64_t q = dis(gen);
        uint64_t r = std::min(q + std::max((uint64_t)0, (uint64_t)std::round(dlen(gen))), max_value);
        //max_seen_value = std::max(r, max_seen_value);
        //uint64_t a = std::min(q, r);
        //uint64_t b = std::max(q, r);
        //std::cerr << q << ", " << r << std::endl;
        //bb.add(intpair(q, r));
        //intervals.emplace_back();
        //interval i = intervals.back();
        //i.l = q;
        //i.r = r;
        intervals.push_back(interval(q, r));
    }
    //tree.index();

    faststabbing db(intervals, intervals.size(), max_value);
    //p_iitii db = bb.build(n_domains);
    //p_iitii db = bb.build();
//#pragma omp parallel for
    for (int n=0; n<max_value; ++n) {
        std::vector<interval*> ovlp = db.query(n);
        if (n % 1000 == 0) std::cerr << n << "\r";
        //std::cerr << n << " has " << ovlp.size() << " overlaps" << std::endl;
        for (auto& s : ovlp) {
            if (s->l > n || s->r < n) {
                std::cerr << "tree broken at " << n << std::endl;
            }
        }
    }
    
    //std::vector<intpair> results = db.overlap(22, 25);
    // alternative: db.overlap(22, 25, results);

    //for (const auto& p : results)
    //std::cout << p.first << "\t" << p.second << std::endl;
    return 0;
}
