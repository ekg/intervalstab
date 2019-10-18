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

// define OGDF_DEBUG in the debug mode of your IDE

#include "FastStabbing.h"
#include <vector>
#include <ogdf/basic/basic.h>
#include <fstream>
#include <time.h>
#include <math.h>


using namespace ogdf;

// main method
int main(int argc, char* argv[]) {
	cout.precision(3);
	cout.setf(cout.fixed);
//	srand(unsigned(time(NULL))); // init randomizer
	
	int randomMin = 10, randomMax = 100, randomStep = 1, rangemultiplicator=2, iterations = 500, numberQueries = 50; // presets
	double chazelle = 0; // invoke chazelle algorithm on value>0 with delta=value
	bool onlySearch = false; // returns the first interval stabbed, if there is some, an empty list otherwise
	bool exponentialLength = false; // length of random intervals is exponentially distributed (expected value=1000)
	char logFile[255]; // name of logFile
	sprintf_s(logFile, "%s", "Results/output.txt");

	for (int i = 0; i < argc; i++) {
		const char* arg = argv[i];

		if (strcmp(arg,"-n") == 0) {
			if (i == argc-3) {
				cerr << "Missing additional arguments for -n!" << endl;
			} else {
				sscanf(argv[++i], "%d", &randomMin);
				sscanf(argv[++i], "%d", &randomMax);
				sscanf(argv[++i], "%d", &randomStep);
			}
		} else if (strcmp(arg,"-rangemultiplicator") == 0) {
			if (i == argc-1) {
				cerr << "Missing additional arguments for -rangemultiplicator!" << endl;
			} else {
				sscanf(argv[++i], "%d", &rangemultiplicator);
			}
		} else if (strcmp(arg,"-chazelle") == 0) {
			if (i == argc-1) {
				cerr << "Missing additional arguments for -chazelle!" << endl;
			} else {
				sscanf(argv[++i], "%lf", &chazelle);
			}
		} else if (strcmp(arg,"-onlySearch") == 0) {
			onlySearch = true;
		} else if (strcmp(arg,"-exponentialLength") == 0) {
			exponentialLength = true;
		} else if (strcmp(arg,"-queries") == 0) {
			if (i == argc-1) {
				cerr << "Missing additional arguments for -queries!" << endl;
			} else {
				sscanf(argv[++i], "%d", &numberQueries);
			}
		} else if (strcmp(arg,"-iterations") == 0) {
			if (i == argc-1) {
				cerr << "Missing additional arguments for -iterations!" << endl;
			} else {
				sscanf(argv[++i], "%d", &iterations);
			}
		} else if (strcmp(arg,"-log") == 0) {
			if (i == argc-1) {
				cerr << "Missing additional argument for -log!" << endl;
			} else {
				sprintf_s(logFile, "%s", argv[++i]);
			}
		}
	}

	std::ofstream log(logFile);
	Array<interval> a(randomMax);
	Array<int> queries(numberQueries);
	std::vector<interval*> output;

	// test with n intervals in the range [bigN]
	int i,n,iter,j=0;
	for (n = randomMin; n <= randomMax; n = n+randomStep) {
		for (iter = 0; iter<iterations; ++iter) {
			// generate random input
			cout << "Chazelle=" << chazelle << ", exp=" << exponentialLength << ", generate random input: " << flush;
			bool identical;
			do {
				for (i = 0; i < n; ++i) {
					do {
						if (exponentialLength) { // exponential distributed with expected value 1000
							do {
								a[i].l = randomNumber(1,n*rangemultiplicator);
								a[i].r = a[i].l+(int)(-1000*::log(randomDouble(0,1)));
							} while (a[i].r > n*rangemultiplicator);
						} else {
							a[i].l = randomNumber(1,n*rangemultiplicator);
							a[i].r = randomNumber(1,n*rangemultiplicator);
						}
					} while (a[i].l > a[i].r); // intervals [a,a] are allowed
					a[i].parent = NULL;
					a[i].leftsibling = NULL;
					a[i].rightchild = NULL;
					a[i].pIt = NULL;
					a[i].smaller = NULL;
					#ifdef OGDF_DEBUG
					a[i].stabbed = false;
					#endif
				}
				//cout << a << "\n";
				// sort input
				cout << "sort input... " << flush;
				if (chazelle < 0.00001) {
					IntervalComparerInv comp;
					a.quicksortCT(0,n-1,comp);
				} else {
					IntervalComparer comp;
					a.quicksortCT(0,n-1,comp);
				}
				identical = false;
				for (i = 1; i < n; ++i) {
					if (a[i].l == a[i-1].l && a[i].r == a[i-1].r) {
						identical = true;
						break;
					}
				}
			} while (identical && !exponentialLength); // while identical intervals exist (not on exponentialLength!)
			//cout << "\n" << a << "\n";

			// create random query array
			cout << "\ncreate random queries...\n" << flush;
			for (j = 0; j < numberQueries; j++) {
				queries[j] = randomNumber(1,n*rangemultiplicator);
			}

			cout << "preprocessing...\t" << n << "\t" << n*rangemultiplicator << "\t" << iter << "\n" << flush;
			double time,secPre,secQueries = 0;
			long outputSize = 0;
			long numComparisons = 0;

			// preprocessing
			usedTime(time);
			if (chazelle < 0.00001) {
				// my algorithm
				FastStabbing stabbing(a,n,n*rangemultiplicator);
				secPre = usedTime(time);
				cout << numberQueries << " queries...\n\n" << flush;

				// queries
				usedTime(time);
				for (j = 0; j < numberQueries; ++j) {
					stabbing.query(queries[j],output,onlySearch,numComparisons);
					outputSize = outputSize+(int)output.size();
				}
			} else {
				// chazelle algorithm
				ChazelleStabbing stabbing(a,n,n*rangemultiplicator,chazelle);
				secPre = usedTime(time);
				cout << numberQueries << " queries...\n\n" << flush;

				// queries
				usedTime(time);
				for (j = 0; j < numberQueries; ++j) {
					stabbing.query(queries[j],output,onlySearch,numComparisons);
					outputSize = outputSize+(int)output.size();
				}
			}
			secQueries = usedTime(time);
			log << n << "\t" << iter << "\t" << n*rangemultiplicator << "\t" << secPre << "\t" << secQueries
				<< "\t" << numberQueries << "\t" << outputSize << "\t" << numComparisons << endl << flush;
		}
	}

	#ifdef OGDF_DEBUG
	cout << "\npress any key...";
	getchar();
	#endif
	return 0;
}
