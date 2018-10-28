/*
 * Numerical solution of the integral equation to approximate f(alpha)
 * for adaptive Quickselect variants.
 *
 * There is not real user interface; change the code below to use the
 * algorithm variant of interest and recompile...
 * The output is written to file.
 * Algorithms are defined in sesquickselect.h and proportion_from_k.h.
 */

#include <iostream>
#include <utility>
#include <vector>
#include "integral_helper.h"
#include "proportion_from_k.h"
#include "sesquickselects.h"

//typedef boost::rational<long> num;
typedef long double num;

int main() {
	unsigned int nSamples = 200;
	unsigned int nIterations = 50;
	double initValue = 1;

	// SELECT ALGORITHM HERE
	// default: sesquickselect with optimal cutoff, costs are scanned elements
	QuickselectParams<num> * pParams = new SesquickselectK2<num> {};

	string filename = pParams->name()
	                  + "-samples" + to_string(nSamples)
	                  + "-iter" + to_string(nIterations)
	                  + "-init" + to_string(initValue);
	vector<num> fOld(nSamples, initValue), f(nSamples);
	std::cout << filename << std::endl;
	for (int r = 1; r <= nIterations; ++r) {
		cout << "Running Iteration " << r << " ..." << endl;
		integralEquationIteration(pParams, fOld, f);
		swap(fOld, f); // yes, expensive copy, but s dominated by above computations anyway ...
		printXYAsMathematicaList(f);
	}
	writeXYtoTableFile(f, filename + ".tab");

	return 0;
}
