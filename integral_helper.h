#ifndef INTEGRAL_EQ_COMPUTATIONS_H
#define INTEGRAL_EQ_COMPUTATIONS_H

#include <iostream>
#include <cmath>
#include <vector>
#include <memory>
#include <cassert>
#include <fstream>

using namespace std;



template<typename Num>
string to_string(vector<Num> v) {
	string res = "[";
	for (unsigned int i = 0; i < v.size(); ++i)
		res += (i == 0 ? "" : ",") + to_string(v[i]);
	return res + "]";
}



template<typename Num>
Num scannedElements(vector<unsigned int> t) {
//	cout << to_string(t) << endl;
	switch (t.size()) {
		case 2: return 1;
		case 3: {
			unsigned int k = t[0] + t[1] + t[2] + 2;
			return 1. + min((t[0] + 1.) / (k + 1), (t[2] + 1.) / (k + 1));
		}
		default:
			throw exception();
	}
}


/** interface to implement by algorithm variants */
template<typename Num>
struct QuickselectParams {
	virtual Num partitioningCosts(Num alpha) const {
		return scannedElements<Num>(samplingParameter(alpha));
	}
	virtual vector<unsigned int> samplingParameter(Num alpha) const = 0;
	virtual string name() { return "some-quickselect"; }
};





template<typename Num>
Num pow(const Num x, unsigned int i) {
	if (i == 0) return 1;
	if (i == 1) return x;
	return x * pow(x, i-1);
}

template<typename Num>
Num beta(unsigned int a, unsigned int b) {
	if (a > b) return beta<Num>(b,a);
	// a <= b
	Num res = 1.0 / b;
	unsigned int denom = b+1;
	for (unsigned int i = 1; i < a; ++i) {
		Num num = i;
		res *= (num / denom++);
	}
	return res;
}

template<typename Num>
Num beta(unsigned int const a, unsigned int const b, unsigned int const c) {
	unsigned int C = max(max(a,b),c);
	unsigned int A = min(min(a,b),c);
	unsigned int B = a+b+c - A - C;
	// A <= B <= C
	Num res = 1.0 / C;
	res *= 1.0 / ++C;
	unsigned int denom = C+1;
	for (unsigned int i = 1; i < B; ++i) {
		Num num = i;
		res *= (num / denom++);
	}
	for (unsigned int i = 1; i < A; ++i) {
		Num num = i;
		res *= (num / denom++);
	}
	return res;
}




unsigned int tLeft(vector<unsigned int> t, int l) {
	assert(l >= 0);
	unsigned int res = 0;
	for (int i = 0; i < l-1; ++i) {
		res += t[i] + 1;
	}
	return res - 1;
}

unsigned int tRight(vector<unsigned int> t, int l) {
	assert (l < t.size());
	unsigned int res = 0;
	for (int i = l; i < t.size(); ++i) {
		res += t[i] + 1;
	}
	return res - 1;
}


template<typename Num>
Num alphaForIndex(const long i, const long nSamples) {
	Num num = 2*i+1;
	return num / (2*nSamples);
}

template<typename Num>
long indexForAlpha(const Num alpha, const long nSamples) {
	// alpha in [0,1]
	// alpha * (L-1) in [0,L-1]
	Num x = alpha * nSamples;
	auto l = (long) floor((double) x);
	return l >= nSamples ? nSamples-1 : l;
}



template<typename Num>
void printXYAsMathematicaList(const vector<Num> &values) {
	long L = values.size();
	cout << "{";
	for (long i = 0; i < L; ++i) {
		if (i > 0) cout << ",";
		cout << "{" << alphaForIndex<double>(i, L) << "," << values[i] << "}";
	}
	cout << "}" << endl;
}

template<typename Num>
void writeXYtoTableFile(const vector<Num> &values, const string &filename) {
	ofstream fout;
	fout.open(filename, fstream::out );
	fout.precision(20);
	long L = values.size();
	fout << "alpha\tf(alpha)" << endl;
	for (long i = 0; i < L; ++i) {
		fout << alphaForIndex<double>(i, L) << "\t" << values[i] << endl;
	}
	fout.close();
}




template<typename Num>
void integralEquationIteration(const QuickselectParams<Num> * quickselectParams,
                               const vector<Num> & fOld, vector<Num> & fNewOut) {
	const long L = fNewOut.size();
	const long LOld = fOld.size();
	const Num Delta = 1.0 / L;

	for (int i = 0; i < L; ++i) {
		auto alpha = alphaForIndex<Num>(i, L);
		const vector<unsigned int> t = quickselectParams->samplingParameter(alpha);
		auto s = static_cast<int>(t.size());

		// partitioning costs
		Num newValue = quickselectParams->partitioningCosts(alpha);
		// first integral (left call)
		{
			const unsigned int t1 = t[0], t2 = tRight(t, 1) ;
			Num step = Delta * (1 - alpha);
			const Num commonFactor = step / beta<Num>(t1 + 1, t2 + 1);
			Num sum = 0;
			for (int j = 0; j < (1-alpha) / step; ++j) {
				Num u = alpha + step/2 + j*step;
//			}
//			for (Num u = alpha + step / 2; u < 1; u += step) {
				Num summand = fOld[indexForAlpha(alpha / u, LOld)]
				              * pow(u,t1+1) * pow(1-u,t2);
				sum += summand;
			}
			newValue += sum * commonFactor;
		}
		// second integral (right call)
		{
			const unsigned int t1 = t[s-1], t2 = tLeft(t, s) ;
			Num step = Delta * alpha;
			const Num commonFactor = step / beta<Num>(t1 + 1, t2 + 1);
			Num sum = 0;
			for (int j = 0; j < alpha / step; ++j) {
				Num v = 0 + step/2 + j*step;
//			for (Num v = 0 + step / 2; v < alpha; v += step) {
				Num summand = fOld[indexForAlpha((alpha-v) / (1-v), LOld)]
				              * pow(1-v,t1+1) * pow(v,t2);
				sum += summand;
			}
			newValue += sum * commonFactor;
		}
		// third integrals (middle calls)
		for (int l = 2; l <= s-1; ++l) {
			const unsigned int t1 = tLeft(t,l), t2 = t[l-1], t3 = tRight(t,l) ;
			Num step1 = Delta * alpha, step2 = Delta * (1-alpha);
			const Num commonFactor = step1 * step2 / beta<Num>(t1 + 1, t2 + 1, t3 + 1);
			Num sum = 0;
			for (int j1 = 0; j1 < alpha / step1; ++j1) {
				for (int j2 = 0; j2 < (1 - alpha) / step2; ++j2) {
					Num u = 0 + step1 / 2 + j1 * step1;
					Num v = alpha + step2 / 2 + j2 * step2;
					Num summand = fOld[indexForAlpha((alpha - u) / (v - u), LOld)]
					              * pow(u,t1) * pow(v - u, t2 + 1) * pow(1-v, t3) ;
					sum += summand;
				}
			}
			newValue += sum * commonFactor;
		}

		fNewOut[i] = newValue;
	}

}



#endif //INTEGRAL_EQ_COMPUTATIONS_H
