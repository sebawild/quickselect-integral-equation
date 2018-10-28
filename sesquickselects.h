//
// Created by seb on 10/28/18.
//

#ifndef INTEGRAL_EQ_SESQUICKSELECTS_H
#define INTEGRAL_EQ_SESQUICKSELECTS_H


template<typename Num>
struct ClassicQuickselect : public QuickselectParams<Num> {
	Num partitioningCosts(const Num alpha) const override {
		return 1;
	}
	vector<unsigned int> samplingParameter(const Num alpha) const override {
		return {0, 0};
	}
	string name() override { return "classic-quickselect"; }
};

template<typename Num>
struct YBBSelectSE : public QuickselectParams<Num> {
	Num partitioningCosts(const Num alpha) const override {
		return 4./3;
	}
	vector<unsigned int> samplingParameter(const Num alpha) const override {
		return {0, 0, 0};
	}
	string name() override { return "ybb-select-t0-se"; }
};

template<typename Num>
struct WaterlooSelect : public QuickselectParams<Num> {
	Num partitioningCosts(const Num alpha) const override {
		return 1.5;
	}
	vector<unsigned int> samplingParameter(const Num alpha) const override {
		return {0,0,0,0};
	}
	string name() override { return "waterloo-quickselect-t0"; }
};



template<typename Num>
struct NonadaptiveQuickselect : public QuickselectParams<Num> {
	NonadaptiveQuickselect(vector<unsigned int> t, const Num a) : a(a), t(std::move(t)) {}
	const Num a;
	const vector<unsigned int> t;
	Num partitioningCosts(const Num alpha) const override { return a; }
	vector<unsigned int> samplingParameter(const Num alpha) const override { return t; }
	string name() override { return "nonadaptive-t"+to_string(t) + "-a"+to_string(a); }
};







template<typename Num>
struct SesquickselectK2 : public QuickselectParams<Num> {
	const Num nuStar = 0.265716848;
	vector<unsigned int> samplingParameter(const Num alpha) const override {
		if (alpha < nuStar) return {0, 1};
		if (alpha > 1 - nuStar) return {1, 0};
		// else
		return {0, 0, 0};
	}
	string name() override { return "sesquickselect-k2-optimal-nu"; }
};


template<typename Num> struct SesquickselectK3 : public QuickselectParams<Num> {
	explicit SesquickselectK3(vector<Num> nu = {0.1035,0.5}) : nu(std::move(nu)) {}
	const vector<Num> nu;
	vector<unsigned int> samplingParameter(const Num alpha) const override {
		if (alpha < nu[0])      return  {0,2};
		if (alpha < nu[1])      return {0,0,1};
		if (alpha <= 1 - nu[1]) return  {1,1}; // not used!
		if (alpha <= 1 - nu[0]) return {1,0,0};
		if (alpha > 1 - nu[0])  return  {2,0};
	}
	string name() override { return "sesquickselect-k3-nu" + to_string(nu); }
};

template<typename Num> struct SesquickselectK4 : public QuickselectParams<Num> {
	explicit SesquickselectK4(vector<Num> nu ={.06,.28,.5,.5}) : nu(std::move(nu)) {}
	const vector<Num> nu;
	vector<unsigned int> samplingParameter(const Num alpha) const override {
		if (alpha < nu[0])      return  {0,3};
		if (alpha < nu[1])      return {0,0,2};
		if (alpha < nu[2])      return {0,1,1};
		if (alpha < nu[3])      return {1,0,1}; // not used!
		if (alpha <= 1 - nu[3]) return {0,2,0}; // not used!
		if (alpha <= 1 - nu[2]) return {1,0,1}; // not used!
		if (alpha <= 1 - nu[1]) return {1,1,0};
		if (alpha <= 1 - nu[0]) return {2,0,0};
		if (alpha > 1 - nu[0])  return  {3,0};
	}
	string name() override { return "sesquickselect-k4-nu" + to_string(nu); }
};

template<typename Num> struct SesquickselectK5 : public QuickselectParams<Num> {
	explicit SesquickselectK5(vector<Num> nu = {0.036, 0.036, 0.153, 0.153, 0.5}) : nu(std::move(nu)) {}
	const vector<Num> nu;
	vector<unsigned int> samplingParameter(const Num alpha) const override {
		if (alpha < nu[0])      return  {0,4};
		if (alpha < nu[1])      return  {1,3}; // not used!
		if (alpha < nu[2])      return {0,0,3};
		if (alpha < nu[3])      return {1,0,2}; // not used!
		if (alpha < nu[4])      return {0,1,2};
		if (alpha <= 1 - nu[4]) return {1,1,1}; // not used! also {2,2} not better
		if (alpha <= 1 - nu[3]) return {2,1,0};
		if (alpha <= 1 - nu[2]) return {2,0,1}; // not used!
		if (alpha <= 1 - nu[1]) return {3,0,0};
		if (alpha <= 1 - nu[0]) return  {3,1}; // not used!
		if (alpha > 1 - nu[0])  return  {4,0};
	}
	string name() override { return "sesquickselect-k5-nu" + to_string(nu); }
};

template<typename Num> struct SesquickselectK6 : public QuickselectParams<Num> {
	explicit SesquickselectK6(vector<Num> nu = {0.025, 0.09, 0.38, 0.38, 0.5}) : nu(std::move(nu)) {}
	const vector<Num> nu;
	vector<unsigned int> samplingParameter(const Num alpha) const override {
		if (alpha < nu[0])      return  {0,5};  // min
		if (alpha < nu[1])      return {0,0,4}; // min
		if (alpha < nu[2])      return {0,1,3}; // min
		if (alpha < nu[3])      return {0,2,2}; // min // not used!
		if (alpha < nu[4])      return {1,1,2};
		if (alpha <= 1 - nu[4]) return {2,0,2}; // not used!
//		if (alpha <= 1 - nu[4]) return {1,2,1}; // not used!
		if (alpha <= 1 - nu[3]) return {2,1,1};
		if (alpha <= 1 - nu[2]) return {2,2,0}; // not used!
		if (alpha <= 1 - nu[1]) return {3,1,0};
		if (alpha <= 1 - nu[0]) return {4,0,0};
		if (alpha > 1 - nu[0])  return  {5,0};
	}
	string name() override { return "sesquickselect-k6-nu" + to_string(nu); }
};

template<typename Num> struct SesquickselectK7 : public QuickselectParams<Num> {
	explicit SesquickselectK7(vector<Num> nu = {0.02, 0.06, 0.2875, 0.2875, 0.465, 0.465, 0.5}) : nu(std::move(nu)) {}
	const vector<Num> nu;
	vector<unsigned int> samplingParameter(const Num alpha) const override {
		if (alpha < nu[0])      return  {0,6};
		if (alpha < nu[1])      return {0,0,5};
		if (alpha < nu[2])      return {0,1,4};
		if (alpha < nu[3])      return {0,2,3}; // not used!
		if (alpha < nu[4])      return {1,1,3}; // used! // guess not
		if (alpha < nu[5])      return {0,3,2}; // not used & guess not
		if (alpha < nu[6])      return {1,2,2};
		if (alpha <= 1 - nu[6])      return {3,3}; // not used! & guess not
//		if (alpha <= 1 - nu[6])      return {2,1,2}; // not used! & guess not
		if (alpha <= 1 - nu[5]) return {2,2,1};
		if (alpha <= 1 - nu[4]) return {2,3,0};
		if (alpha <= 1 - nu[3]) return {3,1,1};
		if (alpha <= 1 - nu[2]) return {3,2,0};
		if (alpha <= 1 - nu[1]) return {4,1,0};
		if (alpha <= 1 - nu[0]) return {5,0,0};
		if (alpha > 1 - nu[0])  return  {6,0};
	}
	string name() override { return "sesquickselect-k7-nu" + to_string(nu); }
};



typedef vector<unsigned int> tvec_t;

const tvec_t allTsK3[] = {{0, 2}, {1, 1}, {2, 0}, {0, 0, 1}, {0, 1, 0}, {1, 0, 0}};
const tvec_t allTsK4[] = {
		{0, 3}, {1, 2}, {2, 1}, {3, 0},
		{0, 0, 2}, {0, 1, 1}, {0, 2, 0}, {1, 0, 1}, {1, 1, 0}, {2, 0, 0}
};

const tvec_t allTsK5[] =
		{
				{0, 4}, {1, 3}, {2, 2}, //{3, 1}, {4, 0},
				{0, 0, 3}, {0, 1, 2}, {0, 2, 1}, {0, 3, 0}, {1, 0, 2},
				{1, 1, 1}, {1, 2, 0}, //{2, 0, 1}, {2, 1, 0}, {3, 0, 0}
		};
const tvec_t allTsK6[] = {
		{0, 5}, {1, 4}, {2, 3}, {3, 2}, {4, 1}, {5, 0},
		{0, 0, 4}, {0, 1, 3}, {0, 2, 2}, {0, 3, 1}, {0, 4, 0}, {1, 0, 3},
		{1, 1, 2}, {1, 2, 1}, {1, 3, 0}, {2, 0, 2}, {2, 1, 1}, {2, 2, 0},
		{3, 0, 1}, {3, 1, 0}, {4, 0, 0}};
const tvec_t allTsK7[] = {
		{0, 6}, {1, 5}, {2, 4}, {3, 3}, //{4, 2}, {5, 1}, {6, 0},
		{0, 0, 5}, {0, 1, 4}, {0, 2, 3}, {0, 3, 2}, {0, 4, 1}, {0, 5, 0},
		{1, 0, 4}, {1, 1, 3}, {1, 2, 2}, {1, 3, 1}, {1, 4, 0}, {2, 0, 3},
		{2, 1, 2}, {2, 2, 1}, {2, 3, 0},
//		{3, 0, 2}, {3, 1, 1}, {3, 2, 0}, {4, 0, 1}, {4, 1, 0}, {5, 0, 0}
};
const tvec_t allTsK9[] = {
		{0,0,7},{0,1,6},{0,2,5},{0,3,4},{0,4,3},{0,5,2},{0,6,1},{0,7,0},{1,0,6},{1,1,5},{1,2,4},{1,3,3},{1,4,2},{1,5,1},{2,0,5},{2,1,4},{2,2,3},{2,3,2},{3,0,4},{3,1,3},{0,8},{1,7},{2,6},{3,5},{4,4}
};



#endif //INTEGRAL_EQ_SESQUICKSELECTS_H
