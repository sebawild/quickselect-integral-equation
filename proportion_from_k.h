//
// Created by seb on 7/26/18.
//

#ifndef INTEGRAL_EQ_PROPORTION_FROM_K_H
#define INTEGRAL_EQ_PROPORTION_FROM_K_H


template<typename Num> struct ProportionFromK2 : public QuickselectParams<Num> {
	vector<unsigned int> samplingParameter(const Num alpha) const override {
		if (alpha < 0.5)   return {0,1};
		if (alpha >= 0.5)  return {1,0};
	}
	Num partitioningCosts(const Num alpha) const override { return 1; }
	string name() override { return string("proportion-from-k2"); }
};


template<typename Num> struct ProportionFromK3 : public QuickselectParams<Num> {
	ProportionFromK3(const Num nu1 = 0.182000001) : nu1(nu1) {}
	const Num nu1;
	vector<unsigned int> samplingParameter(const Num alpha) const override {
		if (alpha < nu1)      return {0,2};
		if (alpha <= 1 - nu1) return {1,1};
		if (alpha > 1 - nu1)  return {2,0};
	}
	Num partitioningCosts(const Num alpha) const override {
		return 1;
	}
	string name() override {
		return string("proportion-from-k3")
		       + "-nu1" + to_string(nu1)
				;
	}
};

template<typename Num> struct ProportionFromK4 : public QuickselectParams<Num> {
	ProportionFromK4(const Num nu1 = 0.09) : nu1(nu1) {}
	const Num nu1;
	vector<unsigned int> samplingParameter(const Num alpha) const override {
		if (alpha < nu1)      return {0,3};
		if (alpha < 0.5)      return {1,2};
		if (alpha <= 1 - nu1) return {2,1};
		if (alpha > 1 - nu1)  return {3,0};
	}
	Num partitioningCosts(const Num alpha) const override {
		return 1;
	}
	string name() override {
		return string("proportion-from-k4")
		       + "-nu1" + to_string(nu1)
				;
	}
};

template<typename Num> struct ProportionFromK5 : public QuickselectParams<Num> {
	ProportionFromK5(const Num nu1 = 0.0510000001, const Num nu2 = 0.3050000001) : nu1(nu1), nu2(nu2) {}

	const Num nu1, nu2;
	vector<unsigned int> samplingParameter(const Num alpha) const override {
		if (alpha < nu1)      return {0,4};
		if (alpha < nu2)      return {1,3};
		if (alpha <= 1 - nu2) return {2,2};
		if (alpha <= 1 - nu1) return {3,1};
		if (alpha > 1 - nu1)  return {4,0};
	}
	Num partitioningCosts(const Num alpha) const override {
		return 1;
	}
	string name() override {
		return string("proportion-from-k5")
		       + "-nu1" + to_string(nu1)
		       + "-nu2" + to_string(nu2)
				;
	}
};

template<typename Num> struct ProportionFromK6 : public QuickselectParams<Num> {
	ProportionFromK6(const Num nu1 = 0.033, const Num nu2 = 0.2045) : nu1(nu1), nu2(nu2) {}
	const Num nu1, nu2;
	vector<unsigned int> samplingParameter(const Num alpha) const override {
		if (alpha < nu1)      return {0,5};
		if (alpha < nu2)      return {1,4};
		if (alpha < 0.5)      return {2,3};
		if (alpha <= 1 - nu2) return {3,2};
		if (alpha <= 1 - nu1) return {4,1};
		if (alpha > 1 - nu1)  return {5,0};
	}
	Num partitioningCosts(const Num alpha) const override {
		return 1;
	}
	string name() override {
		return string("proportion-from-k6")
		       + "-nu1" + to_string(nu1)
		       + "-nu2" + to_string(nu2)
				;
	}
};



template<typename Num> struct ProportionFromK7 : public QuickselectParams<Num> {
	ProportionFromK7(const Num nu1 = 0.027, const Num nu2 = 0.15, const Num nu3 = 0.365) :
			nu1(nu1), nu2(nu2), nu3(nu3) {}
	const Num nu1, nu2, nu3;
	vector<unsigned int> samplingParameter(const Num alpha) const override {
		if (alpha < nu1)      return {0,6};
		if (alpha < nu2)      return {1,5};
		if (alpha < nu3)      return {2,4};
		if (alpha <= 1 - nu3) return {3,3};
		if (alpha <= 1 - nu2) return {4,2};
		if (alpha <= 1 - nu1) return {5,1};
		if (alpha > 1 - nu1)  return {6,0};
	}
	string name() override {
		return string("proportion-from-k7")
		       + "-nu1" + to_string(nu1)
		       + "-nu2" + to_string(nu2)
		       + "-nu3" + to_string(nu3)
				;
	}
};



template<typename Num> struct ProportionFromK9 : public QuickselectParams<Num> {
	ProportionFromK9(const Num nu1 = 0.014500001, const Num nu2 = 0.0900001, const Num nu3 = 0.2200001,
	                 const Num nu4 = 0.400001) :
			nu1(nu1), nu2(nu2), nu3(nu3), nu4(nu4) {}
	const Num nu1, nu2, nu3, nu4;
	vector<unsigned int> samplingParameter(const Num alpha) const override {
		if (alpha < nu1)      return  {0,8};
		if (alpha < nu2)      return {1,7};
		if (alpha < nu3)      return {2,6};
		if (alpha < nu4)      return {3,5};
		if (alpha <= 1 - nu4) return {4,4};
		if (alpha <= 1 - nu3) return {5,3};
		if (alpha <= 1 - nu2) return {6,2};
		if (alpha <= 1 - nu1) return {7,1};
		if (alpha > 1 - nu1)  return  {8,0};
	}
	Num partitioningCosts(const Num alpha) const override {
		return 1;
	}
	string name() override {
		return string("proportion-from-k9")
		       + "-nu1" + to_string(nu1)
		       + "-nu2" + to_string(nu2)
		       + "-nu3" + to_string(nu3)
		       + "-nu4" + to_string(nu4)
				;
	}
};



#endif //INTEGRAL_EQ_PROPORTION_FROM_K_H
