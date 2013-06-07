#ifndef THETA_H
#define THETA_H

#include <map>
#include <armadillo>

class Theta {

public:
	Theta(double x);

private:
	std::map<std::string, double> mp;

};

class TH_Constant : public Theta {

public:
	TH_Constant(double x);

private:
	std::map<std::string, double> mp;

};

class TH_Uniform : public Theta {
	TH_Uniform(double minimum, double maximum, double x);

};

class TH_Missing_t_walk : public Theta {

public:
	TH_Uniform(double x1, double x2, double x3, 
			double s1, double s2, double p1, double p2);

private:
	std::map<std::string, double> mp;
	
};


#endif
