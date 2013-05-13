#include "location.hpp"
#include <math.h>
#include <boost/math/special_functions/gamma.hpp>
#include <iostream>

using boost::math::lgamma;
const double pi = boost::math::constants::pi<double>();

t_walk_Posterior::t_walk_Posterior(
		double const & x1_,
		double 			 & x2_,
		double const & x3_,
		double const & p1_,   // degrees of freedom 1
		double const & p2_,   // degrees of freedom 2
		double const & s1_,   // scale 1
		double const & s2_,   // scale 2
		trng::yarn2 & R_
) :	x1(x1_), x2(x2_), x3(x3_),
		p1(p1_), p2(p2_), s1(s1_), s2(s2_), R(R_) 
{
	std::cout << x1 << std::endl;
	std::cout << x2 << std::endl;
}

std::map<std::string, double> t_walk_Posterior::state() const {
	std::map<std::string, double> out;
	out["x1"] = x1;
	out["x2"] = x2;
	out["x3"] = x3;
	out["p1"] = p1;
	out["p2"] = p2;
	out["s1"] = s1;
	out["s2"] = s2;
	return out;
}

void t_walk_Posterior::jump(double x2_) {
	x2 = x2_;
}
// Don't divide doubles by integers, etc...
double t_walk_Posterior::lpdf(double x2_) {
	double lpdf;
	lpdf = 	lgamma((p1+1.0)/2.0) - lgamma(p1/2.0) -
		0.5 * log(p1*pi*pow(s1,2)) - 
		(p1+1)/2 * log(1 + (pow(x2_-x1,2))/(p1*pow(s1,2)) ) +
					lgamma((p2+1)/2) - lgamma(p2/2) -
		0.5 * log(p2*pi*pow(s2,2)) - 
		(p2+1)/2 * log(1 + (pow(x3-x2_,2))/(p2*pow(s2,2)) );
	return lpdf;
}

double t_walk_Posterior::lpdf() { return lpdf(x2); }

