#include "rv_t_walk.hpp"

#include <boost/math/special_functions/gamma.hpp>
#include <iostream>
#include <math.h>
#include <limits>
#include <stdexcept>

const double pi = boost::math::constants::pi<double>();


RV_t_walk::RV_t_walk(
		double const & x1_,
		double			 & X,
		double const & os_,
		double const & p1_,
		double const & s1_,
		trng::yarn2  & R_
) : x1(x1_), x2(X), os(os_), p1(p1_), s1(s1_), 
		R(R_) 
{ 
	if (p1 <= 0.0) {
		std::stringstream msg;
		msg << "Given value for 'p1' is " << p1;
		msg << " degrees of freedom must be a positive number.";
		throw std::domain_error(msg.str());
	}
	if (s1 <= 0.0) {
		std::stringstream msg;
		msg << "Given value for 's1' is " << s1;
		msg << " scale must be a positive number.";
		throw std::domain_error(msg.str());
	}
}

std::map<std::string, double> RV_t_walk::state() const {
	std::map<std::string, double> out;
	out["x1"] = x1;
	out["X"] = x2;
  out["os"] = os,
	out["p1"] = p1;
	out["s1"] = s1;
	return out;
}

void RV_t_walk::jump(double X) { 
	x2 = X; 
}

double RV_t_walk::draw() {
	double S, V, W, T;
	do {
		S = 2.0 * U(R) - 1.0;
		V = 2.0 * U(R) - 1.0;
		W = S*S + V*V;
	} while (W > 1.0);
	T = S * sqrt(p1*(pow(W,-2.0/p1)-1.0)/W);
  x2 = T * s1 + x1 + os;
	return x2; 
}

double RV_t_walk::lpdf(double X) {
	double lpdf;
	lpdf = boost::math::lgamma((p1+1.0)/2.0) - boost::math::lgamma(p1/2.0) -
		0.5 * log(p1*pi*pow(s1,2)) -
		(p1+1.0)/2.0 * log(1.0 + (pow(X-(x1+os),2))/(p1*pow(s1,2)) );
	return lpdf;
}

double RV_t_walk::lpdf() { return lpdf(x2); }

