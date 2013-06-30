#include "random.hpp"
#include <math.h>
#include <boost/math/special_functions/gamma.hpp>
#include <iostream>

#include <limits>
#include <trng/uniform01_dist.hpp>

using boost::math::lgamma;
const double pi = boost::math::constants::pi<double>();

RV_Constant::RV_Constant(double & X) : x(X) {}

void RV_Constant::jump(double X) { }

double RV_Constant::draw() {return x;}

double RV_Constant::lpdf(double X) {
	if (X == x)
		return std::numeric_limits<double>::infinity();
	else 
		return 0.0;
}

double RV_Constant::lpdf() { return std::numeric_limits<double>::infinity(); }



RV_Uniform::RV_Uniform(
		double & X, 
		double const & minimum, double const & maximum, 
		trng::yarn2 & R_
) : x(X), min(minimum), max(maximum), R(R_) { }

void RV_Uniform::jump(double X) {
	if ((X > min) && (X < max))
		x = X;
}

double RV_Uniform::draw() { 
	x = min + (U(R) * (max - min)); 
/*	std::cout << "x:   " << x << ", ";
	std::cout << "min: " << min << ", ";
	std::cout << "max: " << max << ", " << std::endl; */
	return x;
}

double RV_Uniform::lpdf(double X) {
	if ((X > min) && (X < max)) 
		return (1/(max-min));
	else
		return 0.0;
}

double RV_Uniform::lpdf() { 
	if ((x > min) && (x < max))
		return (1/(max-min)); 
	else
		return 0;
}



RV_Missing_t_walk::RV_Missing_t_walk(
		double const & x1_,
		double 			 & X,
		double const & x3_,
		double const & p1_,   // degrees of freedom 1
		double const & p2_,   // degrees of freedom 2
		double const & s1_,   // scale 1
		double const & s2_,   // scale 2
		trng::yarn2 & R_
) :	x1(x1_), x2(X), x3(x3_),
		p1(p1_), p2(p2_), s1(s1_), s2(s2_), 
		companion(3,3), R(R_), eigvalues(3)
{
	companion(1,0) = 1.0;
	companion(2,1) = 1.0;
	find_peaks();
}

std::map<std::string, double> RV_Missing_t_walk::state() const {
	std::map<std::string, double> out;
	out["x1"] = x1;
	out["X"] = x2;
	out["x3"] = x3;
	out["p1"] = p1;
	out["p2"] = p2;
	out["s1"] = s1;
	out["s2"] = s2;
	return out;
}

void RV_Missing_t_walk::jump(double X) { x2 = X; }

double RV_Missing_t_walk::draw() { 
	// Slice sampler here, complicaterated a little by the fact that the
	// slice will usually/often consist of two parts.	
	
	return x2; 
}

double RV_Missing_t_walk::lpdf(double X) {
	double lpdf;
	lpdf = 	lgamma((p1+1.0)/2.0) - lgamma(p1/2.0) -
		0.5 * log(p1*pi*pow(s1,2)) - 
		(p1+1.0)/2.0 * log(1.0 + (pow(X-x1,2))/(p1*pow(s1,2)) ) +
					lgamma((p2+1.0)/2.0) - lgamma(p2/2.0) -
		0.5 * log(p2*pi*pow(s2,2)) - 
		(p2+1.0)/2.0 * log(1.0 + (pow(x3-X,2))/(p2*pow(s2,2)) );
	return lpdf;
}

double RV_Missing_t_walk::lpdf() { return lpdf(x2); }

void find_slice() {
	find_peaks()
	// set l_bound/r_bound 1 and 2 starting from peak1 and peak2!!!
}

void RV_Missing_t_walk::find_peaks() {
	companion(0,2) = 0.5 * (x1*x1*x3 + x1*p2*s2*s2 + x3*p1*s1*s1 + x3*x3*x1);
	companion(1,2) = -0.5 * (p1*s1*s1 + x1*x1 + 4.0*x1*x3 + x3*x3 + p2*s2*s2);
	companion(2,2) = 1.5 * (x1+x3);
	arma::eig_gen(cx_eigval, cx_eigvec, companion);	
	eigvalues = arma::real(cx_eigval);
	eigvalues = arma::sort(eigvalues);
	peak1 = eigvalues[0];
	peak2 = eigvalues[2];
}


