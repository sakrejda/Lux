#include "random.hpp"
#include <math.h>
#include <cmath>
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
		companion(3,3), R(R_), eigvalues(3),
		bounds_pk1(2), bounds_pk2(2), EXPO(1.0)
{
	companion(1,0) = 1.0;
	companion(2,1) = 1.0;
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
  ly = lpdf() - EXPO(R);
	find_slice();
	while(true) {
		x_new = choose();	
		ly_new = lpdf(x_new);
		if (ly_new > ly) // && (x_new > min) && (x_new < max))  truncation.
			break; 
		else 
			trim();
	}
	x2 = x_new;
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

void RV_Missing_t_walk::find_slice() {
	find_peaks();
	bounds_pk1 = step_out(peak1);
	bounds_pk2 = step_out(peak2);
	if (
		(bounds_pk1[0] < std::min(bounds_pk1[1], bounds_pk2[1])) &&
		(bounds_pk2[0] < std::min(bounds_pk1[1], bounds_pk2[1]))
	) {
		double x = (bounds_pk1[0]+bounds_pk2[1])/2.0;
		bounds_pk1[1] = x;
		bounds_pk2[0] = x;
	}
}

std::vector<double> RV_Missing_t_walk::step_out(double peak) {
	std::vector<double> bounds(2);
	int m = 10;
	double w = (s1+s2)/2;
	bounds[0] = peak - w * U(R);  // w = use (s1+s2)/2
	bounds[1] = bounds[0] + w;
	int j = std::floor(m * U(R));    // m needed...
	int k = (m-1) - j;
	while ((j>0) && lpdf() < lpdf(bounds[0]) ) { // trunc. can be added here.
		bounds[0] = bounds[0] - w;
		j = j - 1;
	}
	while ((k>0) && lpdf() < lpdf(bounds[1]) ) { // truncation can be added here.
		bounds[1] = bounds[1] + w;
		k = k - 1;
	}
	return bounds;
}

void RV_Missing_t_walk::trim() {
	if (x_new <= bounds_pk1[1]) {
		if (x_new <= x2) {
			bounds_pk1[0] = x_new;
		} else {
			bounds_pk1[1] = x_new;
		}
	} else {
		if (x_new <= x2) {
			bounds_pk2[0] = x_new;
		} else {
			bounds_pk2[1] = x_new;
		}
	}
}

double RV_Missing_t_walk::choose() {
	double d = (bounds_pk2[1] - bounds_pk2[0]) + (bounds_pk1[1] - bounds_pk1[0]);
	double q = (bounds_pk2[0] - bounds_pk1[1]);
	double l = U(R) * d;
	if ((bounds_pk1[0] + l) > bounds_pk1[1]) {
		return (l + q + bounds_pk1[0]);
	} else {
		return (l + bounds_pk1[0]);
	}
}

void RV_Missing_t_walk::find_peaks() {
	companion(0,2) = 0.5 * (x1*x1*x3 + x1*p2*s2*s2 + x3*p1*s1*s1 + x3*x3*x1);
	companion(1,2) = -0.5 * (p1*s1*s1 + x1*x1 + 4.0*x1*x3 + x3*x3 + p2*s2*s2);
	companion(2,2) = 1.5 * (x1+x3);
	arma::eig_gen(cx_eigval, cx_eigvec, companion);	
	eigvalues = arma::real(cx_eigval);
	eigvalues = arma::sort(eigvalues);
	// Calculate condition number, if it's < 10^(-6) or so, 
	// set both peaks to the average value... 
	peak1 = eigvalues[0];
	peak2 = eigvalues[2];
	//	std::cout << "Peak 1: " << peak1 << ", peak 2: " << peak2 << std::endl;
}


