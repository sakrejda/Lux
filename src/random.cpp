#include "random.hpp"
//#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
#include <iostream>
#include <math.h>

//#include <slicer-continuous.hpp>

#include <limits>
//#include <Eigen/Dense>
//#include <Eigen/Eigenvalues>

//using boost::math::lgamma;
const double pi = boost::math::constants::pi<double>();

RV_Constant::RV_Constant(double & X) : x(X) {}

void RV_Constant::jump(double X) { }

double RV_Constant::draw() {return x;}

double RV_Constant::lpdf(double X) {
	if (X == x)
		return 0.0;
	else 
		return -1.0 * std::numeric_limits<double>::infinity();
}

double RV_Constant::lpdf() { return lpdf(x); }



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

RV_Normal::RV_Normal(
		double       & X_,
		double const & mu_,
		double const & s_,
		trng::yarn2  & R_
) : X(X_), mu(mu_), s(s_), R(R_), NORMAL(mu_, s_) { }

void RV_Normal::jump(double X) { X = X; }

double RV_Normal::draw() {
	X = NORMAL(R);
	return X;
}

double RV_Normal::lpdf(double X) {
	NORMAL.mu(mu);
	NORMAL.sigma(s);
	return log(NORMAL.pdf(X));
}

double RV_Normal::lpdf() { lpdf(X); }

RV_t_walk::RV_t_walk(
		double const & x1_,
		double			 & X,
		double const & os_,
		double const & p1_,
		double const & s1_,
		trng::yarn2  & R_
) : x1(x1_), x2(X), os(os_), p1(p1_), s1(s1_), 
		R(R_), EXPO(1.0), bounds(2) { }

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
		
//		lgamma((p1+1.0)/2.0) - lgamma(p1/2.0) -
//		0.5 * log(p1*pi*pow(s1,2)) - 
//		(p1+1.0)/2.0 * log(1.0 + (pow(X-x1,2))/(p1*pow(s1,2)) ) +
//					lgamma((p2+1.0)/2.0) - lgamma(p2/2.0) -
//		0.5 * log(p2*pi*pow(s2,2)) - 
//		(p2+1.0)/2.0 * log(1.0 + (pow(x3-X,2))/(p2*pow(s2,2)) );
	return lpdf;
}

double RV_t_walk::lpdf() { return lpdf(x2); }


RV_Missing_t_walk::RV_Missing_t_walk(
		double const & x1_,
		double 			 & X,
		double const & x3_,
		double const & os1_,
		double const & os2_,
		double const & p1_,   // degrees of freedom 1
		double const & p2_,   // degrees of freedom 2
		double const & s1_,   // scale 1
		double const & s2_,   // scale 2
		trng::yarn2 & R_
) :	x1(x1_), x2(X), x3(x3_),
		os1(os1_), os2(os2_),
		p1(p1_), p2(p2_), s1(s1_), s2(s2_), 
		companion(3,3), R(R_), eigvalues(3),
		bounds_pk1(2), bounds_pk2(2), EXPO(1.0)//,
//		ecompanion(3,3)
{
	companion.zeros();
	companion(1,0) = 1.0;
	companion(2,1) = 1.0;
// 	ecompanion(0,0) = 0.0; ecompanion(0,1) = 0.0;
//	ecompanion(1,0) = 1.0; ecompanion(1,1) = 0.0;
//	ecompanion(2,0) = 0.0; ecompanion(2,1) = 1.0;
}

std::map<std::string, double> RV_Missing_t_walk::state() const {
	std::map<std::string, double> out;
	out["x1"] = x1;
	out["X"] = x2;
	out["x3"] = x3;
	out["os1"] = os1;
	out["os2"] = os2;
	out["p1"] = p1;
	out["p2"] = p2;
	out["s1"] = s1;
	out["s2"] = s2;
	return out;
}

void RV_Missing_t_walk::jump(double X) { x2 = X; }

double RV_Missing_t_walk::draw() {  // Slice sampler w/ two peaks.
  ly = lpdf() - EXPO(R);
	find_slice();
	double ii = 0;
	while(true) {
		//ii++; // if (ii > 10000) throw std::runtime_error("Too many steps.");
		ii++; if (ii > 10000) { std::cout << "Too many steps." << std::endl; }
		x_new = choose();	
		ly_new = lpdf(x_new);
		if (ly_new >= ly) // && (x_new > min) && (x_new < max))  truncation.
			break; 
		else 
			trim();
	}
	x2 = x_new;
	return x2; 
}

double RV_Missing_t_walk::lpdf(double X) {
	double lpdf;
	lpdf = 	boost::math::lgamma((p1+1.0)/2.0) - boost::math::lgamma(p1/2.0) -
		0.5 * log(p1*pi*pow(s1,2)) - 
		(p1+1.0)/2.0 * log(1.0 + (pow(X -(x1+os1),2))/(p1*pow(s1,2)) ) +
					boost::math::lgamma((p2+1.0)/2.0) - boost::math::lgamma(p2/2.0) -
		0.5 * log(p2*pi*pow(s2,2)) - 
		(p2+1.0)/2.0 * log(1.0 + (pow(x3-(X +os2),2))/(p2*pow(s2,2)) );
	return lpdf;
}

double RV_Missing_t_walk::lpdf() { return lpdf(x2); }

void RV_Missing_t_walk::find_slice() {
	find_peaks();
	if (ly <= lpdf(peak1)) bounds_pk1 = step_out(peak1);
	if (ly <= lpdf(peak2)) bounds_pk2 = step_out(peak2);
	if (ly > lpdf(peak1)) {   // In case peak 1 does not reach slice level.
		bounds_pk1[0] = bounds_pk2[0];
		bounds_pk1[1] = bounds_pk2[0];
	}
	if (ly > lpdf(peak2)) {		// In case peak 2 does not reach slice level.
		bounds_pk2[0] = bounds_pk1[0];
		bounds_pk2[1] = bounds_pk1[0];
	}
	if (  // In case the two parts of the slice overlap, split it:
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
	int m = 200;
	double w = (s1+s2)/2.0;
	bounds[0] = peak - w * U(R);  // w = use (s1+s2)/2
	bounds[1] = bounds[0] + w;
	int j = std::floor(m * U(R));    // m needed...
	int k = (m-1) - j;
	while ((j>0) && (ly < lpdf(bounds[0])) ) { // trunc. can be added here.
		bounds[0] = bounds[0] - w;
		j = j - 1;
	}
	while ((k>0) && (ly < lpdf(bounds[1])) ) { // truncation can be added here.
		bounds[1] = bounds[1] + w;
		k = k - 1;
	}
	return bounds;
}

void RV_Missing_t_walk::trim() {
	if (x_new <= bounds_pk1[1]) {
		if (x_new <= peak1) {
			bounds_pk1[0] = x_new;
		} else {
			bounds_pk1[1] = x_new;
		}
	} else {
		if (x_new <= peak2) {
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
	companion(0,2) = 0.5 * (
										(x1+os1)*(x1+os1)*(x3-os2) + 
										(x1+os1)*p2*s2*s2 + 
										(x3-os2)*p1*s1*s1 + 
										(x3-os2)*(x3-os2)*(x1+os1));
	companion(1,2) = -0.5 * (
										p1*s1*s1 + 
										(x1+os1)*(x1+os1) + 
										4.0*(x1+os1)*(x3-os2) + 
										(x3-os2)*(x3-os2) + p2*s2*s2);
	companion(2,2) = 1.5 * ((x1+os1)+(x3-os2));
	if (!arma::eig_gen(cx_eigval, cx_eigvec, companion)) 
		throw std::runtime_error("Failed eigenvalue decomposition.");
	eigvalues = arma::real(cx_eigval);
	eigvalues = arma::sort(eigvalues);
	// Calculate condition number, if it's < 10^(-6) or so, 
	// set both peaks to the average value... 
	peak1 = eigvalues[0];
	peak2 = eigvalues[2];
//	ecompanion(0,2) = companion(0,2);
//	ecompanion(1,2) = companion(1,2);
//	ecompanion(2,2) = companion(2,2);
//	Eigen::EigenSolver<Eigen::MatrixXd> es(ecompanion);
//	std::cout << "Peaks: \n" << es.eigenvalues() << std::endl;
	
}

// observed_normal after here.
RV_Missing_t_walk_observed_normal::RV_Missing_t_walk_observed_normal(
		double const & x1_,
		double 			 & X,
		double const & x3_,
		double const & os1_,
		double const & os2_,
		double const & p1_,   // degrees of freedom 1
		double const & p2_,   // degrees of freedom 2
		double const & s1_,   // scale 1
		double const & s2_,   // scale 2
		double const & Xobs_, // observed location
		double const & so2_,
		trng::yarn2 & R_
) :	x1(x1_), x2(X), x3(x3_),
		os1(os1_), os2(os2_),
		p1(p1_), p2(p2_), s1(s1_), s2(s2_), 
		Xobs(Xobs_), so2(so2_),
		companion(3,3), R(R_), eigvalues(3),
		bounds_pk1(2), bounds_pk2(2), EXPO(1.0)//,
//		ecompanion(3,3)
{
	companion.zeros();
	companion(1,0) = 1.0;
	companion(2,1) = 1.0;
// 	ecompanion(0,0) = 0.0; ecompanion(0,1) = 0.0;
//	ecompanion(1,0) = 1.0; ecompanion(1,1) = 0.0;
//	ecompanion(2,0) = 0.0; ecompanion(2,1) = 1.0;
}

std::map<std::string, double> RV_Missing_t_walk_observed_normal::state() const {
	std::map<std::string, double> out;
	out["x1"] = x1;
	out["X"] = x2;
	out["x3"] = x3;
	out["os1"] = os1;
	out["os2"] = os2;
	out["p1"] = p1;
	out["p2"] = p2;
	out["s1"] = s1;
	out["s2"] = s2;
	out["so2"] = so2;
	out["Xobs"] = Xobs;
	return out;
}

void RV_Missing_t_walk_observed_normal::jump(double X) { x2 = X; }

double RV_Missing_t_walk_observed_normal::draw() {  // Slice sampler w/ two peaks.
  ly = lpdf() - EXPO(R);
	find_slice();
	double ii = 0;
	while(true) {
		//ii++; // if (ii > 10000) throw std::runtime_error("Too many steps.");
		ii++; if (ii > 10000) { std::cout << "Too many steps." << std::endl; }
		x_new = choose();	
		ly_new = lpdf(x_new);
		if (ly_new >= ly) // && (x_new > min) && (x_new < max))  truncation.
			break; 
		else 
			trim();
	}
	x2 = x_new;
	return x2; 
}

double RV_Missing_t_walk_observed_normal::lpdf(double X) {
	double lpdf;
	lpdf = 	boost::math::lgamma((p1+1.0)/2.0) - boost::math::lgamma(p1/2.0) -
		0.5 * log(p1*pi*pow(s1,2)) - 
		(p1+1.0)/2.0 * log(1.0 + (pow(X -(x1+os1),2))/(p1*pow(s1,2)) ) +
					boost::math::lgamma((p2+1.0)/2.0) - boost::math::lgamma(p2/2.0) -
		0.5 * log(p2*pi*pow(s2,2)) - 
		(p2+1.0)/2.0 * log(1.0 + (pow(x3-(X +os2),2))/(p2*pow(s2,2)) );
	lpdf += -0.5 * log(2*PI*pow(so2,2)) - (pow(X - Xobs,2) / (2*pow(so2,2)));
	return lpdf;
}

double RV_Missing_t_walk_observed_normal::lpdf() { return lpdf(x2); }

void RV_Missing_t_walk_observed_normal::find_slice() {
	find_peaks();
	if (ly <= lpdf(peak1)) bounds_pk1 = step_out(peak1);
	if (ly <= lpdf(peak2)) bounds_pk2 = step_out(peak2);
	if (ly > lpdf(peak1)) {   // In case peak 1 does not reach slice level.
		bounds_pk1[0] = bounds_pk2[0];
		bounds_pk1[1] = bounds_pk2[0];
	}
	if (ly > lpdf(peak2)) {		// In case peak 2 does not reach slice level.
		bounds_pk2[0] = bounds_pk1[0];
		bounds_pk2[1] = bounds_pk1[0];
	}
	if (  // In case the two parts of the slice overlap, split it:
		(bounds_pk1[0] < std::min(bounds_pk1[1], bounds_pk2[1])) &&
		(bounds_pk2[0] < std::min(bounds_pk1[1], bounds_pk2[1]))
	) {
		double x = (bounds_pk1[0]+bounds_pk2[1])/2.0;
		bounds_pk1[1] = x;
		bounds_pk2[0] = x;
	}
}


std::vector<double> RV_Missing_t_walk_observed_normal::step_out(double peak) {
	std::vector<double> bounds(2);
	int m = 200;
	double w = (s1+s2)/2.0;
	bounds[0] = peak - w * U(R);  // w = use (s1+s2)/2
	bounds[1] = bounds[0] + w;
	int j = std::floor(m * U(R));    // m needed...
	int k = (m-1) - j;
	while ((j>0) && (ly < lpdf(bounds[0])) ) { // trunc. can be added here.
		bounds[0] = bounds[0] - w;
		j = j - 1;
	}
	while ((k>0) && (ly < lpdf(bounds[1])) ) { // truncation can be added here.
		bounds[1] = bounds[1] + w;
		k = k - 1;
	}
	return bounds;
}

void RV_Missing_t_walk_observed_normal::trim() {
	if (x_new <= bounds_pk1[1]) {
		if (x_new <= peak1) {
			bounds_pk1[0] = x_new;
		} else {
			bounds_pk1[1] = x_new;
		}
	} else {
		if (x_new <= peak2) {
			bounds_pk2[0] = x_new;
		} else {
			bounds_pk2[1] = x_new;
		}
	}
}

double RV_Missing_t_walk_observed_normal::choose() {
	double d = (bounds_pk2[1] - bounds_pk2[0]) + (bounds_pk1[1] - bounds_pk1[0]);
	double q = (bounds_pk2[0] - bounds_pk1[1]);
	double l = U(R) * d;
	if ((bounds_pk1[0] + l) > bounds_pk1[1]) {
		return (l + q + bounds_pk1[0]);
	} else {
		return (l + bounds_pk1[0]);
	}
}

void RV_Missing_t_walk_observed_normal::find_peaks() {
	// NEEDS MODIFICATION TO INCLUDE OBSERVATION PROCESS
	companion(0,2) = 0.5 * (
										(x1+os1)*(x1+os1)*(x3-os2) + 
										(x1+os1)*p2*s2*s2 + 
										(x3-os2)*p1*s1*s1 + 
										(x3-os2)*(x3-os2)*(x1+os1));
	companion(1,2) = -0.5 * (
										p1*s1*s1 + 
										(x1+os1)*(x1+os1) + 
										4.0*(x1+os1)*(x3-os2) + 
										(x3-os2)*(x3-os2) + p2*s2*s2);
	companion(2,2) = 1.5 * ((x1+os1)+(x3-os2));
	if (!arma::eig_gen(cx_eigval, cx_eigvec, companion)) 
		throw std::runtime_error("Failed eigenvalue decomposition.");
	eigvalues = arma::real(cx_eigval);
	eigvalues = arma::sort(eigvalues);
	// Calculate condition number, if it's < 10^(-6) or so, 
	// set both peaks to the average value... 
	peak1 = eigvalues[0];
	peak2 = eigvalues[2];
//	ecompanion(0,2) = companion(0,2);
//	ecompanion(1,2) = companion(1,2);
//	ecompanion(2,2) = companion(2,2);
//	Eigen::EigenSolver<Eigen::MatrixXd> es(ecompanion);
//	std::cout << "Peaks: \n" << es.eigenvalues() << std::endl;
}	

