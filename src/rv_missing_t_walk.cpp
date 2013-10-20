#include "rv_missing_t_walk.hpp"

#include <boost/math/special_functions/gamma.hpp>
#include <math.h>
#include <limits>

const double pi = boost::math::constants::pi<double>();



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
{
	companion.zeros();
	companion(1,0) = 1.0;
	companion(2,1) = 1.0;
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
	double B1 = ( -1*(ptm1+1) - (pt__+1) );
  double B2 = ( (ptm1+1)*(xtm1+deltatm1+2*xtp1-2*deltat__) +
          (pt__+1)*(xtp1-deltat__+2*xtm1+2*deltatm1)  );
  double B3 = (
	  -1*(ptm1+1)*(
		pt__*pow(sigmat__,2) + pow(xtp1-deltat__,2) + 
			2*(xtm1+deltatm1)*(xtp1-deltat__)
		) -
		 1*(pt__+1)*(
		ptm1*pow(sigmatm1,2) + pow(xtm1+deltatm1,2) + 
			2*(xtp1-deltat__)*(xtm1+deltatm1)) 
	);
  double B4 = ( 
		(ptm1+1)*(xtm1+deltatm1)*
			(pt__*pow(sigmat__,2)+pow(xtp1-deltat__,2)) +
    (pt__+1)*(xtp1-deltat__)*
			(ptm1*pow(sigmatm1,2)+pow(xtm1+deltatm1,2)) 
	);

  companion(0,2) = -B4/B1;
	companion(1,2) = -B3/B1;
	companion(2,2) = -B2/B1;
	
	if (!arma::eig_gen(cx_eigval, cx_eigvec, companion)) 
		throw std::runtime_error("Failed eigenvalue decomposition.");
	arma::cx_vec::iterator re_eigval_end = 
		std::remove_if(cx_eigval.begin(), cx_eigval.end(), has_imaginary);
	std:sort(cx_eigval.begin(), re_eigval_end);
	num_real_eigval = re_eigval_end - cx_eigval.begin() - 1;
	if (num_real_eigval == 0) {
		throw std::runtime_error("No real eigenvalues available to find peaks.");
	} else if (num_real_eigval == 1) {
		peak1 = arma::real(cx_eigval)[0];
		peak2 = arma::real(cx_eigval)[0];
	} else if (num_real_eigval == 3) {
		peak1 = arma::real(cx_eigval)[0];
		peak2 = arma::real(cx_eigval)[2];
	} else {
		throw std::runtime_error("Number of real eigenvalues not equal to 0, 1, or 3.");
	}
	
}

