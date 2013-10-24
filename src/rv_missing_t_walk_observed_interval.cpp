#include "rv_missing_t_walk_observed_interval.hpp"

#include <boost/math/special_functions/gamma.hpp>
#include <math.h>
#include <limits>

const double pi = boost::math::constants::pi<double>();


RV_Missing_t_walk_observed_interval::RV_Missing_t_walk_observed_interval(
		double const & x1_,
		double 			 & X,
		double const & x3_,
		double const & os1_,
		double const & os2_,
		double const & p1_,   // degrees of freedom 1
		double const & p2_,   // degrees of freedom 2
		double const & s1_,   // scale 1
		double const & s2_,   // scale 2
		double const & Xmin_,
		double const & Xmax_,
		trng::yarn2 & R_
) :	RV_Missing_t_walk_core(x1_, X, x3_, os1_, os2_, p1_, p2_, s1_, s2_, R_),
		Xmin(Xmin_), Xmax(Xmax_)
{
	std::cout << "Xmin: " << Xmin << ", Xmax: " << Xmax << std::endl;

	// Companion matrix for eigenvalue peak-finding.
	companion.set_size(3,3);
	companion.zeros();
	companion(1,0) = 1.0;
	companion(2,1) = 1.0;

	// Reserve sizes, is this important?
	peaks.reserve(2);
	peak_bound_lr.reserve(2);
	intervals.reserve(1);
}

std::map<std::string, double> RV_Missing_t_walk_observed_interval::state() const {
	std::map<std::string, double> state = RV_Missing_t_walk_core::state();
	state["Xmin"] = Xmin;
	state["Xmax"] = Xmax;
	return state;
}


double RV_Missing_t_walk_observed_interval::lpdf(double X) {
	std::cout << "lpdf of RV_Missing_t_walk_observed_normal" << std::endl;
	std::cout << "x1: " << x1 << std::endl;
	std::cout << "X:  " << X  << std::endl;
	std::cout << "x3: " << x3 << std::endl;
	std::cout << "os1: " << os1 << std::endl;
	std::cout << "os2: " << os2 << std::endl;
	std::cout << "p1: " << p1 << std::endl;
	std::cout << "p2: " << p2 << std::endl;
	std::cout << "s1: " << s1 << std::endl;
	std::cout << "s2: " << s2 << std::endl;
	std::cout << "Xmin: " << Xmin << ", Xmax: " << Xmax << std::endl;
	if ((X < Xmin) || (X > Xmax)) {		
		return -1.0 * std::numeric_limits<double>::infinity();
	}
	double lpdf;
	lpdf = 	boost::math::lgamma((p1+1.0)/2.0) - boost::math::lgamma(p1/2.0) -
		0.5 * log(p1*pi*pow(s1,2)) - 
		(p1+1.0)/2.0 * log(1.0 + (pow(X -(x1+os1),2))/(p1*pow(s1,2)) ) +
					boost::math::lgamma((p2+1.0)/2.0) - boost::math::lgamma(p2/2.0) -
		0.5 * log(p2*pi*pow(s2,2)) - 
		(p2+1.0)/2.0 * log(1.0 + (pow(x3-(X +os2),2))/(p2*pow(s2,2)) );
	return lpdf;
}

double RV_Missing_t_walk_observed_interval::lpdf() { return lpdf(x2); }


void RV_Missing_t_walk_observed_interval::derivative_poly() {
	double B1 = ( -1*(p1+1) - (p2+1) );
  double B2 = ( (p1+1)*(x1+os1+2*x3-2*os2) +
          (p2+1)*(x3-os2+2*x1+2*os1)  );
  double B3 = (
	  -1*(p1+1)*(
		p2*pow(s2,2) + pow(x3-os2,2) + 
			2*(x1+os1)*(x3-os2)
		) -
		 1*(p2+1)*(
		p1*pow(s1,2) + pow(x1+os1,2) + 
			2*(x3-os2)*(x1+os1)) 
	);
  double B4 = ( 
		(p1+1)*(x1+os1)*
			(p2*pow(s2,2)+pow(x3-os2,2)) +
    (p2+1)*(x3-os2)*
			(p1*pow(s1,2)+pow(x1+os1,2)) 
	);

  companion(0,2) = -B4/B1;
	companion(1,2) = -B3/B1;
	companion(2,2) = -B2/B1;
}	



