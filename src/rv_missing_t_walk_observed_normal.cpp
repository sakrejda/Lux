#include "rv_missing_t_walk_observed_normal.hpp"

#include <boost/math/special_functions/gamma.hpp>
#include <math.h>
#include <limits>

#include <iostream>

const double pi = boost::math::constants::pi<double>();


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
) :	RV_Missing_t_walk_core(x1_, X, x3_, os1_, os2_, p1_, p2_, s1_, s2_, R_),
		Xobs(Xobs_), so2(so2_)
{
	// Companion matrix for eigenvalue peak-finding.
	companion.set_size(5,5);
	companion.zeros();
	companion(1,0) = 1.0;
	companion(2,1) = 1.0;
	companion(3,2) = 1.0;
	companion(4,3) = 1.0;

	// Reserve sizes, is this important?
	peaks.reserve(3);
	peak_bound_lr.reserve(3);
	intervals.reserve(2);
}

std::map<std::string, double> RV_Missing_t_walk_observed_normal::state() const {
	std::map<std::string, double> state = RV_Missing_t_walk_core::state();
	state["so2"] = so2;
	state["Xobs"] = Xobs;
	return state;
}


double RV_Missing_t_walk_observed_normal::lpdf(double X) {
	double lpdf;
	lpdf = 	boost::math::lgamma((p1+1.0)/2.0) - boost::math::lgamma(p1/2.0) -
		0.5 * log(p1*pi*pow(s1,2)) - 
		(p1+1.0)/2.0 * log(1.0 + (pow(X -(x1+os1),2))/(p1*pow(s1,2)) ) +
					boost::math::lgamma((p2+1.0)/2.0) - boost::math::lgamma(p2/2.0) -
		0.5 * log(p2*pi*pow(s2,2)) - 
		(p2+1.0)/2.0 * log(1.0 + (pow(x3-(X +os2),2))/(p2*pow(s2,2)) );
	lpdf += -0.5 * log(2*pi*pow(so2,2)) - (pow(X - Xobs,2) / (2*pow(so2,2)));
	return lpdf;
}

double RV_Missing_t_walk_observed_normal::lpdf() { return lpdf(x2); }


void RV_Missing_t_walk_observed_normal::derivative_poly() {
	std::cout << "In poly: " << std::endl;
  std::cout << "x1: " << x1 << ", x2: " << x2 << ", x3: " << x3 << std::endl;
	std::cout << "os1: " << os1 << ", os2: " << os2 << std::endl;
  std::cout << "s1: " << s1 << ", s2: " << s2 << std::endl;
  std::cout << "p1: " << p1 << ", p2: " << p2 << std::endl;
	std::cout << "Xobs: " << Xobs << "so2: " << so2 << std::endl;

	double A1 = 2*( (os2-x3) - (os1+x1) );
  double A2 = ( 
		pow(x3-os2,2) + p2*pow(s2,2) +
    pow(x1+os1,2) + p1*pow(s1,2) +
   -4*(os1+x1)*(os2-x3) 
	);
  double A3 = (
		2*( (os2-x3)*(pow(x1+os1,2)+p1*pow(s1,2)) - 
        (os1+x1)*(pow(x3-os2,2)+p2*pow(s2,2)) 
		)
	);
  double A4 = ( 
		(pow(x3-os2,2) + p2*pow(s2,2)) * 
    (pow(x1+os1,2) + p1*pow(s1,2)) 
	);

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

	companion(0,4) = A4*Xobs           + B4*so2;
	companion(1,4) = A3*Xobs - A4      + B3*so2;
	companion(2,4) = A2*Xobs - A3      + B2*so2;
	companion(3,4) = A1*Xobs - A2      + B1*so2;
	companion(4,4) =    Xobs - A1							 ;
	std::cout << companion << std::endl;
}


