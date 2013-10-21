#include "rv_missing_t_walk_observed_normal.hpp"

#include <boost/math/special_functions/gamma.hpp>
#include <math.h>
#include <limits>

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
) :	RV_Missing_t_walk_observed_core(x1_, X, x3_, os1_, os2_, p1_, p2_, s1_, s2_, R_),
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
	double A1 = 2*( (deltat__-xtp1) - (deltatm1+xtm1) );
  double A2 = ( 
		pow(xtp1-deltat__,2) + pt__*pow(sigmat__,2) +
    pow(xtm1+deltatm1,2) + ptm1*pow(sigmatm1,2) +
   -4*(deltatm1+xtm1)*(deltat__-xtp1) 
	);
  double A3 = (
		2*( (deltat__-xtp1)*(pow(xtm1+deltatm1,2)+ptm1*pow(sigmatm1,2)) - 
        (deltatm1+xtm1)*(pow(xtp1-deltat__,2)+pt__*pow(sigmat__,2)) 
		)
	);
  double A4 = ( 
		(pow(xtp1-deltat__,2) + pt__*pow(sigmat__,2)) * 
    (pow(xtm1+deltatm1,2) + ptm1*pow(sigmatm1,2)) 
	);

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

	companion(0,4) = A4*Xobs           + B4*so2;
	companion(1,4) = A3*Xobs - A4      + B3*so2;
	companion(2,4) = A2*Xobs - A3      + B2*so2;
	companion(3,4) = A1*Xobs - A2      + B1*so2;
	companion(4,4) =    Xobs - A1							 ;
}


