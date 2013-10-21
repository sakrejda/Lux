#ifndef RV_MISSING_T_WALK_H
#define RV_MISSING_T_WALK_H

#include "rv_missing_t_walk_observed_interval.hpp"

class RV_Missing_t_walk : public RV_Missing_t_walk_observed_interval {

public:
	RV_Missing_t_walk(
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
	);

};


#endif
