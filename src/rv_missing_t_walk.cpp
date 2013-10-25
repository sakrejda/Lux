#include "rv_missing_t_walk.hpp"

#include <limits>

RV_Missing_t_walk::NEG_INF = -1.0 * std::numeric_limits<double>::infinity(),
RV_Missing_t_walk::POS_INF =        std::numeric_limits<double>::infinity(), 

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
) :	RV_Missing_t_walk_observed_interval(
			x1_, X, x3_, os1_, os2_, p1_, p2_, s1_, s2_, 
			NEG_INF, POS_INF, R(R_) { }
