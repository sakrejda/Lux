#ifndef RV_MISSING_T_WALK_OBSERVED_INTERVAL_H
#define RV_MISSING_T_WALK_OBSERVED_INTERVAL_H

#include "rv_missing_t_walk_core.hpp"

class RV_Missing_t_walk_observed_interval : public RV_Missing_t_walk_core {

public:
	RV_Missing_t_walk_observed_interval(
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
	);

	std::map<std::string, double> state() const;
	double lpdf(double X);
	double lpdf();

private:
	void derivative_poly();

	double const & Xmin;
	double const & Xmax;

};

#endif
