#ifndef RV_MISSING_T_WALK_OBSERVED_NORMAL_H
#define RV_MISSING_T_WALK_OBSERVED_NORMAL_H

#include "random.hpp"

#include <vector>
#include <map> 

#include <armadillo>
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>
#include <trng/exponential_dist.hpp>


class RV_Missing_t_walk_observed_normal : public Random {

public:
	RV_Missing_t_walk_observed_normal(
		double const & x1_,
		double 			 & X,
		double const & x3_,
		double const & os1_,
		double const & os2_,
		double const & p1_,   // degrees of freedom 1
		double const & p2_,   // degrees of freedom 2
		double const & s1_,   // scale 1
		double const & s2_,   // scale 2
		double const & Xobs_,
		double const & so2_,
		trng::yarn2 & R_
	);

	std::map<std::string, double> state() const;
	void jump(double X);
	double draw();
	double lpdf(double X);
	double lpdf();

private:
	void find_peaks();
	void find_slice();
	std::vector<double> step_out(double peak);
	double choose();
	void trim();

	double const & x1; 
	double       & x2;
	double const & x3;
	double const & os1;
	double const & os2;
	double const & p1;
	double const & p2;
	double const & s1;
	double const & s2;
	double const & Xobs;
	double const & so2;

	double ly; // cached log pdf at current value...
	double x_new; // temporary new sampled value...
	double ly_new; // temporary pdf at sampled value...

	trng::yarn2 & R;  
	trng::uniform01_dist<double> U;
	trng::exponential_dist<double> EXPO;

	double A1, A2, A3, A4, B1, B2, B3, B4;
	arma::Mat<double> companion;
	arma::cx_vec cx_eigval;
	arma::cx_mat cx_eigvec;
	void derivative_poly();

	std::vector<double> peaks;
	std::vector<std::vector<double> > peak_bounds_lr;
	std::vector<double> intervals;
  double total_slice_length;

};

#endif
