#ifndef RV_MISSING_T_WALK_H
#define RV_MISSING_T_WALK_H

#include "random.hpp"

#include <vector>
#include <map> 

#include <armadillo>
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>
#include <trng/exponential_dist.hpp>


class RV_Missing_t_walk : public Random {

public:
	typedef double (*PDF)(double);

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

	std::map<std::string, double> state() const;
	void jump(double X);
	double draw();
	double lpdf(double X);
	double lpdf();

private:
	void find_peaks();
	double peak1, peak2;
	arma::Mat<double> companion;
	arma::cx_vec cx_eigval;
	arma::vec    eigvalues;
	arma::cx_mat cx_eigvec;


	void print_slice(std::string s);

	void find_slice();
	std::vector<double> bounds_pk1;
	std::vector<double> bounds_pk2;
	std::vector<double> step_out(double peak);
	double choose();
	void trim();
	double ly; // cached log pdf at current value...
	double x_new; // temporary new sampled value...
	double ly_new; // temporary pdf at sampled value...

	trng::yarn2 & R;  
	trng::uniform01_dist<double> U;
	trng::exponential_dist<double> EXPO;

	double const & x1; 
	double       & x2;
	double const & x3;
	double const & os1;
	double const & os2;
	double const & p1;
	double const & p2;
	double const & s1;
	double const & s2;

};


#endif
