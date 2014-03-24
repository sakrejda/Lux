#ifndef LOCATION_H
#define LOCATION_H

#include <vector>
#include <map>
#include <memory>
#include <armadillo>
#include <trng/yarn2.hpp>

// This would lend itself to cleaner user code if the Locations class
// took ownership of the data---somebody has to and you don't want to
// hack that at the R/C++ interface (so I've learned)... everything else
// could be the same, except the Rcpp/R interface code would have to act
// like it doesn't have ownership...

class Random;

class Time_Series_State { // State in that observed and unobserved sense.

public:
	Time_Series_State(); // required if used by reference...
	Time_Series_State(
		std::vector<double> times_,
		std::vector<double> y_at_times_,
		std::vector<double> minima_at_times_,
		std::vector<double> maxima_at_times_
	);

const arma::Col<double> & get_times() const;
const arma::Col<double> & get_x_at_times() const;
const arma::Col<double> & get_y_at_times() const;
const arma::Col<double> & get_minima_at_times() const;
const arma::Col<double> & get_maxima_at_times() const;

void set_x_at_times(arma::Col<double> x_at_times_);

private:
	arma::Col<double> times;
	arma::Col<double> x_at_times;
	arma::Col<double> y_at_times;
	arma::Col<double> minima_at_times;
	arma::Col<double> maxima_at_times;

};

class Time_Series_Parameters { // Parameters in that unobservable sense... (?) 

public:
	Time_Series_Parameters();
	Time_Series_Parameters(
			Time_Series_State const & state_,
	);

	const arma::Col<double> get_drift;
	const arma::Col<double> get_scales;
	const arma::Col<double> get_tails;
	const arma::Col<double> get_obs_scales;

	void set_drift(arma::Col<double> drift_);
	void set_scales(arma::Col<double> scales_);
	void set_tails(arma::Col<double> tails_);
	void set_obs_scales(arma::Col<double> obs_scales_);

private:
	Time_Series_State const & state;
	arma::Col<double> drift;
	arma::Col<double> scales;
	arma::Col<double> tails;
	arma::Col<double> obs_scales;

};


class Locations {

public:
	Locations(
			arma::vec locations_, arma::vec drift_, 
			arma::vec tails_, arma::vec scales_, 
			arma::vec obs_scales_,
			arma::vec minima_, arma::vec maxima_, arma::vec draws_, trng::yarn2 & R_);
	const arma::vec & state() const;
	const double & state(arma::uword which) const;

	// Available distributions:
	void bind_constant_distribution	(unsigned int which);
	void bind_uniform_distribution	(
			unsigned int which, trng::yarn2 & R);
	void bind_ordered_uniform_distribution (
			unsigned int which, trng::yarn2 & R);
	void bind_normal_distribution (
			unsigned int which, trng::yarn2 & R);
	void bind_t_walk_distribution_open(
			unsigned int which, trng::yarn2 & R);
	void bind_t_walk_distribution_open_reverse(
			unsigned int which, trng::yarn2 & R);
	void bind_t_walk_distribution (
			unsigned int which, trng::yarn2 & R);
	void bind_t_walk_observed_normal_distribution (
			unsigned int which, trng::yarn2 & R);
	void bind_t_walk_observed_interval_distribution (
			unsigned int which, trng::yarn2 & R);

	// Drop distribution:
	void drop_distribution(unsigned int which);

	void draw();
	arma::vec lpdf(arma::vec X);

	~Locations();  

private:
	trng::yarn2 & R;

	std::vector<std::unique_ptr<Random> > distributions;
	std::map<unsigned int, std::vector<unsigned int> > sample_order;
	std::map<unsigned int, std::vector<unsigned int> > sample_order_bk;

};

#endif
