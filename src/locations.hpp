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

class Locations {

public:
	Locations(
			arma::vec & locations_, arma::vec & drift_, 
			arma::vec & tails_, arma::vec & scales_, 
			arma::vec & obs_scales_,
			arma::vec & minima_, arma::vec & maxima_, arma::vec & draws_, trng::yarn2 & R_);
	arma::vec & state() const;
	double & state(arma::uword which) const;

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
	arma::vec & locations;
	arma::vec & drift;
	arma::vec & tails;
	arma::vec & scales;
	arma::vec & obs_scales;
	arma::vec & minima;
	arma::vec & maxima;

	std::vector<std::unique_ptr<Random> > distributions;
	std::map<unsigned int, std::vector<unsigned int> > sample_order;
	std::map<unsigned int, std::vector<unsigned int> > sample_order_bk;
	arma::vec & draws;

};

#endif
