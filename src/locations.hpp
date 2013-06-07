#ifndef LOCATION_H
#define LOCATION_H

#include "random.hpp"

#include <vector>
#include <map>
#include <memory>

#include <armadillo>
#include <trng/yarn2.hpp>

class Locations {

public:
	Locations(arma::vec & locations_, arma::vec & tails_, arma::vec & scales_, 
			trng::yarn2 & R_);
	arma::vec & state() const;
	double & state(arma::uword which) const;

	// Available distributions:
	void bind_constant_distribution	(unsigned int which);
	void bind_uniform_distribution	(
			unsigned int which, 
			double const & minimum,
			double const & maximum, trng::yarn2 & R);
	void bind_uniform_distribution (
			unsigned int which, trng::yarn2 & R);
	void bind_t_walk_distribution		(
			unsigned int which, 
			double const & p1, double const & p2, 
			double const & s1, double const & s2, trng::yarn2 & R);
	void bind_t_walk_distribution (
			unsigned int which, trng::yarn2 & R);

	~Locations();  

private:
	trng::yarn2 & R;
	arma::vec & locations;
	arma::vec & tails;
	arma::vec & scales;
	std::vector<std::unique_ptr<Random> > distributions;

};

#endif
