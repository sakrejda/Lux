
#include "locations.hpp"

#include <armadillo>
#include <trng/yarn2.hpp>

Locations::Locations(
		arma::vec & locations_, arma::vec & tails_, arma::vec & scales_, 
		arma::vec & minima_, arma::vec & maxima_, trng::yarn2 & R_
) : locations(locations_), 
		tails(tails_),
		scales(scales_),
		minima(minima_),
		maxima(maxima_),
		R(R_),
		distributions(locations_.size()) { }

arma::vec & Locations::state() const { return locations; }
double & Locations::state(arma::uword which) const { return locations[which]; }

// Available distributions:
void Locations::bind_constant_distribution	(
		unsigned int which
) {
	distributions[which] = 
		std::unique_ptr<Random>(new RV_Constant(locations[which]));
}


void Locations::bind_uniform_distribution (
		unsigned int which, trng::yarn2 & R
) {
	distributions[which] = 
		std::unique_ptr<Random>(
				new RV_Uniform(locations[which], minima[which], maxima[which], R));
}

void Locations::bind_ordered_uniform_distribution (
		unsigned int which, 
		trng::yarn2 & R
) {
	distributions[which] = 
		std::unique_ptr<Random>(new RV_Uniform(
				locations[which], locations[which-1], locations[which+1], R));
}

void Locations::bind_t_walk_distribution (
		unsigned int which,
		double const & p1, double const & p2,
		double const & s1, double const & s2, trng::yarn2 & R
) {
	distributions[which] = 
		std::unique_ptr<Random>(new RV_Missing_t_walk(
			locations[which-1], locations[which], locations[which+1], 
			p1, p2, s1, s2, R));
}

void Locations::bind_t_walk_distribution (
		unsigned int which, trng::yarn2 & R
) {
	distributions[which] = 
		std::unique_ptr<Random>(new RV_Missing_t_walk(
			locations[which-1], locations[which], locations[which+1], 
			tails[which-1], tails[which],
			scales[which-1], scales[which], R));
}


Locations::~Locations() { distributions.clear(); }
