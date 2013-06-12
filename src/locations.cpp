
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
	if (distributions[which] == NULL) {
		distributions[which] = 
			std::unique_ptr<Random>(new RV_Constant(locations[which]));
	} else {
		std::stringstream msg;
		msg << "The location " << which << " (" << (which+1) << ")"
					 " already has a distribution.  Not adding.\n";
		throw(std::logic_error(msg.str()));
	}
}


void Locations::bind_uniform_distribution (
		unsigned int which, trng::yarn2 & R
) {
	if (distributions[which] == NULL) {
		distributions[which] = 
			std::unique_ptr<Random>(
				new RV_Uniform(locations[which], minima[which], maxima[which], R));
	} else {
		std::stringstream msg;
		msg << "The location " << which << " (" << (which+1) << ")"
					 " already has a distribution.  Not adding.\n";
		throw(std::logic_error(msg.str()));
	}
}

void Locations::bind_ordered_uniform_distribution (
		unsigned int which, 
		trng::yarn2 & R
) {
	if (distributions[which] == NULL) {
		distributions[which] = 
			std::unique_ptr<Random>(new RV_Uniform(
				locations[which], locations[which-1], locations[which+1], R));
	} else {
		std::stringstream msg;
		msg << "The location " << which << " (" << (which+1) << ")"
					 " already has a distribution.  Not adding.\n";
		throw(std::logic_error(msg.str()));
	}
}

void Locations::bind_t_walk_distribution (
		unsigned int which,
		double const & p1, double const & p2,
		double const & s1, double const & s2, trng::yarn2 & R
) {
	if (distributions[which] == NULL) {
		distributions[which] = 
			std::unique_ptr<Random>(new RV_Missing_t_walk(
				locations[which-1], locations[which], locations[which+1], 
			p1, p2, s1, s2, R));
	} else {
		std::stringstream msg;
		msg << "The location " << which << " (" << (which+1) << ")"
					 " already has a distribution.  Not adding.\n";
		throw(std::logic_error(msg.str()));
	}
}

void Locations::bind_t_walk_distribution (
		unsigned int which, trng::yarn2 & R
) {
	if (distributions[which] == NULL) {
		distributions[which] = 
			std::unique_ptr<Random>(new RV_Missing_t_walk(
				locations[which-1], locations[which], locations[which+1], 
				tails[which-1], tails[which],
				scales[which-1], scales[which], R));
	} else {
		std::stringstream msg;
		msg << "The location " << which << " (" << (which+1) << ")"
					 " already has a distribution.  Not adding.\n";
		throw(std::logic_error(msg.str()));
	}
}

void Locations::drop_distribution(unsigned int which) {
	if (distributions[which] == NULL) {
		std::stringstream msg;
		msg << "The location " << which << " (" << (which+1) << ")"
					 " does not have a distribution.  Not deleting.\n";
		throw(std::logic_error(msg.str()));
	} else {
		distributions[which].reset(NULL);
	}	
	
}

Locations::~Locations() { distributions.clear(); }
