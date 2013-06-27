
#include "locations.hpp"

#include <armadillo>
#include <trng/yarn2.hpp>

Locations::Locations(
		arma::vec & locations_, arma::vec & tails_,  arma::vec & scales_, 
		arma::vec & minima_,    arma::vec & maxima_, arma::vec & draws_, 
		trng::yarn2 & R_
) : locations(locations_), 
		tails(tails_),
		scales(scales_),
		minima(minima_),
		maxima(maxima_),
		R(R_),
		distributions(locations_.size()),
		sample_order(),
		draws(draws_) { }

arma::vec & Locations::state() const { return locations; }
double & Locations::state(arma::uword which) const { return locations[which]; }

// Available distributions:
void Locations::bind_constant_distribution	(
		unsigned int which
) {
	if (distributions[which] == NULL) {
		sample_order[0].push_back(which);
		distributions[which] = 
			std::unique_ptr<Random>(new RV_Constant(draws[which]));
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
		sample_order[1].push_back(which);
		distributions[which] = 
			std::unique_ptr<Random>(
				new RV_Uniform(draws[which], minima[which], maxima[which], R));
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
	// ONLY works if distributions[-1/+1] are drawn on the
	// prior "level".
	if (distributions[which] == NULL) {
		sample_order[2].push_back(which);
		distributions[which] = 
			std::unique_ptr<Random>(new RV_Uniform(
				draws[which], draws[which-1], draws[which+1], R));
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
		sample_order[2].push_back(which);
		distributions[which] = 
			std::unique_ptr<Random>(new RV_Missing_t_walk(
				draws[which-1], draws[which], draws[which+1], 
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
		sample_order[2].push_back(which);
		distributions[which] = 
			std::unique_ptr<Random>(new RV_Missing_t_walk(
				draws[which-1], draws[which], draws[which+1], 
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
		// Reverse lookup in sample_order to remove the value from the
		// correct key...
		for (unsigned int i = 0; i < sample_order.size(); ++i) {
			for (unsigned int j = 0; j < sample_order[i].size(); ++i) {
				if (sample_order[i][j] == which) 
					sample_order[i].erase(sample_order[i].begin()+j);
			}
		}
		distributions[which].reset(NULL);
	}	
	
}

void Locations::draw() {
	// Modify this loop to respect the sample_order vector...
	// Need to also add correct removal order.
	unsigned int which=0;
	for (unsigned int level=0; level < sample_order.size(); ++level) {
		for (unsigned int i = 0; i < sample_order[level].size(); ++i) { 
			which = sample_order[level][i];
			if (distributions[which] == NULL) {
				std::stringstream msg;
				msg << "The location " << which << " (" << (which+1) << ")"
							 " lacks has a distribution.  Not drawing.\n";
				throw(std::logic_error(msg.str()));
			} else {
				draws[which] = distributions[which]->draw();
			}
		}
	}
}

Locations::~Locations() { distributions.clear(); }
