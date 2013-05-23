
#include "locations.hpp"

#include <armadillo>
#include <trng/yarn2.hpp>

Locations::Locations(arma::vec & vec, trng::yarn2 & R_) : locations(vec), R(R_) { }

arma::vec Locations::state() const {
	arma::vec out = locations;
	return out;
}

Locations::~Locations() { delete vec; }
