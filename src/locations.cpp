
#include "locations.hpp"

#include <armadillo>
#include <trng/yarn2.hpp>

Location::Location(arma::vec & vec, trng::yarn2 & R_) : locations(vec), R(R_) { }

arma::vec Location::state() const {
	arma::vec out = locations;
	return out;
}

