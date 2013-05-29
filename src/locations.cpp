
#include "locations.hpp"

#include <armadillo>
#include <trng/yarn2.hpp>

Locations::Locations(arma::vec & vec, trng::yarn2 & R_) : locations(vec), R(R_) { }

arma::vec & Locations::state() const { return locations; }
double & Locations::state(arma::uword which) const { return locations[which]; }

Locations::~Locations() { }
