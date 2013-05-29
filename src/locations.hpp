#ifndef LOCATION_H
#define LOCATION_H

#include <vector>
#include <map>

#include <armadillo>
#include <trng/yarn2.hpp>

class Locations {

public:
	Locations(arma::vec & vec, trng::yarn2 & R_);
	arma::vec & state() const;
	double & state(arma::uword which) const;
	~Locations();  // delete vec, which is just a wrapper for 
								// some R memory, but should not leak.

private:
	trng::yarn2 & R;
	arma::vec & locations;

};

#endif
