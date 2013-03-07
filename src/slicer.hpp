#ifndef SLICER_H
#define SLICER_H

#include <armadillo>
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>
#include <trng/discrete_dist.hpp>


// This could be templated for return type?
class Slicer_Discrete {

public:
	Slicer_Discrete();
	Slicer_Discrete(
		const arma::Row<int> * values,
		const arma::Row<double> * pmf,
		trng::yarn2 * pR
	);
	unsigned int draw();

private:
	trng::discrete_dist CH;
	const arma::Row<double> * ppmf; // Pointer to the PMF.
	const arma::Row<int> * val;
	
};

class Slicer_Continuous {

public:
	Slicer_Continuous();
	Slicer_Continuous( 
		double (*pdf)(double x)
		trng::yarn2 * pR
	);
	double draw();

private:
	trng::uniform01_dist<double> U;
	double sample;

};


#endif
