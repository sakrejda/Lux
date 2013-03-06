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
		const arma::Col<int> * values,
		const arma::Col<double> * pmf,
		trng::yarn2 * pR
	);
	unsigned int draw();

private:
	trng::discrete_dist CH;
	const arma::Col<double> * ppmf; // Pointer to the PMF.
	const arma::Col<int> * val;
	
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
