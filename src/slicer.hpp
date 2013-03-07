#ifndef SLICER_H
#define SLICER_H

#include <armadillo>
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>
#include <trng/discrete_dist.hpp>


// This could be templated for return type?
template <class T_VAL, class T_PMF, class T_RET> 
class Slicer_Discrete {

public:
	Slicer_Discrete();
	Slicer_Discrete(
		const T_VAL * values,
		const T_PMF * pmf,
		trng::yarn2 * pR
	);
	T_RET draw();

private:
	trng::discrete_dist CH;
	const T_PMF * ppmf; // Pointer to the PMF.
	const T_VAL * val;
	
};

class Slicer_Continuous {

public:
	Slicer_Continuous();
	Slicer_Continuous( 
		double (*pdf)(double x),
		trng::yarn2 * pR
	);
	double draw();

private:
	trng::uniform01_dist<double> U;
	double sample;

};


#endif
