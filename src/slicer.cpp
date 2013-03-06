#include "slicer.hpp"

Slicer_Discrete::Slicer_Discrete() : 
	ppmf(0), val(0), CH(0), pR(0) {}

Slicer_Discrete::Slicer_Discrete( 
	const arma::Col<int> * values,
	const arma::Col<double> * pmf,
	trng::yarn2 * pR;
	) : 
	ppmf(pmf), val(values), CH((*pmf).n_elem) {} 

unsigned int Slicer_Discrete::draw() {
	// for_each?
	for (unsigned int i=0; i < (*ppmf).n_elem; ++i) CH.param(i,(*ppmf)[i]);
	return (*val)[CH(*pR)];	
}

Slicer_Continuous::Slicer_Continuous() {}

Slicer_Continuous::Slicer_Continuous( double (*pdf)(double x) ) {}

double Slicer_Continuous::draw() {}

