#ifndef SLICER_H
#define SLICER_H

#include <armadillo>
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>
#include <trng/discrete_dist.hpp>
#include <trng/exponential_dist.hpp>
#include <iostream>

template <class T_VAL, class T_PMF, class T_RET> 
class Slicer_Discrete {

public:
	Slicer_Discrete();
	Slicer_Discrete(
		unsigned int start,
		const T_VAL * values,
		const T_PMF * lpmf,
		trng::yarn2 * pR
	);
	T_RET draw();

private:
	trng::discrete_dist CH;
	trng::exponential_dist<double> EXPO;
	trng::uniform01_dist<> U;
	const T_PMF * plpmf; // Pointer to the PMF.
	const T_VAL * val;
	trng::yarn2 * R;

	std::vector<unsigned int> indexes;
	unsigned int x, x_new;
	double ly;
	
};



template <class T_VAL, class T_PMF, class T_RET> 
Slicer_Discrete<T_VAL, T_PMF, T_RET>::Slicer_Discrete() : 
	x(0), x_new(0), indexes(0),
	plpmf(0), val(0), CH(0), EXPO(0), R(0), EXPO(1.0) {}

template <class T_VAL, class T_PMF, class T_RET> 
Slicer_Discrete<T_VAL, T_PMF, T_RET>::Slicer_Discrete( 
	unsigned int start,
	const T_VAL * values,
	const T_PMF * lpmf,
	trng::yarn2 * pR
	) : 
	x(start), x_new(0), indexes(0),
	plpmf(lpmf), val(values), CH((*lpmf).n_elem), R(pR), EXPO(1.0) 
{
	indexes.reserve((*lpmf).n_elem);

} 

template <class T_VAL, class T_PMF, class T_RET> 
T_RET Slicer_Discrete<T_VAL, T_PMF, T_RET>::draw() {
	ly = (*plpmf)(x) - EXPO(*R);
	indexes.clear();
	for (unsigned int i=0; i < (*plpmf).n_elem; ++i) {
		if ((*plpmf)[i] > ly) indexes.push_back(i);
	}
	int ch = static_cast<int>(std::floor(U(*R) * indexes.size()));
	x = indexes[ch];
	T_RET out = (*val)(x);
	return out;	
}

#endif
