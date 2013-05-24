#ifndef SLICER_CONTINUOUS_H
#define SLICER_CONTINUOUS_H

#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>
#include <trng/exponential_dist.hpp>
#include <cmath>
#include <algorithm>

template <class T_PDF, class T_DOMAIN> 
class Slicer_Continuous {

public:
	Slicer_Continuous();
	Slicer_Continuous(
		const T_PDF lpdf,
		T_DOMAIN x_,
		T_DOMAIN domain_min,
		T_DOMAIN domain_max,
		trng::yarn2 * pR
	);
	T_DOMAIN draw();
	void jump(T_DOMAIN x_);

private:
	const T_PDF plpdf;
	trng::yarn2 * R;
	trng::uniform01_dist<double> U;
	trng::exponential_dist<double> EXPO;
	T_DOMAIN min, max;

	T_DOMAIN w, l_bound, r_bound;  // width of stepping-out steps, and bounds
	std::vector<T_DOMAIN> w_vector;
	unsigned int m;  // max number of step-out steps to take
	T_DOMAIN x, x_new;    // location on density
	double ly, ly_new;    // runif(0,f(x))
	void step_out();

};


template <class T_PDF, class T_DOMAIN> 
Slicer_Continuous<T_PDF, T_DOMAIN>::Slicer_Continuous() : 
	plpdf(0), min(0), max(0), R(0), EXPO(1.0), w(1.0), m(10) {}

template <class T_PDF, class T_DOMAIN> 
Slicer_Continuous<T_PDF, T_DOMAIN>::Slicer_Continuous(
	const T_PDF lpdf,
	T_DOMAIN x_,
	T_DOMAIN domain_min, 
	T_DOMAIN domain_max,
	trng::yarn2 * pR
) : plpdf(lpdf), min(domain_min), max(domain_max), R(pR), EXPO(1.0),
		w(1.0), m(10), x(x_), w_vector(100) {

	for (unsigned int i=0; i < 100; ++i) {
		draw();
		w_vector[i] = r_bound - l_bound;
	}
	std::nth_element(w_vector.begin(), w_vector.begin() + w_vector.size()/2, w_vector.end());
	double width_estimate = w_vector[w_vector.size()/2];
	w = width_estimate/5;
	jump(x_);

}

template <class T_PDF, class T_DOMAIN> 
T_DOMAIN
Slicer_Continuous<T_PDF, T_DOMAIN>::draw() {
//	y = U(*R) * plpdf(x);
  ly = plpdf(x) - EXPO(*R);
	step_out();
	while(true) {
		x_new = l_bound + (U(*R) * std::abs(r_bound - l_bound));
		ly_new = plpdf(x_new);
		if ((ly_new > ly) && (x_new > min) && (x_new < max)) 
			break; 
		else {
			if (x_new <= x) {
				l_bound = x_new;
			} else {
				r_bound = x_new;
			}
		}
	}
	x = x_new;
	return x;
}

template <class T_PDF, class T_DOMAIN> 
void Slicer_Continuous<T_PDF, T_DOMAIN>::step_out() {
	l_bound = x - w * U(*R);
	r_bound = l_bound + w;
	int j = std::floor(m * U(*R)); 
	int k = (m-1) - j;
	while ((j>0) && ly < plpdf(l_bound) && (min < l_bound)) {
		l_bound = l_bound - w;
		j = j - 1;
	}
	while ((k>0) && ly < plpdf(r_bound) && (max > r_bound)) {
		r_bound = r_bound + w;
		k = k - 1;
	}
}

template <class T_PDF, class T_DOMAIN> 
void Slicer_Continuous<T_PDF, T_DOMAIN>::jump(T_DOMAIN x_) {
	if ( (x_ > min) && (x_ < max)) x = x_; 
}

#endif
