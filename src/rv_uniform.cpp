#include "rv_uniform.hpp"
#include <math.h>
#include <limits>
#include <stdexcept>


RV_Uniform::RV_Uniform(
		double & X, 
		double const & minimum, 
		double const & maximum, 
		trng::yarn2 & R_
) : x(X), min(minimum), max(maximum), R(R_) { 
	if ((x < min) || (x > max)) 
		out_of_domain();
}

double RV_Uniform::draw() { 
	x = min + (U(R) * (max - min)); 
	return x;
}

double RV_Uniform::lpdf(double X) {
	if ((X >= min) && (X <= max)) 
		return log(1.0/(max-min));
	else
		return -1.0 * std::numeric_limits<double>::infinity();
}

double RV_Uniform::lpdf() { return lpdf(x); }

void RV_Uniform::out_of_domain() {
		std::stringstream msg;
		msg << "Given value for 'x' is " << x << ", this is "
					 "not in [" << min << "," << max << "].";
		throw std::domain_error(msg.str());
}
