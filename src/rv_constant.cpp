#include "rv_constant.hpp"
#include <math.h>
#include <limits>


RV_Constant::RV_Constant(double & X) : x(X) {}

double RV_Constant::draw() {return x;}

double RV_Constant::lpdf(double X) {
	if (X == x)
		return 0.0;
	else 
		return -1.0 * std::numeric_limits<double>::infinity();
}

double RV_Constant::lpdf() { return lpdf(x); }


