#include "rv_normal.hpp"

#include <math.h>

RV_Normal::RV_Normal(
		double       & X_,
		double const & mu_,
		double const & s_,
		trng::yarn2  & R_
) : X(X_), mu(mu_), s(s_), R(R_), NORMAL(mu_, s_) { }

void RV_Normal::jump(double X) { X = X; }

double RV_Normal::draw() {
	X = NORMAL(R);
	return X;
}

double RV_Normal::lpdf(double X) {
	NORMAL.mu(mu);
	NORMAL.sigma(s);
	return log(NORMAL.pdf(X));
}

double RV_Normal::lpdf() { lpdf(X); }

