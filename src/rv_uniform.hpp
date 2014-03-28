#ifndef RV_UNIFORM_H
#define RV_UNIFORM_H

#include "random.hpp"

#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>

class RV_Uniform : public Random {

public:
	RV_Uniform(
		double 			 & X,
		double const & minimum,
		double const & maximum,
		trng::yarn2 & R_);

	void draw();
	double lpdf(double X);
	double lpdf();

private:
		trng::yarn2  & R;
		double 			 & x;
		double const & min;
		double const & max;
  	trng::uniform01_dist<> U;
		void out_of_domain();

};

#endif
