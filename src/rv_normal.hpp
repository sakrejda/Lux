#ifndef RV_NORMAL_H
#define RV_NORMAL_H

#include "random.hpp"

#include <map> 

#include <trng/yarn2.hpp>
#include <trng/normal_dist.hpp>

class RV_Normal : public Random {

public:
	RV_Normal(
		double       & X_,
		double const & mu_,
		double const & s_,
		trng::yarn2  & R_
	);
	std::map<std::string, double> state() const;
	double draw();
	double lpdf(double X);
	double lpdf();

private:
	trng::yarn2 & R;  
	trng::normal_dist<> NORMAL;

	double       & X;
	double const & mu;
	double const & s;

};


#endif
