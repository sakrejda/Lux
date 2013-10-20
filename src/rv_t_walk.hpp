#ifndef RV_T_WALK_H
#define RV_T_WALK_H

#include "random.hpp"

#include <vector>
#include <map> 

#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>
#include <trng/exponential_dist.hpp>


class RV_t_walk : public Random {
	friend double LPDF(double x);

public:
	RV_t_walk(
		double const & x1_,
		double			 & X,
		double const & os_,
		double const & p1_,
		double const & s1_,
		trng::yarn2  & R_
	);
	std::map<std::string, double> state() const;
	void jump(double X);
	double draw();
	double lpdf(double X);
	double lpdf();

private:
	trng::yarn2 & R;  
	trng::uniform01_dist<double> U;
	trng::exponential_dist<double> EXPO;

	double const & x1; 
	double       & x2;
	double const & os;
	double const & p1;
	double const & s1;


};

#endif
