#ifndef RV_CONSTANT_H
#define RV_CONSTANT_H

#include "random.hpp"
#include <vector>
#include <map> 


class RV_Constant : public Random {

public:
	RV_Constant(double & X);

	void jump(double X);
	double draw();
	double lpdf(double X);
	double lpdf();

private:
	double const & x;

};

#endif