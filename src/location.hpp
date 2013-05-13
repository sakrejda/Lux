#ifndef LOCATION_H
#define LOCATION_H

#include <vector>
#include <map>

#include <armadillo>
#include <trng/yarn2.hpp>

#include "slicer-discrete.hpp"
#include "slicer-continuous.hpp"

class t_walk_Posterior {

public:
	typedef double (*PDF)(double);

	t_walk_Posterior(
		double const & x1_,
		double 			 & x2_,
		double const & x3_,
		double const & p1_,   // degrees of freedom 1
		double const & p2_,   // degrees of freedom 2
		double const & s1_,   // scale 1
		double const & s2_,   // scale 2
		trng::yarn2 & R_
	);

	std::map<std::string, double> state() const;
	double draw();
	void jump(double x2_);
	void slice(double lpdf_);
	double lpdf(double x2_);
	double lpdf();

private:
	trng::yarn2 & R;
	double const & x1; 
	double const & x2;
	double const & x3;
	double const & p1;
	double const & p2;
	double const & s1;
	double const & s2;

};



class Location_State {

public:
	Location_State();
	Location_State(
		arma::Mat<double> location;
	);



};


























