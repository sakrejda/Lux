#ifndef RANDOM_H
#define RANDOM_H

#include <vector>
#include <map> 

#include <armadillo>
#include <Eigen/Dense>
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>
#include <trng/exponential_dist.hpp>
//#include <trng/student_t_dist.hpp>

//#include "slicer-discrete.hpp"
//#include "slicer-continuous.hpp"

class Random { 

public:
	virtual void jump(double X) = 0;
	virtual double draw() = 0;
	virtual double lpdf(double X) = 0;
	virtual double lpdf() = 0;	

private:

};

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

class RV_Uniform : public Random {

public:
	typedef double (*PDF)(double);

	RV_Uniform(
		double 			 & X,
		double const & minimum,
		double const & maximum,
		trng::yarn2 & R_);

	void jump(double X);
	double draw();
	double lpdf(double X);
	double lpdf();

private:
		trng::yarn2  & R;
		double 			 & x;
		double const & min;
		double const & max;
  	trng::uniform01_dist<> U;

};

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
	std::vector<double> bounds;
	trng::yarn2 & R;  
	trng::uniform01_dist<double> U;
	trng::exponential_dist<double> EXPO;

	double const & x1; 
	double       & x2;
	double const & os;
	double const & p1;
	double const & s1;

	std::vector<double> step_out();
	double ly; // cached log pdf at current value...
	double x_new; // temporary new sampled value...
	double ly_new; // temporary pdf at sampled value...

};

class RV_Missing_t_walk : public Random {

public:
	typedef double (*PDF)(double);

	RV_Missing_t_walk(
		double const & x1_,
		double 			 & X,
		double const & x3_,
		double const & p1_,   // degrees of freedom 1
		double const & p2_,   // degrees of freedom 2
		double const & s1_,   // scale 1
		double const & s2_,   // scale 2
		trng::yarn2 & R_
	);

	std::map<std::string, double> state() const;
	void jump(double X);
	double draw();
	double lpdf(double X);
	double lpdf();

private:
	void find_peaks();
	double peak1, peak2;
	arma::Mat<double> companion;
	arma::cx_vec cx_eigval;
	arma::vec    eigvalues;
	arma::cx_mat cx_eigvec;

	Eigen::MatrixXd ecompanion;

	void print_slice(std::string s);

	void find_slice();
	std::vector<double> bounds_pk1;
	std::vector<double> bounds_pk2;
	std::vector<double> step_out(double peak);
	double choose();
	void trim();
	double ly; // cached log pdf at current value...
	double x_new; // temporary new sampled value...
	double ly_new; // temporary pdf at sampled value...

	trng::yarn2 & R;  
	trng::uniform01_dist<double> U;
	trng::exponential_dist<double> EXPO;

	double const & x1; 
	double       & x2;
	double const & x3;
	double const & p1;
	double const & p2;
	double const & s1;
	double const & s2;

};


#endif
