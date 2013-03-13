#ifndef RECAPTURE_H
#define RECAPTURE_H

#include <vector>
#include <map>

#include <armadillo>
#include <trng/yarn2.hpp>

#include "slicer.hpp"

class Recapture_State {

public:
	Recapture_State();  // required
	Recapture_State(
		std::vector<int> times_of_surveys,  // ts
		std::vector<std::vector<int> > times_of_recaptures,   //  build 'caught'
		std::vector<int> times_of_deaths,
		std::vector<bool> known_deaths
	);

	const	int& get_N() const;
	const int& get_K() const;
	const arma::Col<int> & get_surveys() const;
	const arma::Mat<int> & get_recaptures() const; 
	const arma::Col<int> & get_births() const;
	const arma::Col<int> & get_deaths() const;
	const arma::Col<int> & get_first_obs() const;
	const arma::Col<int> & get_last_obs() const;
	const arma::Col<int> & get_known_deaths() const;

	void set_td(arma::Col<int> times_of_deaths);

private:
	arma::Col<int> ts;
	arma::Col<int> tb;
	arma::Col<int> td;
	arma::Col<int> fo;
	arma::Col<int> lo;
	arma::Mat<int> caught;
	arma::Mat<int> uncaught;
	arma::Mat<double> available;
	int number_of_individuals;
	int number_of_occasions;
	arma::Col<int> known_death;

	void init();

};

class Recapture_Parameters {
	
public:
	Recapture_Parameters();
	Recapture_Parameters(
		Recapture_State const & state_,
		unsigned int scale
	);

	void resize_PHI( unsigned int scale);
	
	const arma::Mat<double>& get_PHI() const;
	const arma::Mat<double>& get_PHI(const arma::uvec & row_nums) const;

	const arma::Mat<double>& get_P() const;
	const arma::Mat<double>& get_P(const arma::uvec & row_nums) const;

	void set_PHI( arma::Mat<double> PHI_);
	void set_P( 	arma::Mat<double> P_ );

private:
	Recapture_State const & state;
	arma::Mat<double> PHI;
	arma::Mat<double> P;

};


class Recapture_Likelihood {

public:
	Recapture_Likelihood();
	Recapture_Likelihood(
		Recapture_State const & state_,
		Recapture_Parameters const & parameters_
	);

	double get_likelihood(bool log=true);

private:
	Recapture_State const & state;
	Recapture_Parameters const & parameters;

	arma::Col<double> log_likelihood_phi;
	arma::Col<double> log_likelihood_p;
	
	void calc_ll_phi(const arma::uword & i);
	void calc_ll_p(const arma::uword & i);

	void init();

};

class Recapture_td_Posterior {	

public:
	Recapture_td_Posterior(
		Recapture_State const & state_,
		Recapture_Parameters const & parameters_,
		trng::yarn2 & R_
	);

	arma::Col<int> draw();
	arma::field<arma::Row<double> > calc_log_mass_function();

private:
	// Acutally, maybe a better strategy is to keep just a reference to
	// the relevant state members?
	trng::yarn2 & R;
	Recapture_Parameters const & parameters;
	Recapture_State const & state;
	arma::Mat<double> S;
	arma::Mat<double> D;
	arma::field< arma::Row<double> > td_PMF;
	arma::Row<int> choices;
	arma::Col<int> td;
	arma::field<Slicer_Discrete<arma::Row<int>, arma::Row<double>, int> * > slicers;
	unsigned int N, K;

};
#endif
