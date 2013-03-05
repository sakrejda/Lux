#ifndef RECAPTURE_H
#define RECAPTURE_H

#include <vector>
#include <map>

#include <armadillo>

class Recapture_State {

public:
	Recapture_State();  // required
	Recapture_State(
		std::vector<int> times_of_surveys,  // ts
		std::vector<std::vector<int> > times_of_recaptures   //  build 'caught'
		std::vector<int> times_of_deaths,
		std::vector<bool> known_deaths
	);

	const	int& get_N() const;
	const int& get_K() const;
	const arma::Col<int>& get_surveys() const;
	const arma::Col<int>& get_recaptures(const arma::uvec& row_nums) const;
	const arma::Col<int>& get_recaptures(const arma::span& span = arma::span::all) const;
 	const arma::Col<int>& get_births(const arma::uvec & i_nums) const;
 	const arma::Col<int>& get_births(const arma::span& span = arma::span::all) const;
 	const arma::Col<int>& get_deaths(const arma::uvec & i_nums) const;
 	const arma::Col<int>& get_deaths(const arma::span& span = arma::span::all) const;
 	const arma::Col<int>& get_known_deaths(const arma::uvec & i_nums) const;
 	const arma::Col<int>& get_known_deaths(const arma::span& span = arma::span::all) const;
 	const arma::Col<int>& get_first_obs(const arma::uvec & i_nums) const;
 	const arma::Col<int>& get_first_obs(const arma::span& span = arma::span::all) const;
 	const arma::Col<int>& get_last_obs(const arma::uvec & i_nums) const;
 	const arma::Col<int>& get_last_obs(const arma::span& span = arma::span::all) const;


	void set_td( arma::Col<arma::uword> indexes, arma::Col<int> times_of_deaths);

private:
	arma::Col<int> ts;
	arma::Col<int> tb;
	arma::Col<int> td;
	arma::Col<int> fo;
	arma::Col<int> lo;
	arma::SpMat<int> caught;
	arma::SpMat<double> caught_double;
	arma::SpMat<int> uncaught;
	arma::SpMat<double> uncaught_double;
	arma::SpMat<double> available;
	int number_of_individuals;
	int number_of_occasions;
	std::vector<bool> known_death;

	void init();

};

class Recapture_Parameters {
	
public:
	Recapture_Parameters();
	Recapture_Parameters(
		Recapture_State const & state,
		unsigned int scale
	);

	void resize_PHI( unsigned int scale);
	
	const arma::Mat<double>& get_PHI();
	const arma::Mat<double>& get_PHI(const arma::uvec & row_nums);

	const arma::Mat<double>& get_P();
	const arma::Mat<double>& get_P(const arma::uvec & row_nums);

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
		Recapture_State const & state,
		Recapture_Parameters const & parameters
	);

	double get_likelihood(bool log=true);

private:
	Recapture_State const & state;
	Recapture_Parameters const & parameters;

	arma::Col<double> log_likelihood_phi;
	arma::Col<double> log_likelihood_p;
	
	void calc_ll_phi(const arma::uword & i);
	void calc_ll_p(const arma::uword & i;

	void init();

};

class Recapture_Priors {

public:
	Recapture_Priors();
	Recapture_Priors(
		Recapture_Parameters const & parameters_,
		double (*priors_)(Recapture_Parameters const & parameters)
	)
	~Recapture_Priors();

	double get_prior_density(bool log=true);

private:
	bool default_par;
	Recapture_Parameters const & parameters;
	double (*priors)(Recapture_Parameters const & parameters);
	


};


#endif
