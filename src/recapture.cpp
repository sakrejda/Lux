#include "survival_2.hpp"
#include <iostream>

//
//  Member functions for Recapture_State.
//
Recapture_State::Recapture_State() : number_of_individuals(0) {}

Recapture_State::Recapture_State(
    std::vector<int> times_of_surveys,
    std::vector<std::vector<int> > times_of_recaptures
		std::vector<int> times_of_deaths,
		std::vector<bool> known_deaths
) : ts(times_of_surveys), 
		number_of_individuals(times_of_recaptures.size()),
		number_of_occasions(*max_element(times_of_surveys.begin(),times_of_surveys.end())+1),
    fo(times_of_recaptures.size()), lo(times_of_recaptures.size()),
    caught(times_of_recaptures.size(),
					 *max_element(times_of_surveys.begin(),times_of_surveys.end())+1),
    uncaught(times_of_recaptures.size(),
					 *max_element(times_of_surveys.begin(),times_of_surveys.end())+1),
		known_death(times_of_recaptures.size()),
    tb(times_of_recaptures.size()),
		td(times_of_deaths),
    available(times_of_recaptures.size(),times_of_surveys.size())
{
    for ( int i=0; i < number_of_individuals; ++i ) {
        for ( int j=0; j < times_of_recaptures[i].size(); ++j ) {
            caught(i,times_of_recaptures[i][j]) = 1;          }
    }
    init();
}


const int& Recapture_State::get_N() const { return number_of_individuals; }
const int& Recapture_State::get_K() const { return number_of_occasions; }
const arma::Col<int>& Recapture_State::get_surveys()   const { return ts; }


const arma::Col<int> & 
Recapture_State::get_recaptures(const arma::uvec& row_nums) const { 
	return caught.rows(row_nums);
}

const arma::Col<int> & 
Recapture_State::get_recaptures(const span& span = arma::span::all) const { 
	return caught.rows(span);
}

const arma::Col<int>& get_births(const arma::uvec & i_nums) const {
	return tb.elem(i_nums);
}
const arma::Col<int>& get_births(const span & span = arma::span::all) const {
	return tb(span);
}

const arma::Col<int>& get_deaths(const arma::uvec & i_nums) const {
	return td.elem(i_nums);
}
const arma::Col<int>& get_deaths(const span & span = arma::span::all) const {
	return td(span);
}

const arma::Col<int>& get_first_obs(const arma::uvec & i_nums) const {
	return fo.elem(i_nums);
}
const arma::Col<int>& get_first_obs(const span & span = arma::span::all) const {
	return fo(span);
}

const arma::Col<int>& get_last_obs(const arma::uvec & i_nums) const {
	return lo.elem(i_nums);
}
const arma::Col<int>& get_last_obs(const span & span = arma::span::all) const {
	return lo(span);
}


void Recapture_State::init() {
    for ( arma::uword i=0; i < caught.n_rows; ++i ) {
    		bool find_fo = true;
        for ( arma::uword j=0; j < caught.n_cols; ++j ) {
            if (find_fo && caught(i,j) == 1 ) { fo[i] = j; find_fo = false;}
            if (           caught(i,j) == 1 ) { lo[i] = j;                  }
						caught_double(i,j) == (double)caught(i,j);
						if ( caught(i,j) == 0 ) uncaught(i,j) = 1;
						uncaught_double(i,j) == (double)caught(i,j);
            if ( j > tb[i] && j < td[i] ) {
                available(i,j) = 1.0;
            }
        }
				known_death[i] = known_deaths[i];
    }
    tb = fo;
}

void Recapture_State::set_td(
        arma::Col<arma::uword> indexes,
        arma::Col<int> times_of_deaths
) {
    for ( int k=0; k < indexes.n_elem; ++k ) {
        td[indexes[k]] = times_of_deaths[k];
    }
}


//
//  Member functions for Recapture_Parameters.
//

Recapture_Parameters::Recapture_Parameters() : PHI(1,1), P(1,1) {}

Recapture_Parameters::Recapture_Parameters(
		Recapture_State const & state,
		unsigned int scale = 5
) : PHI(state.get_N(), state.get_K() ),
    P(state.get_N(), state.get_K())
{
	resize_PHI(scale);
}

void Recapture_Parameters::resize_PHI( unsigned int scale) {
    unsigned int phi_rows = state.get_N();
		unsigned int phi_cols = state.get_K();
		PHI.zeros(phi_rows, scale * phi_cols);
		P.zeros(phi_rows, scale * phi_cols);
}

const arma::Mat<double>& Recapture_Parameters::get_PHI() { return PHI; }
const arma::Mat<double>& Recapture_Parameters::get_P() { return P; }

void Recapture_Parameters::set_PHI( arma::Mat<double> PHI_ ) { PHI = PHI_; }
void Recapture_Parameters::set_P(   arma::Mat<double> P_   ) { P   = P_; }


//
//  Member functions for Recapture_Likelihood;
//
//
Recapture_Likelihood::Recapture_Likelihood() {}

Recapture_Likelihood::Recapture_Likelihood(
		Recapture_State const & state,
		Recapture_Parameters const & parameters
) :
{
	int rows = parameters.get_PHI().n_rows
	log_likelihood_phi.zeros(rows);
	log_likelihood_p.zeros(rows);
}

double get_likelihood(bool log=true) {
	for (unsigned int i=0; i < state.get_N(); ++i) {
		calc_ll_phi(i); calc_ll_p(i);
	}
	double log_likelihood = 
		arma::accu(log_likelihood_phi) + arma::accu(log_likelihood_p);
	if (log) 
		return log_likelihood;
	else
		return exp(log_likelihood);
}

void calc_ll_phi(const arma::uword & i) {
	log_likelihood_phi[i] = 0.0;
	int tb = arma::as.scalar(state.get_births(i));
	int td = arma::as.scalar(state.get_deaths(i));
	for(unsigned int t=tb; t < td-1; ++t) {
		log_likelihood_phi[i] += log(parameters.get_PHI()(i,t));
	}
	if (!state.get_known_deaths(i)) 
		log_likelihood_phi[i] += log(1-parameters.get_PHI()(i,td-1));
}

void calc_ll_p(const arma::uword & i) {
	log_likelihood_p[i] = 0.0;
	int tb = arma::as.scalar(state.get_births(i));
	int td = arma::as.scalar(state.get_deaths(i));
	for (unsigned int t=tb+1; t < td; ++t) {	
		if (caught(i,t) == 1) {
			ll_p_components[i] += log(parameters.get_P()(i,t));
		} else {
			ll_p_components[i] += log(1-parameters.get_P()(i,t));
		}
	}
}











// Posterior functions. 
Recapture_Posterior::Recapture_Posterior() : 
	Recapture_Likelihood() {}

Recapture_Posterior::Recapture_Posterior(
		std::vector<int> times_of_surveys,
		std::vector<std::vector<int> > times_of_recaptures,
		std::vector<int> times_of_deaths,
		std::vector<bool> known_deaths
) : Recapture_Likelihood(times_of_surveys, times_of_recaptures,
			times_of_deaths, known_deaths) {
	init();
}

double Recapture_Posterior::get_lp() { return 0; }   /// Flat priors....


double Recapture_Posterior::get_log_posterior() {
	return get_ll() + get_lp();
}


void Recapture_Posterior::init() {}

