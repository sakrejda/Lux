#include "recapture.hpp"
#include <iostream>

//
//  Member functions for Recapture_State.
//
Recapture_State::Recapture_State() : number_of_individuals(0) {}

Recapture_State::Recapture_State(
    std::vector<int> times_of_surveys,
    std::vector<std::vector<int> > times_of_recaptures,
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
		caught.zeros();
    for ( int i=0; i < number_of_individuals; ++i ) {
      for ( int j=0; j < times_of_recaptures[i].size(); ++j ) {
				caught(i,times_of_recaptures[i][j]) = 1;          
			}
			if (known_deaths[i]) known_death(i) = 1; else known_death(i) = 0;
    }
    init();
}


const int& Recapture_State::get_N() const { return number_of_individuals; }
const int& Recapture_State::get_K() const { return number_of_occasions; }
const arma::Col<int>& Recapture_State::get_surveys()   const { return ts; }

const arma::Mat<int> & 
Recapture_State::get_recaptures() const { return caught; }

const arma::Col<int> & 
Recapture_State::get_births() const { return tb; }

const arma::Col<int> & 
Recapture_State::get_deaths() const { return td; }

const arma::Col<int> &
Recapture_State::get_known_deaths() const { return known_death; }

const arma::Col<int> & 
Recapture_State::get_first_obs() const { return fo; }

const arma::Col<int> & 
Recapture_State::get_last_obs() const { return lo; }


void Recapture_State::init() {
    for ( arma::uword i=0; i < caught.n_rows; ++i ) {
    		bool find_fo = true;
        for ( arma::uword j=0; j < caught.n_cols; ++j ) {
            if (find_fo && caught(i,j) == 1 ) { fo[i] = j; find_fo = false;}
            if (           caught(i,j) == 1 ) { lo[i] = j;                  }
						if ( caught(i,j) == 0 ) uncaught(i,j) = 1;
            if ( j > tb[i] && j < td[i] ) { available(i,j) = 1.0; }
        }
    }
    tb = fo;
}

void Recapture_State::set_td(arma::Col<int> times_of_deaths) {
	td = times_of_deaths;
}


//
//  Member functions for Recapture_Parameters.
//

Recapture_Parameters::Recapture_Parameters(
		Recapture_State const & state_,
		unsigned int scale = 5
) : PHI(state.get_N(), state.get_K() ),
    P(state.get_N(), state.get_K()),
		state(state_)
{
	resize_PHI(scale);
}

void Recapture_Parameters::resize_PHI( unsigned int scale) {
    unsigned int phi_rows = state.get_N();
		unsigned int phi_cols = state.get_K();
		PHI.zeros(phi_rows, scale * phi_cols);
		P.zeros(phi_rows, scale * phi_cols);
}

const arma::Mat<double>& Recapture_Parameters::get_PHI() const { return PHI; }
const arma::Mat<double>& Recapture_Parameters::get_P() const { return P; }

void Recapture_Parameters::set_PHI( arma::Mat<double> PHI_ ) { PHI = PHI_; }
void Recapture_Parameters::set_P(   arma::Mat<double> P_   ) { P   = P_; }


//
//  Member functions for Recapture_Likelihood;
//
//

Recapture_Likelihood::Recapture_Likelihood(
		Recapture_State const & state_,
		Recapture_Parameters const & parameters_
) : state(state_), parameters(parameters_)
{
	int rows = (parameters.get_PHI()).n_rows;
	log_likelihood_phi.zeros(rows);
	log_likelihood_p.zeros(rows);
}

double Recapture_Likelihood::get_likelihood(bool log) {
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

void Recapture_Likelihood::calc_ll_phi(const arma::uword & i) {
	log_likelihood_phi[i] = 0.0;
	int tb = arma::as_scalar(state.get_births()(i));
	int td = arma::as_scalar(state.get_deaths()(i));
	const arma::Mat<double>& PHI = parameters.get_PHI();
	for(unsigned int t=tb; t < td-1; ++t) {
		log_likelihood_phi[i] += log(PHI(i,t));
	}
	if (!state.get_known_deaths()(i)) 
		log_likelihood_phi[i] += log(1-PHI(i,td-1));
}

void Recapture_Likelihood::calc_ll_p(const arma::uword & i) {
	log_likelihood_p[i] = 0.0;
	int tb = arma::as_scalar(state.get_births()(i));
	int td = arma::as_scalar(state.get_deaths()(i));
	const arma::Mat<double>& P = parameters.get_P();
	const arma::Mat<int>& caught = state.get_recaptures();
	for (unsigned int t=tb+1; t < td; ++t) {	
		if (caught(i,t) == 1) {
			log_likelihood_p[i] += log(P(i,t));
		} else {
			log_likelihood_p[i] += log(1-P(i,t));
		}
	}
}


//
//  Member functions for Recapture_td_Posterior
//
//

Recapture_td_Posterior::Recapture_td_Posterior(
	Recapture_State const & state_,
	Recapture_Parameters const & parameters_,
	trng::yarn2 & R_
) : parameters(parameters_), state(state_), R(R_) {
	N = parameters.get_PHI().n_rows;
	K = parameters.get_PHI().n_cols;
	td.set_size(N);
	slicers.set_size(N);
	S.resize(N, K);
	D.resize(N, K);
	td_PMF.set_size(N);
	choices.resize(K);
	for ( unsigned int i=0; i < K; ++i ) {
		choices(i) = i;
	}
	Slicer_Discrete<arma::Row<int>, arma::Row<double>, int> * s;
	for ( unsigned int i=0; i < N; ++i ) {
		td_PMF(i).resize(K);
		s = new Slicer_Discrete<arma::Row<int>, arma::Row<double>, int>(&choices, &td_PMF(i), &R);
		slicers(i) = s;
	}

}

arma::Col<int> Recapture_td_Posterior::draw() {
	calc_log_mass_function();
	for ( unsigned int i=0; i < N; ++i ) {
		td(i) = (*slicers(i)).draw();
	}
	return td;	
}

arma::field<arma::Row<double> > Recapture_td_Posterior::calc_log_mass_function() {
	const arma::Col<int> & lo = state.get_last_obs();
	const arma::Mat<double> & PHI = parameters.get_PHI();
	const arma::Mat<double> & P = parameters.get_P();
	S.zeros();
	D.zeros();

	for ( unsigned int i=0; i < N; ++i ) {
		td_PMF(i).zeros();
		S(i,lo[i]+1) = 0.0;
		D(i,lo[i]+1) = log( 1-PHI(i,lo[i]) );
		td_PMF(i)(lo[i]+1) = exp(  S(i,lo[i]+1) + D(i,lo[i]+1)  ); // exp(
		for ( unsigned int t=lo[i]+2; t < K; ++t ) {
			S(i,t) = log( PHI(i,t-2) ) + S(i,t-1) + log( 1 - P[i,t-1]);
			D(i,t) = log( 1-PHI(i,t-1) );
			td_PMF(i)(t) = exp(  S(i,t) + D(i,t)  ); // exp(
		}
	}
	return td_PMF;
}

