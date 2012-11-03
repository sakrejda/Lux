#include "survival.hpp"
#include <iostream>
// Member functions for Data_FLAT only:
//
Recapture_Data_FLAT::Recapture_Data_FLAT() : number_of_individuals(0) {}

Recapture_Data_FLAT::Recapture_Data_FLAT(
    std::vector<int> times_of_surveys,
    std::vector<std::vector<int> > times_of_recaptures
) : ts(times_of_surveys), 
		number_of_individuals(times_of_recaptures.size()),
		number_of_occasions(*max_element(times_of_surveys.begin(),times_of_surveys.end())+1),
    fo(times_of_recaptures.size()), lo(times_of_recaptures.size()),
    caught(times_of_recaptures.size(),
					 *max_element(times_of_surveys.begin(),times_of_surveys.end())+1),
    uncaught(times_of_recaptures.size(),
					 *max_element(times_of_surveys.begin(),times_of_surveys.end())+1),
		known_death(times_of_recaptures.size()),
    tb(times_of_recaptures.size()) 
{
    for ( int i=0; i < number_of_individuals; ++i ) {
        for ( int j=0; j < times_of_recaptures[i].size(); ++j ) {
            caught(i,times_of_recaptures[i][j]) = 1;          }
    }
    init();
}


int Recapture_Data_FLAT::get_N() const { return number_of_individuals; }
int Recapture_Data_FLAT::get_K() const { return number_of_occasions; }

arma::Row<int> Recapture_Data_FLAT::get_recaptures(int i) const { 
	arma::Row<int> recaptures(caught.n_cols);
	for ( arma::uword j=0; j < caught.n_cols; ++j ) {
		recaptures(j) = caught(i,j);
	}
	return recaptures;
}

arma::Row<int> Recapture_Data_FLAT::get_surveys()   const { return ts; }

arma::Col<int> Recapture_Data_FLAT::get_births()    const { return tb; }
arma::Col<int> Recapture_Data_FLAT::get_first_obs() const { return fo; }
arma::Col<int> Recapture_Data_FLAT::get_last_obs()  const { return lo; }

std::vector<bool> Recapture_Data_FLAT::get_sampled() { 
	std::vector<bool> sampled_vec(sampled.size());
	for ( unsigned int i=0; i < sampled.size(); ++i ) sampled_vec[i] = sampled[i];
	return sampled_vec; 
}

void Recapture_Data_FLAT::init() {
    for ( arma::uword i=0; i < caught.n_rows; ++i ) {
    		bool find_fo = true;
        for ( arma::uword j=0; j < caught.n_cols; ++j ) {
            if (find_fo && caught(i,j) == 1 ) { fo[i] = j; find_fo = false;}
            if (           caught(i,j) == 1 ) { lo[i] = j;                  }
						caught_double(i,j) == (double)caught(i,j);
						if ( caught(i,j) == 0 ) uncaught(i,j) = 1;
						uncaught_double(i,j) == (double)caught(i,j);
        }
				known_death[i] = false;
    }
    tb = fo;
		for (unsigned int t=0; t < get_K(); ++t) sampled[t] = false;
		for (unsigned int i=0; i < ts.size(); ++i) sampled[ts[i]] = true;
}

// Member functions for State_FLAT only:
//
Recapture_State_FLAT::Recapture_State_FLAT() :
    Recapture_Data_FLAT() { init(); }

Recapture_State_FLAT::Recapture_State_FLAT(
    std::vector<int> times_of_surveys,
    std::vector<std::vector<int> > times_of_recaptures,
		std::vector<int> times_of_deaths,
		std::vector<bool> known_deaths
) : Recapture_Data_FLAT(times_of_surveys, times_of_recaptures),
        td(times_of_deaths),
        available(times_of_recaptures.size(),times_of_surveys.size()) {
    init();
}

void Recapture_State_FLAT::set_td(
        arma::Col<arma::uword> indexes,
        arma::Col<int> times_of_deaths
) {
    for ( int k=0; k < indexes.n_elem; ++k ) {
        td[indexes[k]] = times_of_deaths[k];
    }
}

void Recapture_State_FLAT::known_deaths( arma::Col<arma::uword> indexes ) {
	for ( int k=0; k < indexes.n_elem; ++k ) {
		known_death[indexes[k]] = true;
	}
}

arma::Col<int> Recapture_State_FLAT::get_deaths()    const { return td; }


void Recapture_State_FLAT::init() {
    for (int i=0; i < caught.n_rows; ++i) {
        for (int j=0; j < caught.n_cols; ++j) {
            if ( j > tb[i] && j < td[i] ) {
                available(i,j) = 1.0;
            }
        }
    }
}


// Member functions for Likelihood_FLAT only:
//
Recapture_Likelihood_FLAT::Recapture_Likelihood_FLAT(
) : Recapture_State_FLAT(), PHI(0,0), P(0,0),
    log_likelihood(0), SES(0), fresh_ll(false)
{
    init();
}

//number_of_occasions(*max_element(times_of_surveys.begin(),times_of_surveys.end())+1),
Recapture_Likelihood_FLAT::Recapture_Likelihood_FLAT(
	std::vector<int> times_of_surveys,
	std::vector<std::vector<int> > times_of_recaptures,
	std::vector<int> times_of_deaths,
	std::vector<bool> known_deaths
) : Recapture_State_FLAT(times_of_surveys, times_of_recaptures,
			times_of_deaths, known_deaths),
    PHI(times_of_recaptures.size(),
				*max_element(times_of_surveys.begin(),times_of_surveys.end())),
    P(times_of_recaptures.size(),
				*max_element(times_of_surveys.begin(),times_of_surveys.end())+1),
    log_likelihood(0), 
		ll_phi_components(times_of_recaptures.size()), 
		ll_p_components(times_of_recaptures.size()), 
		SES(3), 
		fresh_ll(false),
		fresh_ll_p_components(times_of_recaptures.size()),
		fresh_ll_phi_components(times_of_recaptures.size())
{
    init();
}

void Recapture_Likelihood_FLAT::resize_PHI( unsigned int scale) {
    unsigned int phi_rows = number_of_individuals;
		unsigned int phi_cols = 0;
    for ( int i=0; i < tb.n_elem; ++i ) {
    	if ( td[i] > phi_cols) phi_cols = td[i];
    }
		PHI.set_size(								phi_rows, scale * phi_cols);
    ll_phi_components.set_size(phi_rows);
}

arma::Mat<double> Recapture_Likelihood_FLAT::get_PHI() { return PHI; }
arma::Mat<double> Recapture_Likelihood_FLAT::get_P() { return P; }

void Recapture_Likelihood_FLAT::set_PHI( arma::Mat<double> PHI_ ) { 
	PHI = PHI_; 
	fresh_ll = false;
}
void Recapture_Likelihood_FLAT::set_P(   arma::Mat<double> P_   ) { 
	P   = P_; 
	fresh_ll = false;
}


double Recapture_Likelihood_FLAT::get_ll() { 
	if (fresh_ll) {
		return log_likelihood; 
	} else {
		update_ll();
		return log_likelihood;
	}
}


arma::Col<double> Recapture_Likelihood_FLAT::get_ll_phi_components() {
	if (fresh_ll) {
		return ll_phi_components;
	} else {
		update_ll_phi_components();
		return ll_phi_components;
	}
}

arma::Col<double> Recapture_Likelihood_FLAT::get_ll_p_components() {
	if (fresh_ll) {
		return ll_p_components;
	} else {
		update_ll_p_components();
		return ll_p_components;
	}
}

double Recapture_Likelihood_FLAT::get_part_ll(
	arma::Col<arma::uword> indexes		
) {
	update_part_ll(indexes);
	return part_log_likelihood;
}


arma::Col<double> Recapture_Likelihood_FLAT::get_ll_phi_components(
	arma::Col<arma::uword> indexes
) {
	update_ll_phi_components(indexes);
	return ll_phi_components.elem(indexes);
}

arma::Col<double> Recapture_Likelihood_FLAT::get_ll_p_components(
	arma::Col<arma::uword> indexes
) {
	update_ll_p_components(indexes);
	return ll_p_components.elem(indexes);
}

// Likelihood calculations and getters, the meat.
void Recapture_Likelihood_FLAT::update_ll_phi_components() {
    for ( arma::uword i=0; i < number_of_individuals; ++i ) {
			ll_phi_components[i] = 0.0;
			for ( unsigned int t=tb[i]; t < td[i]-1; ++t ) {
				ll_phi_components[i] += log(PHI(i,t));
			}
			if (!known_death[i]) ll_phi_components[i] += log(1-PHI(i,td[i]-1));
    }
}

void Recapture_Likelihood_FLAT::update_ll_p_components() {
		for ( unsigned int i=0; i < number_of_individuals; ++i ) {
			ll_p_components[i] = 0.0;
			for ( unsigned int t=tb[i]+1; t < td[i]; ++t ) {
				if (sampled[t]) {
					if ( caught(i,t) == 1 ) {
						ll_p_components[i] += log(P(i,t));
					} else {
						ll_p_components[i] += log(1-P(i,t));
					}
				}
			}
		}
}

void Recapture_Likelihood_FLAT::update_ll() {
	update_ll_phi_components();
	update_ll_p_components();
	log_likelihood = arma::accu(ll_phi_components) + arma::accu(ll_p_components);
}


// Indexed likelihood getters, the partialy meat.
void Recapture_Likelihood_FLAT::update_ll_phi_components(
    arma::Col<arma::uword> indexes
) {
		int i=0;
    for ( int k=0; k < indexes.size(); ++k ) {
			i = indexes[k];
			ll_phi_components[i] = 0;
			for ( unsigned int t=tb[i]; t < td[i]-1; ++t ) {
				ll_phi_components[i] += log(PHI(i,t));
			}
			if (!known_death[i]) ll_phi_components[i] += log(1-PHI(i,td[i]-1));
    }
}



void Recapture_Likelihood_FLAT::update_ll_p_components(
	arma::Col<arma::uword> indexes		
) {
		int i=0;
    for ( unsigned int k=0; k < indexes.size(); ++k ) {
			i = indexes[k];
			ll_p_components[i] = 0.0;
			for ( unsigned int t=tb[i]+1; t < td[i]; ++t ) {
				if (sampled[t]) {
					if ( caught(i,t) == 1 ) {
						ll_p_components[i] += log(P(i,t));
					} else {
						ll_p_components[i] += log(1-P(i,t));
					}
				}
			}

    }

}


void Recapture_Likelihood_FLAT::update_part_ll( arma::Col<arma::uword> indexes ) {
	if (!fresh_ll) {
		update_ll_p_components(indexes);
		update_ll_phi_components(indexes);
	}
	part_log_likelihood = arma::accu(ll_phi_components.elem(indexes)) + 
												arma::accu(ll_p_components.elem(indexes));
}


void Recapture_Likelihood_FLAT::Recapture_Likelihood_FLAT::init() {
	resize_PHI(SES);
	PHI.zeros();
	P.zeros();
}



// Posterior functions. 
Recapture_Posterior_FLAT::Recapture_Posterior_FLAT() : 
	Recapture_Likelihood_FLAT() {}

Recapture_Posterior_FLAT::Recapture_Posterior_FLAT(
		std::vector<int> times_of_surveys,
		std::vector<std::vector<int> > times_of_recaptures,
		std::vector<int> times_of_deaths,
		std::vector<bool> known_deaths
) : Recapture_Likelihood_FLAT(times_of_surveys, times_of_recaptures,
			times_of_deaths, known_deaths) {
	init();
}

double Recapture_Posterior_FLAT::get_lp() { return 0; }   /// Flat priors....

double Recapture_Posterior_FLAT::get_log_posterior() {
	double ll = get_ll();
	double lp = get_lp();
	return ll + lp;
};

void Recapture_Posterior_FLAT::init() {}

// Sampler functions:
Recapture_Proposal_FLAT::Recapture_Proposal_FLAT() :
	Recapture_Posterior_FLAT() {}

Recapture_Proposal_FLAT::Recapture_Proposal_FLAT(
		std::vector<int> times_of_surveys,
		std::vector<std::vector<int> > times_of_recaptures,
		std::vector<int> times_of_deaths,
		std::vector<bool> known_deaths
) : Recapture_Posterior_FLAT(
		times_of_surveys, times_of_recaptures,
		times_of_deaths, known_deaths) 
{ 
	init();
}

double Recapture_Proposal_FLAT::propose_td() {
	double log_asymmetry = 0.0;
	for ( arma::uword i=0; i < number_of_individuals; ++i) {
		for ( int t=lo[i]; t < PHI.n_cols; ++t) {
			td[i] = t + 1;
			if ( U(R) > PHI[i,t] ) break; 
		}
	}
	return log_asymmetry;
}

double Recapture_Proposal_FLAT::propose_td( arma::Col<arma::uword> indexes ) {
	double log_asymmetry = 0.0;

	return log_asymmetry;
}

void Recapture_Proposal_FLAT::init() {
	for (arma::uword i=0; i < number_of_individuals; ++i) {
		
	}
}
