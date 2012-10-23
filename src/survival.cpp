#include "survival.hpp"

// Member functions for Data_FLAT only:
//
Recapture_Data_FLAT::Recapture_Data_FLAT() : number_of_individuals(0) {}

Recapture_Data_FLAT::Recapture_Data_FLAT(
    std::vector<int> times_of_surveys,
    std::vector<std::vector<int> > times_of_recaptures
) : ts(times_of_surveys), number_of_individuals(times_of_recaptures.size()),
        fo(times_of_recaptures.size()), lo(times_of_recaptures.size()),
        caught(times_of_recaptures.size(),times_of_surveys.size()),
				known_death(times_of_recaptures.size()),
        tb(times_of_recaptures.size()) {
    for ( int i=0; i < number_of_individuals; ++i ) {
        for ( int j=0; j < times_of_recaptures[i].size(); ++j ) {
            caught[i,times_of_recaptures[i][j]] = 1;
        }
    }
    init();
}


int Recapture_Data_FLAT::get_N() const { return number_of_individuals; }
arma::SpMat<int> Recapture_Data_FLAT::get_recaptures(int i) const { 
	return caught.row(i); 
}
int Recapture_Data_FLAT::get_tb(int i) const { return tb(i); }
arma::Col<int> Recapture_Data_FLAT::get_surveys() const { return ts; }

void Recapture_Data_FLAT::init() {
    bool find_fo = true;
    for ( int i=0; i < caught.n_rows; ++i ) {
        for ( int j=0; j < caught.n_cols; ++j ) {
            if (find_fo && caught(i,j) == 1 ) { fo[i] = j; find_fo = false; }
            if (           caught(i,j) == 1 ) { lo[i] = j;                  }
						caught_double(i,j) == (double)caught(i,j);
						if ( caught(i,j) == 0 ) uncaught(i,j) = 1;
						uncaught_double(i,j) == (double)caught(i,j);
        }
				known_death[i] = false;
    }
    tb = fo;
}

// Member functions for State_FLAT only:
//
Recapture_State_FLAT::Recapture_State_FLAT() :
    Recapture_Data_FLAT() { init(); }

Recapture_State_FLAT::Recapture_State_FLAT(
    std::vector<int> times_of_surveys,
    std::vector<std::vector<int> > times_of_recaptures
) : Recapture_Data_FLAT(times_of_surveys, times_of_recaptures),
        td(times_of_recaptures.size()),
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

void Recapture_State_FLAT::init() {
    td = lo;
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
    log_likelihood(0), SES(0), ONES(0,0), ZEROS(0,0), fresh_ll(false)
{
    init();
}

Recapture_Likelihood_FLAT::Recapture_Likelihood_FLAT(
	std::vector<int> times_of_surveys,
	std::vector<std::vector<int> > times_of_recaptures
) : Recapture_State_FLAT(times_of_surveys, times_of_recaptures),
    PHI(times_of_recaptures.size(),times_of_surveys.size()-1),
    P(times_of_recaptures.size(),times_of_surveys.size()),
    log_likelihood(0), 
		ll_phi_components(times_of_recaptures.size()), 
		ll_p_components(times_of_recaptures.size()), 
		SES(3), ONES(times_of_recaptures.size(), times_of_surveys.size()),
		ZEROS(times_of_recaptures.size(), times_of_surveys.size()),
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


double Recapture_Likelihood_FLAT::get_ll() { 
	if(fresh_ll) {
		return log_likelihood; 
	} else {
		update_ll();
		return log_likelihood;
	}
}

double Recapture_Likelihood_FLAT::get_part_ll(
	arma::Col<arma::uword> indexes		
) {
	update_ll_p_components(indexes);
	update_ll_phi_components(indexes);
	update_part_ll(indexes);
	return part_log_likelihood;

}

void Recapture_Likelihood_FLAT::update_ll_phi_components() {
    for ( arma::uword i=0; i < number_of_individuals; ++i ) {
			ll_phi_components[i] = 0.0;
			for ( unsigned int t=tb[i]; t < td[i]-1; ++t ) {
				ll_phi_components[i] += log(PHI(i,t));
			}
			if (!known_death[i]) ll_phi_components[i] += log(1-PHI[i,td[i]-1]);
    }
}


void Recapture_Likelihood_FLAT::update_ll_phi_components(
    arma::Col<arma::uword> indexes
) {
		int i=0;
    for ( int k=0; k < indexes.size(); ++k ) {
			i = indexes[k];
			ll_phi_components[i] = 0;
			for ( unsigned int t=tb[i]; t < td[i]-1; ++t ) {
				ll_phi_components[i] += log(PHI[i,t]);
			}
			if (!known_death[i]) ll_phi_components[i] += log(1-PHI[i,td[i]-1]);
    }
}

void Recapture_Likelihood_FLAT::update_ll_p_components() {
		for ( unsigned int i=0; i < number_of_individuals; ++i ) {
			ll_p_components[i] = 0.0;
			for ( unsigned int t=tb[i]+1; t < td[i]; ++t ) {
				if ( caught[i,t] == 1 ) {
					ll_p_components[i] += log(P[i,t]);
				} else {
					ll_p_components[i] += log(1-log(P[i,t]));
				}
			}
		}
}


void Recapture_Likelihood_FLAT::update_ll_p_components(
	arma::Col<arma::uword> indexes		
) {
	arma::Col<double> ONES_row = arma::ones< arma::Col<double> >(available.n_cols);
		int i=0;
    for ( unsigned int k=0; k < indexes.size(); ++k ) {
			i = indexes[k];
			ll_p_components[i] = 0.0;
			for ( unsigned int t=tb[i]+1; t < td[i]; ++t ) {
				if ( caught[i,t] == 1 ) {
					ll_p_components[i] += log(P[i,t]);
				} else {
					ll_p_components[i] += log(1-log(P[i,t]));
				}
			}

    }

}

void Recapture_Likelihood_FLAT::update_ll() {
	log_likelihood = arma::accu(ll_phi_components) + arma::accu(ll_p_components);
}

void Recapture_Likelihood_FLAT::update_part_ll( arma::Col<arma::uword> indexes ) {
}

void Recapture_Likelihood_FLAT::Recapture_Likelihood_FLAT::init() {
	resize_PHI(SES);
}






