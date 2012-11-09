#include "proposals_FLAT.hpp"
#include <iostream>
#include <math.h>

// Sampler functions:
Simulation_Proposal_FLAT::Simulation_Proposal_FLAT() {}

Simulation_Proposal_FLAT::Simulation_Proposal_FLAT(
	const Recapture_Posterior_FLAT& theta
) : {
	td_proposed.set_size(theta.number_of_individuals);
	td_proposed.zeros();
	td_proposed = theta.td;
	log_proposal_density.set_size(number_of_individuals);
	log_proposal_density.zeros();
	log_proposal_density = calc_log_proposal_density(theta);
}

arma::Col<int> Simulation_Proposal_FLAT::propose_td(
	const Recapture_Posterior_FLAT& theta		
) {
	for ( arma::uword i=0; i < theta.number_of_individuals; ++i) {
		log_proposal_density[i] = 0.0;
		for ( int t=theta.lo[i]; t < theta.PHI.n_cols; ++t) {
			td_proposed[i] = t + 1;
			if ( U(R) > theta.PHI(i,t) ) {
				break; 
				log_proposal_density[i] += log(1.0-theta.PHI(i,t));
			} else {
				log_proposal_density[i] += log(theta.PHI(i,t));
			}
		}
	}
	return td_proposed;
}

arma::Col<int> Simulation_Proposal_FLAT::propose_td( 
	const Recapture_Posterior_FLAT& theta,
	arma::Col<arma::uword> indexes 
) {
	// Non-indexed portions of log proposal densities are left
	// non-updated, so they have non-sense values.  User must deal
	// with this properly (use matching call to get_pd.
	arma::uword i;
	for ( arma::uword k=0; k < indexes.n_elem; ++k ) {
		i = indexes[k];
		log_proposal_density[i] = 0.0;
		for ( int t=theta.lo[i]; t < theta.PHI.n_cols; ++t) {
			td_proposed[i] = t + 1;
			if ( U(R) > theta.PHI(i,t) ) {
				break; 
				log_proposal_density[i] += log(1.0-theta.PHI(i,t));
			} else {
				log_proposal_density[i] += log(    theta.PHI(i,t));
			}
		}
	}
	return td_proposed;
}


double Simulation_Proposal_FLAT::get_pd() const {
	return arma::accu(log_proposal_density);
}

double Simulation_Proposal_FLAT::get_pd( arma::Col<arma::uword> indexes) const {
	return arma::accu(log_proposal_density.elem(indexes));

arma::Col<double> Simulation_Proposal_FLAT::calc_log_proposal_density(
	const Recapture_Posterior_FLAT& theta		
) {
	for (arma::uword i=0; i < theta.number_of_individuals; ++i) { 
		// td[i]-1 is the interval # of death.
		log_proposal_density[i] = log(1-theta.PHI(i,theta.td[i]-1));		
		for ( int t=theta.lo[i]; t < theta.td[i]-1; ++t ) {
			log_proposal_density[i] += log(theta.PHI(i,t));
		}
	}
	return log_proposal_density;
}


//////

// Sampler functions:
Slice_Proposal_FLAT::Slice_Proposal_FLAT() : {}

Slice_Proposal_FLAT::Slice_Proposal_FLAT(
	const Recapture_Posterior_FLAT& theta
) : {
	td_proposed.set_size(theta.number_of_individuals);
	td_proposed.zeros();
	td_proposed = theta.td;
	log_proposal_density.set_size(number_of_individuals);
	log_proposal_density.zeros();
	log_proposal_density = calc_log_proposal_density(theta);
}

arma::Col<int> Slice_Proposal_FLAT::propose_td(
	const Recapture_Posterior_FLAT& theta		
) {
	double h;
	int tmax;
	if (!theta.fresh_ll) theta.calc_td_pdf();
	for ( arma::uword i=0; i < theta.number_of_individuals; ++i) {
		h = U(R) * theta.td_pdf(i,theta.td[i]);
		tmax = theta.lo[i]+1;
		for( unsigned int t=theta.lo[i]+2; h < theta.td_pdf(i,t); ++t ) {
			tmax = t;
		}
		td_proposed[i] = int(U(R) * (double(tmax) + 1.0));
	}
	return td_proposed;
}



double Slice_Proposal_FLAT::get_pd() const {
	return arma::accu(log_proposal_density);
}

arma::Col<double> Slice_Proposal_FLAT::calc_log_proposal_density(
	const Recapture_Posterior_FLAT& theta		
) {
	for (arma::uword i=0; i < theta.number_of_individuals; ++i) { 
		// td[i]-1 is the interval # of death.
		log_proposal_density[i] = log(1-theta.PHI(i,theta.td[i]-1));		
		for ( int t=theta.lo[i]; t < theta.td[i]-1; ++t ) {
			log_proposal_density[i] += log(theta.PHI(i,t));
		}
	}
	return log_proposal_density;
}


/// To calculate horizontal extent of slice, for each i, I need
/// the density for death at t (therefore 2d).  Run this once
/// for everyone will let me reuse the calculations
arma::Mat<double> Slice_Proposal_FLAT::calc_log_proposal_density(
	const Recapture_Posterior_FLAT& theta		
) {
	log_proposal_density[i] = log(1-theta.PHI(i,td[i]-1));		
	for ( int t=theta.lo[i]; t < td[i]-1; ++t ) {
		log_proposal_density[i] += log(theta.PHI(i,t));
	}
	return log_proposal_density;
}
