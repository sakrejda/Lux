#include "proposals_FLAT.hpp"
#include <iostream>
#include <math.h>

// Sampler functions:
//Simulation_Proposal_FLAT::Simulation_Proposal_FLAT(
//		) : theta() {}

Simulation_Proposal_FLAT::Simulation_Proposal_FLAT(
	const Recapture_Posterior_FLAT& theta_
) : theta(theta_) {
	td_proposed.set_size(theta.number_of_individuals);
	td_proposed.zeros();
	td_proposed = theta.td;
	log_proposal_density.set_size(theta.number_of_individuals);
	log_proposal_density.zeros();
	log_proposal_density = calc_log_proposal_density();
}

arma::Col<int> Simulation_Proposal_FLAT::propose_td() {
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
}

arma::Col<double> Simulation_Proposal_FLAT::calc_log_proposal_density() {
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
//Slice_Proposal_FLAT::Slice_Proposal_FLAT() {}

Slice_Proposal_FLAT::Slice_Proposal_FLAT(
	const Recapture_Posterior_FLAT& theta_
) : theta(theta_) {
	td_proposed.set_size(theta.number_of_individuals);
	td_proposed.zeros();
	td_proposed = theta.td;
	S.set_size(theta.PHI.n_rows, theta.PHI.n_cols);
	S.zeros();
	D.set_size(theta.PHI.n_rows, theta.PHI.n_cols);
	D.zeros();
	td_pdf.set_size(theta.PHI.n_rows, theta.PHI.n_cols);
	td_pdf.zeros();
	calc_td_pdf();
}

arma::Col<int> Slice_Proposal_FLAT::propose_td() {
	double h;
	int tmax;
	if (!theta.fresh_ll) calc_td_pdf();
	for ( arma::uword i=0; i < theta.number_of_individuals; ++i) {
		std::cout << "i: " << i << ", td[i]: " << theta.td[i] << std::endl;
		h = U(R) * td_pdf(i,theta.td[i]);
		tmax = theta.lo[i]+1;
		for( unsigned int t=theta.lo[i]+2; h < td_pdf(i,t); ++t ) {
			std::cout << ", tmax: " << tmax << std::endl;
			tmax = t;
		}
		td_proposed[i] = int(U(R) * (double(tmax) + 1.0));
	}
	return td_proposed;
}


void Slice_Proposal_FLAT::calc_td_pdf() {
	std::cout << "Recalc td pdf!" << std::endl;
	for ( unsigned int i=0; i < theta.number_of_individuals; ++i ) {
		S(i,theta.lo[i]+1) = 0.0;
		D(i,theta.lo[i]+1) = log( 1-theta.PHI(i,theta.lo[i]) );
		td_pdf(i,theta.lo[i]+1) = S(i,theta.lo[i]+1) + D(i,theta.lo[i]+1);
		for ( unsigned int t=theta.lo[i]+2; t < theta.number_of_occasions; ++t ) {
			S(i,t) = log(   theta.PHI(i,t-2) ) + S(i,t-1);
			D(i,t) = log( 1-theta.PHI(i,t-1) );
			td_pdf(i,t) = S(i,t) + D(i,t);
		}
	}
}


