#include "proposals_FLAT.hpp"
#include <iostream>
#include <math.h>

// Sampler functions:
//Simulation_td_Proposal_FLAT::Simulation_td_Proposal_FLAT(
//		) : theta() {}

Simulation_td_Proposal_FLAT::Simulation_td_Proposal_FLAT(
	const Recapture_Posterior_FLAT& theta_
) : theta(theta_) {
	td_proposed.set_size(theta.PHI.n_rows);
	td_proposed.zeros();
	td_proposed = theta.td;
	log_proposal_density.set_size(theta.PHI.n_rows);
	log_proposal_density.zeros();
	log_proposal_density = calc_log_proposal_density();
}

arma::Col<int> Simulation_td_Proposal_FLAT::propose_td() {
	for ( arma::uword i=0; i < theta.PHI.n_rows; ++i) {
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

arma::Col<int> Simulation_td_Proposal_FLAT::propose_td( 
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


double Simulation_td_Proposal_FLAT::get_pd() const {
	return arma::accu(log_proposal_density);
}

double Simulation_td_Proposal_FLAT::get_pd( arma::Col<arma::uword> indexes) const {
	return arma::accu(log_proposal_density.elem(indexes));
}

arma::Col<double> Simulation_td_Proposal_FLAT::calc_log_proposal_density() {
	for (arma::uword i=0; i < theta.PHI.n_rows; ++i) { 
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
//Slice_td_Proposal_FLAT::Slice_td_Proposal_FLAT() {}

Slice_td_Proposal_FLAT::Slice_td_Proposal_FLAT(
	const Recapture_Posterior_FLAT& theta_
) : theta(theta_) {
	td_proposed.set_size(theta.PHI.n_rows);
	td_proposed.zeros();
	td_proposed = theta.td;
	S.set_size(theta.PHI.n_rows, theta.PHI.n_cols);
	S.zeros();
	D.set_size(theta.PHI.n_rows, theta.PHI.n_cols);
	D.zeros();
	td_pdf.set_size(theta.PHI.n_rows, theta.PHI.n_cols);
	td_pdf.zeros();
	for (unsigned int i=0; i < td_pdf.n_rows; ++i) {
		CH.push_back(new trng::discrete_dist(int(td_pdf.n_cols)));
		for (unsigned int t=0; t <= theta.lo[i]; ++t) {
			CH[i]->param(t,0.0);
		}
	}
	calc_td_pdf();
}

Slice_td_Proposal_FLAT::~Slice_td_Proposal_FLAT() {
	for (unsigned int i=0; i < td_pdf.n_rows; ++i) {
		delete CH[i];
	}
}

arma::Col<int> Slice_td_Proposal_FLAT::propose_td() {
	if (!theta.fresh_ll) calc_td_pdf();
	for ( arma::uword i=0; i < theta.PHI.n_rows; ++i) {
		td_proposed[i] = (*CH[i])(R);
		std::cout << td_proposed[i] << ", " << U(R) << std::endl;
	}
	return td_proposed;
}
//arma::Col<int> Slice_td_Proposal_FLAT::propose_td() {
//	double h;
//	int tmin, tmax;
//	if (!theta.fresh_ll) calc_td_pdf();
//	for ( arma::uword i=0; i < theta.PHI.n_rows; ++i) {
//		h = U(R) * exp(td_pdf(i,theta.td[i]));
//		tmin = theta.lo[i]+1;
//		tmax = 0;
//		for( unsigned int t=1; 
//				h < exp(td_pdf(i,t+tmin)) && (t+tmin) < theta.PHI.n_cols; ++t ) {
//			tmax = t;
//		}
//		td_proposed[i] = tmin + int(U(R) * (double(tmax) + 1.0));
//	}
//	return td_proposed;
//}


void Slice_td_Proposal_FLAT::calc_td_pdf() {
	for ( unsigned int i=0; i < theta.PHI.n_rows; ++i ) {
		S(i,theta.lo[i]+1) = 0.0;
		D(i,theta.lo[i]+1) = log( 1-theta.PHI(i,theta.lo[i]) );
		td_pdf(i,theta.lo[i]+1) = S(i,theta.lo[i]+1) + D(i,theta.lo[i]+1);
		for ( unsigned int t=theta.lo[i]+2; t < theta.PHI.n_cols; ++t ) {
			if ( t > theta.P.n_cols ) {
				S(i,t) = log(   theta.PHI(i,t-2) ) + S(i,t-1);
			} else {
				S(i,t) = log(   theta.PHI(i,t-2) ) + S(i,t-1) + log( 1 - theta.P[i,t-1]);
			}
			D(i,t) = log( 1-theta.PHI(i,t-1) );
			td_pdf(i,t) = S(i,t) + D(i,t);
			CH[i]->param(t,td_pdf(i,t));
		}
	}
}


