#ifndef PROPOSALS_FLAT_H
#define PROPOSALS_FLAT_H


#include <vector>
#include <map>

#include <armadillo>
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>

class Simulation_Proposal_FLAT {

public:
	Simulation_Proposal_FLAT();
	Simulation_Proposal_FLAT( const Recapture_Posterior_FLAT& theta );

public:
	arma::Col<int> propose_td( 
			const Recapture_Posterior_FLAT& theta 
	);
	double get_pd() const;

	arma::Col<int> propose_td( 
			const Recapture_Posterior_FLAT& theta,
			arma::Col<arma::uword> indexes 
	);
	double get_pd( arma::Col<arma::uword> indexes) const;

	arma::Col<double> calc_log_proposal_density(
		const Recapture_Posterior_FLAT& theta
	);  // For current state.


protected:
	trng::yarn2 R;
	trng::uniform01_dist<double> U;
	arma::Col<int> td_proposed;

private:
	void init();
	arma::Col<double> log_proposal_density;
};

/////

class Slice_Proposal_FLAT {

public:
	Slice_Proposal_FLAT();
	Slice_Proposal_FLAT( const Recapture_Posterior_FLAT& theta );

public:
	arma::Col<int> propose_td( 
			const Recapture_Posterior_FLAT& theta 
	);
	double get_pd() const;

	arma::Col<double> calc_log_proposal_density(
		const Recapture_Posterior_FLAT& theta
	);  // For current state.


protected:
	trng::yarn2 R;
	trng::uniform01_dist<double> U;
	arma::Col<int> td_proposed;

private:
	void init();
	arma::Col<double> log_proposal_density;
};
#endif
