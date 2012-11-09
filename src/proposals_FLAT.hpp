#ifndef PROPOSALS_FLAT_H
#define PROPOSALS_FLAT_H


#include <vector>
#include <map>

#include <armadillo>
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>

#include "survival_FLAT.hpp"

class Simulation_Proposal_FLAT {

public:
	Simulation_Proposal_FLAT();
	Simulation_Proposal_FLAT( const Recapture_Posterior_FLAT& theta_ );

public:
	arma::Col<int> propose_td();
	double get_pd() const;

	arma::Col<int> propose_td( 
			arma::Col<arma::uword> indexes 
	);
	double get_pd( arma::Col<arma::uword> indexes) const;

	arma::Col<double> calc_log_proposal_density();  // For current state.


protected:
	trng::yarn2 R;
	trng::uniform01_dist<double> U;
	arma::Col<int> td_proposed;
	const Recapture_Posterior_FLAT& theta;

private:
	void init();
	arma::Col<double> log_proposal_density;
};

/////

class Slice_Proposal_FLAT {

public:
	Slice_Proposal_FLAT();
	Slice_Proposal_FLAT( const Recapture_Posterior_FLAT& theta_ );

public:
	arma::Col<int> propose_td();

protected:
	trng::yarn2 R;
	trng::uniform01_dist<double> U;
	arma::Col<int> td_proposed;
	const Recapture_Posterior_FLAT& theta;

private:
	void init();
	arma::Mat<double> S;
	arma::Mat<double> D;
	arma::Mat<double> td_pdf;
	void calc_td_pdf();
};



#endif
