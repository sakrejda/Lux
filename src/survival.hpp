#ifndef SURVIVAL_H
#define SURVIVAL_H

#include <armadillo>
#include <vector>
#include <map>

class Recapture_Data_FLAT {

public:
	Recapture_Data_FLAT();  // required
	Recapture_Data_FLAT(
		std::vector<int> times_of_surveys,  // ts
		std::vector<std::vector<int> > times_of_recaptures   //  build 'caught'
	);

	int get_N() const;
	int get_K() const;
	arma::Row<int> get_recaptures(int i) const;
	arma::Row<int> get_surveys()   const;
	arma::Col<int> get_births()    const;
	arma::Col<int> get_first_obs() const; 
	arma::Col<int> get_last_obs()  const; 
	std::vector<bool> get_sampled();

protected:
	arma::Col<int> ts;
	arma::Col<int> tb;
	arma::Col<int> fo;
	arma::Col<int> lo;
	arma::SpMat<int> caught;
	arma::SpMat<double> caught_double;
	arma::SpMat<int> uncaught;
	arma::SpMat<double> uncaught_double;
	int number_of_individuals;
	int number_of_occasions;
	std::vector<bool> known_death;
	std::map<unsigned int, bool> sampled;

private:
	void init();

};



class Recapture_State_FLAT : public Recapture_Data_FLAT {
	
public:
	Recapture_State_FLAT();
	Recapture_State_FLAT(
			std::vector<int> times_of_surveys,
			std::vector<std::vector<int> > times_of_recaptures,
			std::vector<int> times_of_deaths,
			std::vector<bool> known_deaths
	);

	void set_td(
		arma::Col<arma::uword> indexes,
		arma::Col<int> times_of_deaths
	);

	void known_deaths( arma::Col<arma::uword> indexes );

	arma::Col<int> get_deaths()    const;

protected:
	arma::Col<int> td;
	arma::SpMat<double> available;

private:
	void init();

};

class Recapture_Likelihood_FLAT : public Recapture_State_FLAT {

public:
	Recapture_Likelihood_FLAT();
	Recapture_Likelihood_FLAT(
		std::vector<int> times_of_surveys,
		std::vector<std::vector<int> > times_of_recaptures,
		std::vector<int> times_of_deaths,
		std::vector<bool> known_deaths
	);

	void resize_PHI( unsigned int scale);
	
	arma::Mat<double> get_PHI();
	arma::Mat<double> get_P();

	void set_PHI( arma::Mat<double> );
	void set_P( arma::Mat<double> );

	double get_ll();
	arma::Col<double> get_ll_phi_components();
	arma::Col<double> get_ll_p_components();

	double get_part_ll( arma::Col<arma::uword> indexes );

protected:
	arma::Mat<double> PHI;
	arma::Mat<double> P;

	double log_likelihood;
	double part_log_likelihood;
	arma::Col<double> ll_phi_components;
	arma::Col<double> ll_p_components;

	void update_ll_phi_components();
	void update_ll_p_components();
	void update_ll();

	void update_ll_phi_components( arma::Col<arma::uword> indexes );
	void update_ll_p_components( arma::Col<arma::uword> indexes );
	void update_part_ll( arma::Col<arma::uword> indexes );


private:
	void init();
	unsigned int SES;
	bool fresh_ll;
	arma::Col<int> fresh_ll_p_components;
	arma::Col<int> fresh_ll_phi_components;

};



#endif
