#ifndef SURVIVAL_H
#define SURVIVAL_H

#include <armadillo>
#include <vector>

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
	arma::Col<int> get_surveys() const;
	int get_birth(int i) const;

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

private:
	void init();

};



class Recapture_State_FLAT : public Recapture_Data_FLAT {
	
public:
	Recapture_State_FLAT();
	Recapture_State_FLAT(
			std::vector<int> times_of_surveys,
			std::vector<std::vector<int> > times_of_recaptures
	);

	void set_td(
		arma::Col<arma::uword> indexes,
		arma::Col<int> times_of_deaths
	);

	void known_deaths( arma::Col<arma::uword> indexes );

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
			std::vector<std::vector<int> > times_of_recaptures
	);

	void resize_PHI( unsigned int scale);
	

protected:
	arma::Mat<double> PHI;
	arma::Mat<double> P;
	arma::Mat<double> ONES;
	arma::Mat<double> ZEROS;

	double log_likelihood;
	double part_log_likelihood;
	arma::Col<double> ll_phi_components;
	arma::Col<double> ll_p_components;

	double get_ll();
	double get_part_ll( arma::Col<arma::uword> indexes );
	void update_ll_phi_components();
	void update_ll_phi_components( arma::Col<arma::uword> indexes );
	void update_ll_p_components();
	void update_ll_p_components( arma::Col<arma::uword> indexes );
	void update_ll();
	void update_part_ll( arma::Col<arma::uword> indexes );


private:
	void init();
	unsigned int SES;
	bool fresh_ll;
	arma::Col<int> fresh_ll_p_components;
	arma::Col<int> fresh_ll_phi_components;

};



#endif
