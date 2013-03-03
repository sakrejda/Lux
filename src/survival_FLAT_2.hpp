#ifndef SURVIVAL_H
#define SURVIVAL_H

#include <vector>
#include <map>

#include <armadillo>

class Recapture_Theta {

public:
	Recapture_Theta_FLAT();  // required
	Recapture_Theta_FLAT(
		std::vector<int> times_of_surveys,  // ts
		std::vector<std::vector<int> > times_of_recaptures   //  build 'caught'
	);

	int get_N() const;
	int get_K() const;
	arma::Col<int> get_recaptures(int i) const;
	arma::Col<int> get_surveys()   const;
	arma::Col<int> get_births()    const;
	arma::Col<int> get_first_obs() const; 
	arma::Col<int> get_last_obs()  const; 
	std::vector<bool> get_sampled();
	arma::Col<int> get_deaths()    const;

	void set_td(
		arma::Col<arma::uword> indexes,
		arma::Col<int> times_of_deaths
	);
	void known_deaths( arma::Col<arma::uword> indexes );

protected:
	arma::Col<int> ts;
	arma::Col<int> tb;
	arma::Col<int> fo;
	arma::Col<int> lo;
	arma::SpMat<int> caught;
	arma::SpMat<double> caught_double;
	arma::SpMat<int> uncaught;
	arma::SpMat<double> uncaught_double;
	arma::Col<int> td;
	arma::SpMat<double> available;
	int number_of_individuals;
	int number_of_occasions;
	std::vector<bool> known_death;
	std::map<unsigned int, bool> sampled;

private:
	void init();

};



#endif
