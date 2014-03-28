#ifndef TIME_SERIES_H
#define TIME_SERIES_H

#include <vector>
#include <map>
#include <memory>
#include <armadillo>
#include <trng/yarn2.hpp>

class Random;

class Time_Series_Data {

public:
    Time_Series_Data(); // required if used by reference...
    Time_Series_Data(
        std::vector<double> times_,
        std::vector<double> y_at_times_,
        std::vector<double> minima_at_times_,
        std::vector<double> maxima_at_times_
    );

    const arma::Col<double> & get_times() const;
    const arma::Col<double> & get_y_at_times() const;
    const arma::Col<double> & get_minima_at_times() const;
    const arma::Col<double> & get_maxima_at_times() const;

private:
    arma::Col<double> times;
    arma::Col<double> y_at_times;
    arma::Col<double> minima_at_times;
    arma::Col<double> maxima_at_times;
};


class Time_Series_Parameters { // Parameters in that unobservable sense...

public:
    Time_Series_Parameters();
    Time_Series_Parameters(
            Time_Series_Data & data_,
            arma::Col<double> x_at_times_,
            arma::Col<double> drift_;
            arma::Col<double> scales_;
            arma::Col<double> tails_;
            arma::Col<double> obs_scales_;
    );

    const arma::Col<double> & get_x_at_times() const;
    const arma::Col<double> & get_drift() const;
    const arma::Col<double> & get_scales() const;
    const arma::Col<double> & get_tails() const;
    const arma::Col<double> & get_obs_scales() const;

    void set_x_at_times(arma::Col<double> x_at_times_);
    void set_drift(arma::Col<double> drift_);
    void set_scales(arma::Col<double> scales_);
    void set_tails(arma::Col<double> tails_);
    void set_obs_scales(arma::Col<double> obs_scales_);

    arma::Col<double> & get_x_at_times_handle() ;
    arma::Col<double> & get_drift_handle() ;
    arma::Col<double> & get_scales_handle() ;
    arma::Col<double> & get_tails_handle() ;
    arma::Col<double> & get_obs_scales_handle() ;

    const int size() const;
private:
    const Time_Series_Data & data;
    arma::Col<double> x_at_times;
    arma::Col<double> drift;
    arma::Col<double> scales;
    arma::Col<double> tails;
    arma::Col<double> obs_scales;
};



class Time_Series_Posterior {

public:
    Time_Series_Posterior(
        Time_Series_Data const & data_,
        Time_Series_Parameters & parameters_,
        trng::yarn2 & R_
    );

    // Available distributions:
    void bind_constant_distribution (unsigned int which);
    void bind_uniform_distribution  (int which);
    void bind_ordered_uniform_distribution (int which);
    void bind_normal_distribution (int which);
    void bind_t_walk_distribution_open(int which);
    void bind_t_walk_distribution_open_reverse(int which);
    void bind_t_walk_distribution (int which);
    void bind_t_walk_observed_normal_distribution (int which);
    void bind_t_walk_observed_interval_distribution (int which);

    // Drop distribution:
    void drop_distribution(int which);

        // Main ops!
    void draw();
    arma::vec lpdf(arma::vec X);

        // Destructor needs to delete the smart pointers.
    ~Time_Series_Posterior();

private:
    // Const handles!
    arma::Col<double> const & times;
    arma::Col<double> const & y_at_times;
    arma::Col<double> const & minima_at_times;
    arma::Col<double> const & maxima_at_times;

    // Non-const handles.
    arma::Col<double> & x_at_times;
    arma::Col<double> & drift;
    arma::Col<double> & scales;
    arma::Col<double> & tails;
    arma::Col<double> & obs_scales;

    Time_Series_Data const & data;
    Time_Series_Parameters & parameters;
    trng::yarn2 & R;

    std::vector<std::unique_ptr<Random> > distributions;
    std::map<unsigned int, std::vector<unsigned int> > sample_order;
    std::map<unsigned int, std::vector<unsigned int> > sample_order_bk;

    std::string distribution_already_bound(int which, std::string distr);
    std::string distribution_not_bound(int which, std::string action);
    std::string off_the_end(int which, std::string distr);
};

#endif
