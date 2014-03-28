#ifndef LOCATION_H
#define LOCATION_H

#include <vector>
#include <map>
#include <memory>
#include <armadillo>
#include <trng/yarn2.hpp>

class Random;

class Time_Series_State { // State in that observed and unobserved sense.

public:
    Time_Series_State(); // required if used by reference...
    Time_Series_State(
        std::vector<double> times_,
        std::vector<double> y_at_times_,
        std::vector<double> minima_at_times_,
        std::vector<double> maxima_at_times_
    );

const arma::Col<double> & get_times() const;
const arma::Col<double> & get_x_at_times() const;
const arma::Col<double> & get_y_at_times() const;
const arma::Col<double> & get_minima_at_times() const;
const arma::Col<double> & get_maxima_at_times() const;

void set_x_at_times(arma::Col<double> x_at_times_);

private:
    arma::Col<double> times;
    arma::Col<double> x_at_times;
    arma::Col<double> y_at_times;
    arma::Col<double> minima_at_times;
    arma::Col<double> maxima_at_times;

};

class Time_Series_Parameters { // Parameters in that unobservable sense... (?)

public:
    Time_Series_Parameters();
    Time_Series_Parameters(
            Time_Series_State const & state_
    );

    const arma::Col<double> & get_drift() const;
    const arma::Col<double> & get_scales() const;
    const arma::Col<double> & get_tails() const;
    const arma::Col<double> & get_obs_scales() const;

    void set_drift(arma::Col<double> drift_);
    void set_scales(arma::Col<double> scales_);
    void set_tails(arma::Col<double> tails_);
    void set_obs_scales(arma::Col<double> obs_scales_);

private:
    Time_Series_State const & state;
    arma::Col<double> drift;
    arma::Col<double> scales;
    arma::Col<double> tails;
    arma::Col<double> obs_scales;

};


class Time_Series_Posterior {

public:
    Time_Series_Posterior(
        Time_Series_State const & state_,
        Time_Series_Parameters const & parameters_,
        trng::yarn2 & R_
    );

    // Available distributions:
    void bind_constant_distribution (unsigned int which);
    void bind_uniform_distribution  (
            unsigned int which, trng::yarn2 & R);
    void bind_ordered_uniform_distribution (
            unsigned int which, trng::yarn2 & R);
    void bind_normal_distribution (
            unsigned int which, trng::yarn2 & R);
    void bind_t_walk_distribution_open(
            unsigned int which, trng::yarn2 & R);
    void bind_t_walk_distribution_open_reverse(
            unsigned int which, trng::yarn2 & R);
    void bind_t_walk_distribution (
            unsigned int which, trng::yarn2 & R);
    void bind_t_walk_observed_normal_distribution (
            unsigned int which, trng::yarn2 & R);
    void bind_t_walk_observed_interval_distribution (
            unsigned int which, trng::yarn2 & R);

    // Drop distribution:
    void drop_distribution(unsigned int which);

    void draw();
    arma::vec lpdf(arma::vec X);

    ~Time_Series_Posterior();

private:
    Time_Series_State const & state;
    Time_Series_Parameters const & parameters;
    trng::yarn2 & R;

    std::vector<std::unique_ptr<Random> > distributions;
    std::map<unsigned int, std::vector<unsigned int> > sample_order;
    std::map<unsigned int, std::vector<unsigned int> > sample_order_bk;

};

#endif
