#include "time_series.hpp"
#include "rv_constant.hpp"
#include "rv_uniform.hpp"
#include "rv_normal.hpp"
#include "rv_t_walk.hpp"
#include "rv_missing_t_walk_observed_normal.hpp"
#include "rv_missing_t_walk_observed_interval.hpp"
#include "rv_missing_t_walk.hpp"

#include <iostream>

Time_Series_Data::Time_Series_Data() :
    times(0.0), y_at_times(0.0),
    minima_at_times(0.0), maxima_at_times(0.0) {}

Time_Series_Data::Time_Series_Data(
    std::vector<double> times_,
    std::vector<double> y_at_times_,
    std::vector<double> minima_at_times_,
    std::vector<double> maxima_at_times_
) : times(times_),
    y_at_times(y_at_times_),
    minima_at_times(minima_at_times_),
    maxima_at_times(maxima_at_times_) { }

const arma::Col<double> & Time_Series_Data::get_times() const {return times;}
const arma::Col<double> & Time_Series_Data::get_y_at_times() const {return y_at_times;}
const arma::Col<double> & Time_Series_Data::get_minima_at_times() const {return minima_at_times;}
const arma::Col<double> & Time_Series_Data::get_maxima_at_times() const {return maxima_at_times;}


//Time_Series_Parameters::Time_Series_Parameters() :
//    data(), x_at_times(0.0), drift(0.0),
//    scales(0.0), tails(0.0), obs_scales(0.0) {}

Time_Series_Parameters::Time_Series_Parameters(
    Time_Series_Data & data_,
    std::vector<double> x_at_times_,
    std::vector<double> drift_,
    std::vector<double> scales_,
    std::vector<double> tails_,
    std::vector<double> obs_scales_
) : data(data_), x_at_times(x_at_times_), drift(drift_),
    scales(scales_), tails(tails_), obs_scales(obs_scales_) {}

const arma::Col<double> & Time_Series_Parameters::get_x_at_times() const {return x_at_times;}
const arma::Col<double> & Time_Series_Parameters::get_drift() const {return drift;}
const arma::Col<double> & Time_Series_Parameters::get_scales() const {return scales;}
const arma::Col<double> & Time_Series_Parameters::get_tails() const {return tails;}
const arma::Col<double> & Time_Series_Parameters::get_obs_scales() const {return obs_scales;}

void Time_Series_Parameters::set_x_at_times(arma::Col<double> x_at_times_) {x_at_times = x_at_times_; }
void Time_Series_Parameters::set_drift(arma::Col<double> drift_) {drift = drift_; }
void Time_Series_Parameters::set_scales(arma::Col<double> scales_) {scales = scales_; }
void Time_Series_Parameters::set_tails(arma::Col<double> tails_) {tails = tails_; }
void Time_Series_Parameters::set_obs_scales(arma::Col<double> obs_scales_) {obs_scales = obs_scales_; }

arma::Col<double> & Time_Series_Parameters::get_x_at_times_handle() {return x_at_times;}
arma::Col<double> & Time_Series_Parameters::get_drift_handle() {return drift;}
arma::Col<double> & Time_Series_Parameters::get_scales_handle() {return scales;}
arma::Col<double> & Time_Series_Parameters::get_tails_handle() {return tails;}
arma::Col<double> & Time_Series_Parameters::get_obs_scales_handle() {return obs_scales;}

int Time_Series_Parameters::dimension() const { return x_at_times.size(); }

Time_Series_Posterior::Time_Series_Posterior(
        Time_Series_Data const & data_,
        Time_Series_Parameters & parameters_,
        trng::yarn2 & R_
) : data(data_),
    parameters(parameters_),
    R(R_),
    times(data_.get_times()),
    y_at_times(data_.get_y_at_times()),
    minima_at_times(data_.get_minima_at_times()),
    maxima_at_times(data_.get_maxima_at_times()),
    x_at_times(parameters_.get_x_at_times_handle()),
    drift(parameters_.get_drift_handle()),
    scales(parameters_.get_scales_handle()),
    tails(parameters_.get_tails_handle()),
    obs_scales(parameters_.get_obs_scales_handle()),
    distributions(parameters_.dimension()),
    sample_order()
{
    std::for_each(distributions.begin(), distributions.end(),
        [](std::unique_ptr<Random> & p) { p = nullptr;});
}

// Available distributions:

void Time_Series_Posterior::bind_constant_distribution  (int which) {
    if (distributions[which])
        throw std::logic_error(distribution_already_bound(which, "constant"));

    sample_order[0].push_back(which);
    distributions[which] =
        std::unique_ptr<Random>(new RV_Constant(x_at_times[which]));
}


void Time_Series_Posterior::bind_uniform_distribution (int which) {
    if (distributions[which])
        throw std::logic_error(distribution_already_bound(which, "uniform"));

    sample_order[1].push_back(which);
    distributions[which] = std::unique_ptr<Random>(
        new RV_Uniform(
            x_at_times[which], minima_at_times[which], maxima_at_times[which], R)
    );

}

void Time_Series_Posterior::bind_ordered_uniform_distribution (int which) {
    if ( ( (which-1) < 0) || ((which + 1) == parameters.dimension()) )
        throw std::logic_error(off_the_end(which, "ordered_uniform"));

    if (distributions[which])
        throw std::logic_error(distribution_already_bound(which, "ordered_uniform"));

    sample_order[2].push_back(which);
    distributions[which] = std::unique_ptr<Random>(
        new RV_Uniform(
            x_at_times[which], x_at_times[which-1], x_at_times[which+1], R)
    );

}

void Time_Series_Posterior::bind_normal_distribution (int which) {
    if (distributions[which])
        throw std::logic_error(distribution_already_bound(which, "normal"));

    sample_order[1].push_back(which);
    distributions[which] = std::unique_ptr<Random>(
        new RV_Normal(
            x_at_times[which], y_at_times[which], obs_scales[which], R)
    );

}


void Time_Series_Posterior::bind_t_walk_distribution_open (int which) {
    if ( ( (which-1) < 0)  )
        throw std::logic_error(off_the_end(which, "t_walk_open"));

    if (distributions[which])
        throw std::logic_error(distribution_already_bound(which, "t_walk_open"));

    sample_order[2].push_back(which);
    distributions[which] = std::unique_ptr<Random>(
        new RV_t_walk(
            x_at_times[which-1], x_at_times[which], drift[which-1],
            tails[which-1], scales[which-1], R)
    );

}

void Time_Series_Posterior::bind_t_walk_distribution_open_reverse (int which) {
    if ( ((which + 1) >= parameters.dimension()) )
        throw std::logic_error(off_the_end(which, "t_walk_open_reverse"));

    if (distributions[which])
        throw std::logic_error(distribution_already_bound(which, "t_walk_open_reverse"));

    sample_order[2].push_back(which);
    distributions[which] = std::unique_ptr<Random>(
        new RV_t_walk(
            x_at_times[which+1], x_at_times[which], drift[which],
            tails[which], scales[which], R)
    );

}

void Time_Series_Posterior::bind_t_walk_distribution (int which) {
    if ( ( (which-1) < 0) || ((which + 1) == parameters.dimension()) )
        throw std::logic_error(off_the_end(which, "t_walk"));

    if (distributions[which])
        throw std::logic_error(distribution_already_bound(which, "t_walk"));

    sample_order[2].push_back(which);
    distributions[which] = std::unique_ptr<Random>(
        new RV_Missing_t_walk(
            x_at_times[which-1],    x_at_times[which],  x_at_times[which+1],
            drift[which-1],         drift[which],
            tails[which-1],         tails[which],
            scales[which-1],        scales[which], R)
    );

}

void Time_Series_Posterior::bind_t_walk_observed_normal_distribution (int which) {
    if ( ( (which-1) < 0) || ((which + 1) == parameters.dimension()) )
        throw std::logic_error(off_the_end(which, "t_walk_observed_normal"));

    if (distributions[which])
        throw std::logic_error(distribution_already_bound(which, "t_walk_observed_normal"));

    sample_order[2].push_back(which);
    distributions[which] = std::unique_ptr<Random>(
        new RV_Missing_t_walk_observed_normal(
            x_at_times[which-1],    x_at_times[which], x_at_times[which+1],
            drift[which-1],             drift[which],
            tails[which-1],             tails[which],
            scales[which-1],            scales[which],
            y_at_times[which],      obs_scales[which], R)
        );

}

void Time_Series_Posterior::bind_t_walk_observed_interval_distribution (int which) {
    if ( ( (which-1) < 0) || ((which + 1) == parameters.dimension()) )
        throw std::logic_error(off_the_end(which, "t_walk_observed_interval"));

    if (distributions[which])
        throw std::logic_error(distribution_already_bound(which, "t_walk_observed_interval"));

    sample_order[2].push_back(which);
    distributions[which] = std::unique_ptr<Random>(
                new RV_Missing_t_walk_observed_interval(
            x_at_times[which-1],            x_at_times[which], x_at_times[which+1],
            drift[which-1],                     drift[which],
            tails[which-1],                     tails[which],
            scales[which-1],                    scales[which],
            minima_at_times[which-1], maxima_at_times[which], R)
        );

}

void Time_Series_Posterior::drop_distribution(int which) {
    if ( ( (which) < 0) || (which >= parameters.dimension()) )
        throw std::logic_error(off_the_end(which, "NOT DELETING"));

    //std::for_each(distributions.begin(), distributions.end(),
        //[](std::unique_ptr<Random> & p) {
            //if (p == nullptr)
                //std::cout << "NULL" << std::endl;
            //else
                //std::cout << "NOT NULL" << std::endl;
        //});

    if (distributions[which] == nullptr)
        throw std::logic_error(distribution_not_bound(which, "NOT DELETING"));

    std::cout << "distributions[which] is NOT nullptr." << std::endl;
    // Reverse lookup in sample_order to remove the value from the
    // correct key...
    for (unsigned int i = 0; i < sample_order.size(); ++i) {
        for (unsigned int j = 0; j < sample_order[i].size(); ++i) {
            if (sample_order[i][j] == which)
                sample_order[i].erase(sample_order[i].begin()+j);
        }
    }
    distributions[which].reset(nullptr);


}

void Time_Series_Posterior::draw() {
    // Modify:
    //  - some random objects will own multiple indexes... not yet...

    unsigned int which=0;
    for (unsigned int level=0; level < sample_order.size(); ++level) {
        for (unsigned int i = 0; i < sample_order[level].size(); ++i) {
            which = sample_order[level][i];
            if (distributions[which] == nullptr)
                throw std::logic_error(distribution_not_bound(which, "NOT DRAWING"));
            else
                distributions[which]->draw();
        }
    }
}

arma::vec Time_Series_Posterior::lpdf(arma::vec X) {
    arma::vec lpdfs(parameters.dimension());
    unsigned int which=0;
    for (unsigned int level=0; level < sample_order.size(); ++level) {
        for (unsigned int i = 0; i < sample_order[level].size(); ++i) {
            which = sample_order[level][i];
            if (distributions[which] == nullptr)
                throw std::logic_error(distribution_not_bound(which,"NOT CALCULATING LPDF"));
            else
                lpdfs[which] = distributions[which]->lpdf(X[which]);

        }
    }
    if (sample_order.size() > 0)
        return lpdfs;
    else
        throw std::logic_error("NO distributions set, NOT calculating lpdf.");
}

Time_Series_Posterior::~Time_Series_Posterior() { distributions.clear(); }

std::string Time_Series_Posterior::distribution_already_bound(int which, std::string distr) {
    std::stringstream msg;
    msg << "The index " << which << " (" << (which+1) << ")"
        " already has a distribution.  Not adding " << distr << ".\n";
    return msg.str();
}

std::string Time_Series_Posterior::distribution_not_bound(int which, std::string action) {
    std::stringstream msg;
    msg << "The index " << which << " (" << (which+1) << ")"
        " lacks a distribution: " << action << "."
        "\n";
    return msg.str();
}

std::string Time_Series_Posterior::off_the_end(int which, std::string distr) {
    std::stringstream msg;
    msg << "The index " << which << " (" << (which+1) << ")"
        " is dependent on an off-the-end index: " << distr << ".\n";
    return msg.str();
}
