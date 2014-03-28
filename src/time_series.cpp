#include "time_series.hpp"
#include "rv_constant.hpp"
#include "rv_uniform.hpp"
#include "rv_normal.hpp"
#include "rv_t_walk.hpp"
#include "rv_missing_t_walk_observed_normal.hpp"
#include "rv_missing_t_walk_observed_interval.hpp"
#include "rv_missing_t_walk.hpp"

Time_Series_Posterior::Time_Series_Posterior(
        Time_Series_State const & state_,
        Time_Series_Parameters const & parameters_,
        trng::yarn2 & R_
) : state(state_),
    parameters(parameters_),
    R(R_),
    distributions(locations_.size()),
    sample_order() { }

// Available distributions:

void Time_Series_Posterior::bind_constant_distribution  (int which) {
    if (distributions[which] != NULL)
        throw std::logic_error(distribution_already_bound(which, "constant"));

    sample_order[0].push_back(which);
    distributions[which] =
        std::unique_ptr<Random>(new RV_Constant(draws[which]));
}


void Time_Series_Posterior::bind_uniform_distribution (int which) {
    if (distributions[which] != NULL)
        throw std::logic_error(distribution_already_bound(which, "uniform"));

    sample_order[1].push_back(which);
    distributions[which] =
        std::unique_ptr<Random>(
            new RV_Uniform(draws[which], minima[which], maxima[which], R));

}

void Time_Series_Posterior::bind_ordered_uniform_distribution (int which) {
    if ( ( (which-1) < 0) || ((which + 1) == draws.size()) )
        throw std::logic_error(off_the_end(which, "ordered_uniform"));

    if (distributions[which] != NULL)
        throw std::logic_error(distribution_already_bound(which, "ordered_uniform"));

    sample_order[2].push_back(which);
    distributions[which] =
        std::unique_ptr<Random>(new RV_Uniform(
            draws[which], draws[which-1], draws[which+1], R));

}

void Time_Series_Posterior::bind_normal_distribution (int which, trng::yarn2 & R) {
    if (distributions[which] != NULL)
        throw std::logic_error(distribution_already_bound(which, "normal"));

    sample_order[1].push_back(which);
    distributions[which] =
        std::unique_ptr<Random>(new RV_Normal(
            draws[which], locations[which], obs_scales[which], R));

}


void Time_Series_Posterior::bind_t_walk_distribution_open (int which) {
    if ( ( (which-1) < 0) || ((which + 1) == draws.size()) )
        throw std::logic_error(off_the_end(which, "t_walk_open"));

    if (distributions[which] != NULL)
        throw std::logic_error(distribution_already_bound(which, "t_walk_open"));

    sample_order[2].push_back(which);
    distributions[which] =
        std::unique_ptr<Random>(new RV_t_walk(
            draws[which-1], draws[which], drift[which-1],
            tails[which-1], scales[which-1], R));

}

void Time_Series_Posterior::bind_t_walk_distribution_open_reverse (int which, trng::yarn2 & R) {
    if ( ( (which-1) < 0) || ((which + 1) == draws.size()) )
        throw std::logic_error(off_the_end(which, "t_walk_open_reverse"));

    if (distributions[which] != NULL)
        throw std::logic_error(distribution_already_bound(which, "t_walk_open_reverse"));

    sample_order[2].push_back(which);
    distributions[which] =
        std::unique_ptr<Random>(new RV_t_walk(
            draws[which+1], draws[which], drift[which],
            tails[which], scales[which], R));

}

void Time_Series_Posterior::bind_t_walk_distribution (int which, trng::yarn2 & R) {
    if ( ( (which-1) < 0) || ((which + 1) == draws.size()) )
        throw std::logic_error(off_the_end(which, "t_walk"));

    if (distributions[which] != NULL)
        throw std::logic_error(distribution_already_bound(which, "t_walk"));

    sample_order[2].push_back(which);
    distributions[which] =
        std::unique_ptr<Random>(new RV_Missing_t_walk(
            draws[which-1], draws[which], draws[which+1],
            drift[which-1], drift[which],
            tails[which-1], tails[which],
            scales[which-1], scales[which], R));

}

void Time_Series_Posterior::bind_t_walk_observed_normal_distribution (int which, trng::yarn2 & R) {
    if ( ( (which-1) < 0) || ((which + 1) == draws.size()) )
        throw std::logic_error(off_the_end(which, "t_walk_observed_normal"));

    if (distributions[which] != NULL)
        throw std::logic_error(distribution_already_bound(which, "t_walk_observed_normal"));

    sample_order[2].push_back(which);
    distributions[which] =
        std::unique_ptr<Random>(new RV_Missing_t_walk_observed_normal(
            draws[which-1], draws[which], draws[which+1],
            drift[which-1], drift[which],
            tails[which-1], tails[which],
            scales[which-1], scales[which],
            locations[which], obs_scales[which], R));

}

void Time_Series_Posterior::bind_t_walk_observed_interval_distribution (int which, trng::yarn2 & R) {
    if ( ( (which-1) < 0) || ((which + 1) == draws.size()) )
        throw std::logic_error(off_the_end(which, "t_walk_observed_interval"));

    if (distributions[which] != NULL)
        throw std::logic_error(distribution_already_bound(which, "t_walk_observed_interval"));

    sample_order[2].push_back(which);
    distributions[which] =
        std::unique_ptr<Random>(new RV_Missing_t_walk_observed_interval(
            draws[which-1], draws[which], draws[which+1],
            drift[which-1], drift[which],
            tails[which-1], tails[which],
            scales[which-1], scales[which],
            minima[which-1], maxima[which], R));

}

void Time_Series_Posterior::drop_distribution(int which) {
    if ( ( (which-1) < 0) || ((which + 1) == draws.size()) )
        throw std::logic_error(off_the_end(which, "NOT DELETING"));

    if (distributions[which] == NULL)
        throw std::logic_error(distribution_not_bound(which, "NOT DELETING"));

    // Reverse lookup in sample_order to remove the value from the
    // correct key...
    for (unsigned int i = 0; i < sample_order.size(); ++i) {
        for (unsigned int j = 0; j < sample_order[i].size(); ++i) {
            if (sample_order[i][j] == which)
                sample_order[i].erase(sample_order[i].begin()+j);
        }
    }
    distributions[which].reset(NULL);


}

void Time_Series_Posterior::draw() {
    // Modify:
    //  - some random objects will own multiple indexes... not yet...

    unsigned int which=0;
    for (unsigned int level=0; level < sample_order.size(); ++level) {
        for (unsigned int i = 0; i < sample_order[level].size(); ++i) {
            which = sample_order[level][i];
            if (distributions[which] == NULL)
                throw std::logic_error(distribution_not_bound(which, "NOT DELETING"));
            else
                distributions[which]->draw();
        }
    }
}

arma::vec Time_Series_Posterior::lpdf(arma::vec X) {
    arma::vec lpdfs(distributions.size());
    unsigned int which=0;
    for (unsigned int level=0; level < sample_order.size(); ++level) {
        for (unsigned int i = 0; i < sample_order[level].size(); ++i) {
            which = sample_order[level][i];
            if (distributions[which] == NULL)
                throw std::logic_error(distribution_not_bound(which,"NOT CALCULATING LPDF"));
            else
                lpdfs[which] = distributions[which]->lpdf(X[which]);

        }
    }
    return lpdfs;
}

Time_Series_Posterior::~Time_Series_Posterior() { distributions.clear(); }

static std::string Time_Series_Posterior::distribution_already_bound(int which, std::string distr) {
    std::stringstream msg;
    msg << "The index " << which << " (" << (which+1) << ")"
        " already has a distribution.  Not adding" << distr << ".\n";
    return msg.str();
}

static std::string Time_Series_Posterior::distribution_not_bound(int which, std::string action) {
    std::stringstream msg;
    msg << "The index " << which << " (" << (which+1) << ")"
        " lacks a distribution: " << action << "."
        "\n";
    return msg.str();
}

static std::string Time_Series_Posterior::off_the_end(int which, std::string distr) {
    std::stringstream msg;
    msg << "The location " << which << " (" << (which+1) << ")"
        " is dependent on an off-the-end index.  Not adding " << distr << ".\n";
    return msg.str();
}
