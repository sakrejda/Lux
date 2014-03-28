#include "time_series.hpp"
#include "rv_constant.hpp"
#include "rv_uniform.hpp"
#include "rv_normal.hpp"
#include "rv_t_walk.hpp"
#include "rv_missing_t_walk_observed_normal.hpp"
#include "rv_missing_t_walk_observed_interval.hpp"
#include "rv_missing_t_walk.hpp"

Time_Series::Time_Series(
        Time_Series_State const & state_,
        Time_Series_Parameters const & parameters_,
        trng::yarn2 & R_
) : locations(locations_),
        drift(drift_),
        tails(tails_),
        scales(scales_),
        obs_scales(obs_scales_),
        minima(minima_),
        maxima(maxima_),
        R(R_),
        distributions(locations_.size()),
        sample_order(),
        draws(draws_) { }

const arma::vec & Locations::state() const { return locations; }
const double & Locations::state(arma::uword which) const { return locations[which]; }

// Available distributions:
void Locations::bind_constant_distribution  (
        unsigned int which
) {
    if (distributions[which] == NULL) {
        sample_order[0].push_back(which);
        // Shouldn't this technically be locations not draws?
        // How to resolve?
        distributions[which] =
            std::unique_ptr<Random>(new RV_Constant(draws[which]));
    } else {
        std::stringstream msg;
        msg << "The location " << which << " (" << (which+1) << ")"
                     " already has a distribution.  Not adding.\n";
        throw(std::logic_error(msg.str()));
    }
}


void Locations::bind_uniform_distribution (
        unsigned int which, trng::yarn2 & R
) {
    if (distributions[which] == NULL) {
        sample_order[1].push_back(which);
        distributions[which] =
            std::unique_ptr<Random>(
                new RV_Uniform(draws[which], minima[which], maxima[which], R));
    } else {
        std::stringstream msg;
        msg << "The location " << which << " (" << (which+1) << ")"
                     " already has a distribution.  Not adding.\n";
        throw(std::logic_error(msg.str()));
    }
}

void Locations::bind_ordered_uniform_distribution (
        unsigned int which,
        trng::yarn2 & R
) {
    if ( ( (which-1) < 0) || ((which + 1) == draws.size()) ) {
        std::stringstream msg;
        msg << "The location " << which << " (" << (which+1) << ")"
                     " is dependent on an off-the-end index.\n";
        throw std::logic_error(msg.str());
    }
    if (distributions[which] == NULL) {
        sample_order[2].push_back(which);
        distributions[which] =
            std::unique_ptr<Random>(new RV_Uniform(
                draws[which], draws[which-1], draws[which+1], R));
    } else {
        std::stringstream msg;
        msg << "The location " << which << " (" << (which+1) << ")"
                     " already has a distribution.  Not adding.\n";
        throw(std::logic_error(msg.str()));
    }
}

void Locations::bind_normal_distribution (
        unsigned int which,
        trng::yarn2 & R
) {
    if (distributions[which] == NULL) {
        sample_order[1].push_back(which);
        distributions[which] =
            std::unique_ptr<Random>(new RV_Normal(
                draws[which], locations[which], obs_scales[which], R));
    } else {
        std::stringstream msg;
        msg << "The location " << which << " (" << (which+1) << ")"
                     " already has a distribution.  Not adding.\n";
        throw(std::logic_error(msg.str()));
    }
}

//void Locations::bind_t_walk_distribution (
//      unsigned int which,
//      double const & drift1, double const & drift2,
//      double const & p1, double const & p2,
//      double const & s1, double const & s2, trng::yarn2 & R
//) {
//  if (distributions[which] == NULL) {
//      sample_order[2].push_back(which);
//      distributions[which] =
//          std::unique_ptr<Random>(new RV_Missing_t_walk(
//              draws[which-1], draws[which], draws[which+1],
//              drift1, drift2, p1, p2, s1, s2, R));
//  } else {
//      std::stringstream msg;
//      msg << "The location " << which << " (" << (which+1) << ")"
//                   " already has a distribution.  Not adding.\n";
//      throw(std::logic_error(msg.str()));
//  }
//}

void Locations::bind_t_walk_distribution_open (
        unsigned int which, trng::yarn2 & R
) {
    if ((which-1) < 0 ) {
        std::stringstream msg;
        msg << "The location " << which << " (" << (which+1) << ")"
                     " is dependent on an off-the-end index.\n";
        throw std::logic_error(msg.str());
    }
    if (distributions[which] == NULL) {
        sample_order[2].push_back(which);
        distributions[which] =
            std::unique_ptr<Random>(new RV_t_walk(
                draws[which-1], draws[which], drift[which-1],
                tails[which-1], scales[which-1], R));
    } else {
        std::stringstream msg;
        msg << "The location " << which << " (" << (which+1) << ")"
                     " already has a distribution.  Not adding.\n";
        throw(std::logic_error(msg.str()));
    }
}

void Locations::bind_t_walk_distribution_open_reverse (
        unsigned int which, trng::yarn2 & R
) {
    if ((which) < 0  || ((which+1) >= draws.size()) ) {
        std::stringstream msg;
        msg << "The location " << which << " (" << (which+1) << ")"
                     " is dependent on an off-the-end index.\n";
        throw std::logic_error(msg.str());
    }
    if (distributions[which] == NULL) {
        sample_order[2].push_back(which);
        distributions[which] =
            std::unique_ptr<Random>(new RV_t_walk(
                draws[which+1], draws[which], drift[which],
                tails[which], scales[which], R));
    } else {
        std::stringstream msg;
        msg << "The location " << which << " (" << (which+1) << ")"
                     " already has a distribution.  Not adding.\n";
        throw(std::logic_error(msg.str()));
    }
}

void Locations::bind_t_walk_distribution (
        unsigned int which, trng::yarn2 & R
) {
    if ( ( (which-1) < 0) || ((which + 1) == draws.size()) ) {
        std::stringstream msg;
        msg << "The location " << which << " (" << (which+1) << ")"
                     " is dependent on an off-the-end index.\n";
        throw std::logic_error(msg.str());
    }
    if (distributions[which] == NULL) {
        sample_order[2].push_back(which);
        distributions[which] =
            std::unique_ptr<Random>(new RV_Missing_t_walk(
                draws[which-1], draws[which], draws[which+1],
                drift[which-1], drift[which],
                tails[which-1], tails[which],
                scales[which-1], scales[which], R));
    } else {
        std::stringstream msg;
        msg << "The location " << which << " (" << (which+1) << ")"
                     " already has a distribution.  Not adding.\n";
        throw(std::logic_error(msg.str()));
    }
}

void Locations::bind_t_walk_observed_normal_distribution (
        unsigned int which, trng::yarn2 & R
) {
    if ( ( (which-1) < 0) || ((which + 1) == draws.size()) ) {
        std::stringstream msg;
        msg << "The location " << which << " (" << (which+1) << ")"
                     " is dependent on an off-the-end index.\n";
        throw std::logic_error(msg.str());
    }
    if (distributions[which] == NULL) {
        sample_order[2].push_back(which);
        distributions[which] =
            std::unique_ptr<Random>(new RV_Missing_t_walk_observed_normal(
                draws[which-1], draws[which], draws[which+1],
                drift[which-1], drift[which],
                tails[which-1], tails[which],
                scales[which-1], scales[which],
                locations[which], obs_scales[which], R));
    } else {
        std::stringstream msg;
        msg << "The location " << which << " (" << (which+1) << ")"
                     " already has a distribution.  Not adding.\n";
        throw(std::logic_error(msg.str()));
    }
}

void Locations::bind_t_walk_observed_interval_distribution (
        unsigned int which, trng::yarn2 & R
) {
    if ( ( (which-1) < 0) || ((which + 1) == draws.size()) ) {
        std::stringstream msg;
        msg << "The location " << which << " (" << (which+1) << ")"
                     " is dependent on an off-the-end index.\n";
        throw std::logic_error(msg.str());
    }
    if (distributions[which] == NULL) {
        sample_order[2].push_back(which);
        distributions[which] =
            std::unique_ptr<Random>(new RV_Missing_t_walk_observed_interval(
                draws[which-1], draws[which], draws[which+1],
                drift[which-1], drift[which],
                tails[which-1], tails[which],
                scales[which-1], scales[which],
                minima[which-1], maxima[which], R));
    } else {
        std::stringstream msg;
        msg << "The location " << which << " (" << (which+1) << ")"
                     " already has a distribution.  Not adding.\n";
        throw(std::logic_error(msg.str()));
    }
}

void Locations::drop_distribution(unsigned int which) {
    if (distributions[which] == NULL) {
        std::stringstream msg;
        msg << "The location " << which << " (" << (which+1) << ")"
                     " does not have a distribution.  Not deleting.\n";
        throw(std::logic_error(msg.str()));
    } else {
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

}

void Time_Series_Posterior::draw() {
    // Modify:
    //  - random objects do not return on draws (done...)
    //  - some random objects will own multiple indexes... not yet...

    unsigned int which=0;
    for (unsigned int level=0; level < sample_order.size(); ++level) {
        for (unsigned int i = 0; i < sample_order[level].size(); ++i) {
            which = sample_order[level][i];
            if (distributions[which] == NULL) {
                std::stringstream msg;
                msg << "The location " << which << " (" << (which+1) << ")"
                             " lacks has a distribution.  Not drawing.\n";
                throw(std::logic_error(msg.str()));
            } else {
                distributions[which]->draw();
            }
        }
    }
}

arma::vec Time_Series_Posterior::lpdf(arma::vec X) {
    // Modify this loop to respect the sample_order vector...
    // Need to also add correct removal order.
    arma::vec lpdfs(distributions.size());
    unsigned int which=0;
    for (unsigned int level=0; level < sample_order.size(); ++level) {
        for (unsigned int i = 0; i < sample_order[level].size(); ++i) {
            which = sample_order[level][i];
            if (distributions[which] == NULL) {
                std::stringstream msg;
                msg << "The location " << which << " (" << (which+1) << ")"
                             " lacks has a distribution.  Not calculating.\n";
                throw(std::logic_error(msg.str()));
            } else {
                lpdfs[which] = distributions[which]->lpdf(X[which]);
            }
        }
    }
    return lpdfs;
}

Time_Series_Posterior::~Time_Series_Posterior() { distributions.clear(); }
