#ifndef RV_MISSING_T_WALK_CORE_H
#define RV_MISSING_T_WALK_CORE_H

#include "random.hpp"

#include <vector>
#include <map>
#include <algorithm>
#include <iterator>

#include <armadillo>
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>
#include <trng/exponential_dist.hpp>

// If I use the paradigm from 'recapture.hpp' where the state/parameters
// are classes, maybe this could be a template for a generic 1-D slice
// sampler core... ?  The peak/root finding would have to be broken out
// as well (?) to allow for either eigenvalue based or numerical root
// finding for the slice sampler starting pts...

class RV_Missing_t_walk_core : public Random {

public:
    RV_Missing_t_walk_core(
        double const & x1_,
        double           & X,
        double const & x3_,
        double const & os1_,
        double const & os2_,
        double const & p1_,   // degrees of freedom 1
        double const & p2_,   // degrees of freedom 2
        double const & s1_,   // scale 1
        double const & s2_,   // scale 2
        trng::yarn2 & R_
    );

    virtual std::map<std::string, double> state() const;
    void draw();

protected:
    virtual void derivative_poly() = 0;

    void find_slice();
    void find_peaks();
    std::vector<double> step_out(std::vector<double>::iterator peak_iter);
    double choose();
    void trim();

    double const & x1;
    double       & x2;
    double const & x3;
    double const & os1;
    double const & os2;
    double const & p1;
    double const & p2;
    double const & s1;
    double const & s2;

    double ly; // cached log pdf at current value...
    double x_new; // temporary new sampled value...
    double ly_new; // temporary pdf at sampled value...

    trng::yarn2 & R;
    trng::uniform01_dist<double> U;
    trng::exponential_dist<double> EXPO;

    arma::Mat<double> companion;
    arma::cx_vec cx_eigval;
    arma::cx_mat cx_eigvec;
    std::vector<double> real_eigval;

    std::vector<double> peaks;
    std::vector<double> valleys;
    std::vector<std::vector<double> > peak_bound_lr;
    std::vector<double> intervals;
    std::vector<std::map<std::string, double> > peakery;
  double total_slice_length;

};

#endif
