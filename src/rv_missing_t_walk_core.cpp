#include "rv_missing_t_walk_core.hpp"
#include <iostream>

RV_Missing_t_walk_core::RV_Missing_t_walk_core(
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
) : x1(x1_), x2(X), x3(x3_),
        os1(os1_), os2(os2_),
        p1(p1_), p2(p2_), s1(s1_), s2(s2_),
        ly(0.0), x_new(0.0), ly_new(0.0),
        R(R_), EXPO(1.0),
        total_slice_length(0.0)  { }

std::map<std::string, double> RV_Missing_t_walk_core::state() const {
    std::map<std::string, double> out;
    out["x1"] = x1;
    out["X"] = x2;
    out["x3"] = x3;
    out["os1"] = os1;
    out["os2"] = os2;
    out["p1"] = p1;
    out["p2"] = p2;
    out["s1"] = s1;
    out["s2"] = s2;
    return out;
}

void RV_Missing_t_walk_core::draw() {
  ly = lpdf() - EXPO(R);
    find_slice();

    while(true) {
        x_new = choose();
        if (lpdf(x_new) >= ly)
            break;
        else
            trim();
    }
    x2 = x_new;
}


void RV_Missing_t_walk_core::find_slice() {
    total_slice_length = 0.0;
    intervals.clear();
    find_peaks();

    peaks.erase(
        std::remove_if(peaks.begin(), peaks.end(),
            [=](double x) {
                return lpdf(x) < ly ? true : false;
            }
        ),
        peaks.end()
    );
    if (peaks.size() < 1)
        throw std::runtime_error("No peaks found above slice.");

    // step out from peak.
    peak_bound_lr.clear();
    for (std::vector<double>::iterator i = peaks.begin();
                i != peaks.end(); ++i) {
        peak_bound_lr.push_back(step_out(i));
    }
    // Fix overlap
    for (unsigned int i = 0; i != (peak_bound_lr.size()-1); i++) {
        if (peak_bound_lr[i+1][0] < peak_bound_lr[i][1] ) {
            double mid = (peaks[i+1] + peaks[i])/2.0;
            peak_bound_lr[i+1][0] = mid;
            peak_bound_lr[i][1] = mid;
        }
    }

    for (std::vector<std::vector<double> >::iterator i = peak_bound_lr.begin();
                i != peak_bound_lr.end(); ++i)
    {
        bool lb = (i == (peak_bound_lr.end() - 1 ));
        // Calculate sub-slice to sub-slice distances:
        total_slice_length += (*i)[1] - (*i)[0];
        if (!lb) {
            intervals.push_back((*(i+1))[0] - (*i)[1]);
        }
    }
    if (peak_bound_lr.size() != peaks.size() )
        throw std::runtime_error("Mismatch between number of peaks and number of bounds.");

}


std::vector<double> RV_Missing_t_walk_core::step_out(
        std::vector<double>::iterator peak_iter) {
    bool fp = (peak_iter == (peaks.begin()) );
    bool lp = (peak_iter == (peaks.end()-1) );
    std::vector<double> bounds(2);
    int m = 200;
    double w = (s1+s2)/2.0;
    bounds[0] = *peak_iter - w * U(R);  // w = use (s1+s2)/2
    bounds[1] = bounds[0] + w;
    int j = std::floor(m * U(R));    // m needed...
    int k = (m-1) - j;
    while ((j>0) && (ly < lpdf(bounds[0])) ) {
        bounds[0] = bounds[0] - w;
        j = j - 1;
        if (!fp && (bounds[0] < *(peak_iter-1))) break;
    }
    while ((k>0) && (ly < lpdf(bounds[1])) ) {
        bounds[1] = bounds[1] + w;
        k = k - 1;
        if (!lp && (bounds[1] > *(peak_iter+1))) break;
    }
    if (*peak_iter < bounds[0] || *peak_iter > bounds[1]) {
        throw std::runtime_error("Peak not in bounds.");
    }
    return bounds;
}

void RV_Missing_t_walk_core::trim() {
    for (unsigned int i = 0; i != peak_bound_lr.size(); i++) {
        if (peak_bound_lr[i][0] < x_new && x_new < peak_bound_lr[i][1]) {
            if (x_new < peaks[i] ) {
                total_slice_length -= x_new - peak_bound_lr[i][0];
                peak_bound_lr[i][0] = x_new;
            }   else {
                total_slice_length -= peak_bound_lr[i][1] - x_new;
                peak_bound_lr[i][1] = x_new;
            }
            if ( (peak_bound_lr[i][0] > peaks[i]) || (peak_bound_lr[i][1] < peaks[i]) )
                throw std::runtime_error("Peak not within bounds.");
        }
    }
}

double RV_Missing_t_walk_core::choose() {
    double l = U(R) * total_slice_length + peak_bound_lr[0][0];
    for (std::vector<std::vector<double> >::iterator i = peak_bound_lr.begin();
                i != peak_bound_lr.end(); ++i)
    {
        if ( l <= (*i)[1] ) {
            return l; // + (*i)[0];
        } else
            l = l - ( (*i)[1] - (*i)[0] ) - (*i)[0] + (*(i+1))[0];
    }
    std::runtime_error("RV_Missing_t_walk_core::choose ended out of bounds.");
}

void RV_Missing_t_walk_core::find_peaks() {
    derivative_poly();
    cx_eigval.clear();
    cx_eigvec.clear();
    if (!arma::eig_gen(cx_eigval, cx_eigvec, companion))
        throw std::runtime_error("Failed eigenvalue decomposition.");
    arma::cx_vec::iterator re_eigval_end =
        std::remove_if(cx_eigval.begin(), cx_eigval.end(),
            [](std::complex<double> x) { return x.imag() != 0;});
    std:sort(cx_eigval.begin(), re_eigval_end,
        [](std::complex<double> a, std::complex<double> b) {
            return a.real() < b.real();
        });
    re_eigval_end = std::unique(cx_eigval.begin(), re_eigval_end,
        [](std::complex<double> a, std::complex<double> b) {
            return abs(a.real()-b.real()) < pow(10,-3);
        });

    real_eigval.clear();
    std::transform(cx_eigval.begin(), re_eigval_end,
        std::back_inserter(real_eigval), [](std::complex<double> x) { return x.real();}
    );

    peaks.clear();
    valleys.clear();

    bool toggle = false;
  std::partition_copy(
        real_eigval.begin(), real_eigval.end(),
        std::back_inserter(peaks),
    std::back_inserter(valleys),
    [&toggle](double) { return toggle = !toggle; });

    if (peaks.size() < 1) {
        throw std::runtime_error("No peaks found.");
    }
}

