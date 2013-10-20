#include "rv_missing_t_walk_observed_normal.hpp"

#include <boost/math/special_functions/gamma.hpp>
#include <math.h>
#include <limits>

const double pi = boost::math::constants::pi<double>();


RV_Missing_t_walk_observed_normal::RV_Missing_t_walk_observed_normal(
		double const & x1_,
		double 			 & X,
		double const & x3_,
		double const & os1_,
		double const & os2_,
		double const & p1_,   // degrees of freedom 1
		double const & p2_,   // degrees of freedom 2
		double const & s1_,   // scale 1
		double const & s2_,   // scale 2
		double const & Xobs_, // observed location
		double const & so2_,
		trng::yarn2 & R_
) :	x1(x1_), x2(X), x3(x3_),
		os1(os1_), os2(os2_),
		p1(p1_), p2(p2_), s1(s1_), s2(s2_), 
		Xobs(Xobs_), so2(so2_),
		total_slice_length(0.0),
		companion(5,5), R(R_), 		
		EXPO(1.0)//,
{
	// Companion matrix for eigenvalue peak-finding.
	companion.zeros();
	companion(1,0) = 1.0;
	companion(2,1) = 1.0;
	companion(3,2) = 1.0;
	companion(4,3) = 1.0;

	// Reserve sizes, is this important?
	peaks.reserve(3);
	peak_bound_lr.reserve(3);
	intervals.reserve(2);
}

std::map<std::string, double> RV_Missing_t_walk_observed_normal::state() const {
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
	out["so2"] = so2;
	out["Xobs"] = Xobs;
	return out;
}

void RV_Missing_t_walk_observed_normal::jump(double X) { x2 = X; }

double RV_Missing_t_walk_observed_normal::draw() {  // Slice sampler w/ 1-3 peaks.
  ly = lpdf() - EXPO(R);
	find_slice();
	double ii = 0;
	while(true) {
		ii++; if (ii > 10000) { std::cout << "Too many steps." << std::endl; }
		x_new = choose();	
		ly_new = lpdf(x_new);
		if (ly_new >= ly) 
			break; 
		else 
			trim();
	}
	x2 = x_new;
	return x2; 
}

double RV_Missing_t_walk_observed_normal::lpdf(double X) {
	double lpdf;
	lpdf = 	boost::math::lgamma((p1+1.0)/2.0) - boost::math::lgamma(p1/2.0) -
		0.5 * log(p1*pi*pow(s1,2)) - 
		(p1+1.0)/2.0 * log(1.0 + (pow(X -(x1+os1),2))/(p1*pow(s1,2)) ) +
					boost::math::lgamma((p2+1.0)/2.0) - boost::math::lgamma(p2/2.0) -
		0.5 * log(p2*pi*pow(s2,2)) - 
		(p2+1.0)/2.0 * log(1.0 + (pow(x3-(X +os2),2))/(p2*pow(s2,2)) );
	lpdf += -0.5 * log(2*pi*pow(so2,2)) - (pow(X - Xobs,2) / (2*pow(so2,2)));
	return lpdf;
}

double RV_Missing_t_walk_observed_normal::lpdf() { return lpdf(x2); }

void RV_Missing_t_walk_observed_normal::find_slice() {
	find_peaks();
	std::vector<double>::iterator peaks_end = 
		std::remove_if(peaks.begin(), peaks.end(), 
				[ly](double x) {return lpdf(x) <= ly ? true : false; });

	// step out from peak.
	for (std::vector<double>::iterator i = peaks.begin(); 
				i != peaks_end; i++) {
		peak_bound_lr.push_back(step_out(i));
	}
	for (std::vector<double>::iterator i = peak_bound_lr.begin(); 
				i != peak_bound_lr.end(); i++) 
	{
		// Trim overlap
		bool fb = (i == peak_bound_lr.begin()   );
		bool lb = (i == (peak_bound_lr.end() - 1);
		if (!fb && (*(i-1)[1] > *i[0]) ) {
			double mid = (*(i-1)[1] + *i[0])/2.0;
			*(i-1)[1] = mid;
			*i[0] = mid;
		}
		if (!lb && (*(i+1)[0] < *i[1]) ) {
			double mid = (*(i+1)[0] + *i[1])/2.0;
			*(i+1)[0] = mid;
			*i[1] = mid;
		}

		// Calculate sub-slice to sub-slice distances:
		if (!fb) {
			intervals.push_back(*i[0] - *(i-1)[1]);
		}
		total_slice_length += *i[1] - *i[0];
	}

}


std::vector<double> RV_Missing_t_walk_observed_normal::step_out(
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
	return bounds;
}

void RV_Missing_t_walk_observed_normal::trim() {
	for (std::vector<double>::iterator i = peak_bound_lr.begin(); 
				i != peak_bound_lr.end(); i++) 
	{
		if (*i[0] < x_new && x_new < *i[1]) {
			if ( (x_new - *i[0]) < (*i[1] - x_new) )
				*i[0] = x_new;
			else
				*i[1] = x_new;
		}
	}

}

double RV_Missing_t_walk_observed_normal::choose() {
	double l = U(R) * total_slice_length + peak_bound_lr[0][0]; 
	for (std::vector<double>::iterator i = peak_bound_lr.begin(); 
				i != peak_bound_lr.end(); i++) 
	{
		if ( l < *i[1] )
			return l + *i[0];
		else 
			l = l - *i[1];
	}
}

void RV_Missing_t_walk_observed_normal::derivative_poly() {
	double A1 = 2*( (deltat__-xtp1) - (deltatm1+xtm1) );
  double A2 = ( 
		pow(xtp1-deltat__,2) + pt__*pow(sigmat__,2) +
    pow(xtm1+deltatm1,2) + ptm1*pow(sigmatm1,2) +
   -4*(deltatm1+xtm1)*(deltat__-xtp1) 
	);
  double A3 = (
		2*( (deltat__-xtp1)*(pow(xtm1+deltatm1,2)+ptm1*pow(sigmatm1,2)) - 
        (deltatm1+xtm1)*(pow(xtp1-deltat__,2)+pt__*pow(sigmat__,2)) 
		)
	);
  double A4 = ( 
		(pow(xtp1-deltat__,2) + pt__*pow(sigmat__,2)) * 
    (pow(xtm1+deltatm1,2) + ptm1*pow(sigmatm1,2)) 
	);

	double B1 = ( -1*(ptm1+1) - (pt__+1) );
  double B2 = ( (ptm1+1)*(xtm1+deltatm1+2*xtp1-2*deltat__) +
          (pt__+1)*(xtp1-deltat__+2*xtm1+2*deltatm1)  );
  double B3 = (
	  -1*(ptm1+1)*(
		pt__*pow(sigmat__,2) + pow(xtp1-deltat__,2) + 
			2*(xtm1+deltatm1)*(xtp1-deltat__)
		) -
		 1*(pt__+1)*(
		ptm1*pow(sigmatm1,2) + pow(xtm1+deltatm1,2) + 
			2*(xtp1-deltat__)*(xtm1+deltatm1)) 
	);
  double B4 = ( 
		(ptm1+1)*(xtm1+deltatm1)*
			(pt__*pow(sigmat__,2)+pow(xtp1-deltat__,2)) +
    (pt__+1)*(xtp1-deltat__)*
			(ptm1*pow(sigmatm1,2)+pow(xtm1+deltatm1,2)) 
	);

	companion(0,4) = A4*Xobs           + B4*so2;
	companion(1,4) = A3*Xobs - A4      + B3*so2;
	companion(2,4) = A2*Xobs - A3      + B2*so2;
	companion(3,4) = A1*Xobs - A2      + B1*so2;
	companion(4,4) =    Xobs - A1							 ;
}

void RV_Missing_t_walk_observed_normal::find_peaks() {
	derivative_poly();
	if (!arma::eig_gen(cx_eigval, cx_eigvec, companion)) 
		throw std::runtime_error("Failed eigenvalue decomposition.");
	arma::cx_vec::iterator re_eigval_end = 
		std::remove_if(cx_eigval.begin(), cx_eigval.end(), has_imaginary);
	std:sort(cx_eigval.begin(), re_eigval_end);
	for (arma::cx_vec::iterator i = cx_eigval.begin(); 
				i != re_eigval_end; i += i+2) {
// remove_if to get rid of duplicates.
	num_real_eigval = re_eigval_end - CX_EIGVAl.begin() - 1;
	if (num_real_eigval == 0) {
		throw std::runtime_error("No real eigenvalues available to find peaks.");
	} else if (num_real_eigval == 1) {
		peaks.push_back(arma::real(cx_eigval)[0]);
	} else if (num_real_eigval == 3) {
		peaks.push_back(arma::real(cx_eigval)[0]);
		peaks.push_back(arma::real(cx_eigval)[2]);
	} else if (num_real_eigval == 5) {
		peaks.push_back(arma::real(cx_eigval)[0]);
		peaks.push_back(arma::real(cx_eigval)[2]);
		peaks.push_back(arma::real(cx_eigval)[4]);
	} else {
		throw std::runtime_error("Number of real eigenvalues not equal to 1, 3, or 5.");
	}

}	

