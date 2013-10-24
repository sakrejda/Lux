#include "rv_missing_t_walk_core.hpp"
#include <iostream>

RV_Missing_t_walk_core::RV_Missing_t_walk_core(
		double const & x1_,
		double 			 & X,
		double const & x3_,
		double const & os1_,
		double const & os2_,
		double const & p1_,   // degrees of freedom 1
		double const & p2_,   // degrees of freedom 2
		double const & s1_,   // scale 1
		double const & s2_,   // scale 2
		trng::yarn2 & R_
) :	x1(x1_), x2(X), x3(x3_),
		os1(os1_), os2(os2_),
		p1(p1_), p2(p2_), s1(s1_), s2(s2_), 
		total_slice_length(0.0),
		R(R_), EXPO(1.0) { }

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

void RV_Missing_t_walk_core::jump(double X) { x2 = X; }

double RV_Missing_t_walk_core::draw() {
  ly = lpdf() - EXPO(R);
	find_slice();
	double ii = 0;
	while(true) {
		ii++; if (ii > 10000) { 
			std::cout << "Too many steps." << std::endl; 
			for ( unsigned int i=0; i < peaks.size(); i++ ) {
				std::cout << "Peaks: " << peaks[i] << std::endl;
				std::cout << "LB: " << peak_bound_lr[i][0] << std::endl;
				std::cout << "RB: " << peak_bound_lr[i][1] << std::endl;
			}
		}

		x_new = choose();	
		std::cout << "BOOP did choosing." << std::endl;
		ly_new = lpdf(x_new);
		if (ly_new >= ly) 
			break; 
		else 
			trim();
	}
	std::cout << "BOOP" << std::endl;
	x2 = x_new;
	return x2; 
}


void RV_Missing_t_walk_core::find_slice() {
	find_peaks();
	std::vector<double>::iterator peaks_end = 
		std::remove_if(peaks.begin(), peaks.end(), 
				[=](double x) {return lpdf(x) <= ly ? true : false; });

	// step out from peak.
	for (std::vector<double>::iterator i = peaks.begin(); 
				i != peaks_end; i++) {
		peak_bound_lr.push_back(step_out(i));
		if (peak_bound_lr[peak_bound_lr.size()-1][0] > peaks[i] ||
				peak_bound_lr[peak_bound_lr.size()-1][1] < peaks[i]) {
					throw std::runtime_error("Peak not in bounds.");
		}
	}
	for (std::vector<std::vector<double> >::iterator i = peak_bound_lr.begin(); 
				i != peak_bound_lr.end(); i++) 
	{
		// Trim overlap
		bool fb = (i == (peak_bound_lr.begin()   ));
		bool lb = (i == (peak_bound_lr.end() - 1 ));
		if (!fb && ((*(i-1))[1] > (*i)[0]) ) {
			double mid = ((*(i-1))[1] + (*i)[0])/2.0;
			(*(i-1))[1] = mid;
			(*i)[0] = mid;
		}
		if (!lb && ((*(i+1))[0] < (*i)[1]) ) {
			double mid = ((*(i+1))[0] + (*i)[1])/2.0;
			(*(i+1))[0] = mid;
			(*i)[1] = mid;
		}

		// Calculate sub-slice to sub-slice distances:
		if (!fb) {
			intervals.push_back((*i)[0] - (*(i-1))[1]);
		}
		total_slice_length += (*i)[1] - (*i)[0];
	}

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
	return bounds;
}

void RV_Missing_t_walk_core::trim() {
	for (std::vector<std::vector<double> >::iterator i = peak_bound_lr.begin(); 
				i != peak_bound_lr.end(); i++) 
	{
		if ((*i)[0] < x_new && x_new < (*i)[1]) {
			if ( (x_new - (*i)[0]) < ((*i)[1] - x_new) )
				(*i)[0] = x_new;
			else
				(*i)[1] = x_new;
		}
	}

}

double RV_Missing_t_walk_core::choose() {    
	std::cout << "Total slice length: total_slice_length: ";
	std::cout << total_slice_length << std::endl;
	double l = U(R) * total_slice_length + peak_bound_lr[0][0]; 
	for (std::vector<std::vector<double> >::iterator i = peak_bound_lr.begin(); 
				i != peak_bound_lr.end(); i++) 
	{
			for ( unsigned int j=0; j < peaks.size(); j++ ) {
				std::cout << "Peaks: " << peaks[j] << std::endl;
				std::cout << "LB: " << peak_bound_lr[j][0] << std::endl;
				std::cout << "RB: " << peak_bound_lr[j][1] << std::endl;
			}
		std::cout << "l: " << l << std::endl;
		if ( l <= (*i)[1] )
			return l; // + (*i)[0];
		else 
			l = l - ( (*i)[1] - (*i)[0] ) + (*(i+1))[0];
	}
}

void RV_Missing_t_walk_core::find_peaks() {
	derivative_poly();
	if (!arma::eig_gen(cx_eigval, cx_eigvec, companion)) 
		throw std::runtime_error("Failed eigenvalue decomposition.");
	std::cout << cx_eigval << "\n";
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
	

	std::transform(cx_eigval.begin(), re_eigval_end, 
		std::back_inserter(real_eigval), [](std::complex<double> x) { return x.real();}
	);

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

