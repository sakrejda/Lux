#include "random.hpp"
//#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
#include <iostream>
#include <math.h>

//#include <slicer-continuous.hpp>

#include <limits>
//#include <Eigen/Dense>
//#include <Eigen/Eigenvalues>

//using boost::math::lgamma;
const double pi = boost::math::constants::pi<double>();

RV_Constant::RV_Constant(double & X) : x(X) {}

void RV_Constant::jump(double X) { }

double RV_Constant::draw() {return x;}

double RV_Constant::lpdf(double X) {
	if (X == x)
		return 0.0;
	else 
		return -1.0 * std::numeric_limits<double>::infinity();
}

double RV_Constant::lpdf() { return lpdf(x); }



RV_Uniform::RV_Uniform(
		double & X, 
		double const & minimum, double const & maximum, 
		trng::yarn2 & R_
) : x(X), min(minimum), max(maximum), R(R_) { }

void RV_Uniform::jump(double X) {
	if ((X > min) && (X < max))
		x = X;
}

double RV_Uniform::draw() { 
	x = min + (U(R) * (max - min)); 
	return x;
}

double RV_Uniform::lpdf(double X) {
	if ((X > min) && (X < max)) 
		return (1/(max-min));
	else
		return 0.0;
}

double RV_Uniform::lpdf() { 
	if ((x > min) && (x < max))
		return (1/(max-min)); 
	else
		return 0;
}

RV_Normal::RV_Normal(
		double       & X_,
		double const & mu_,
		double const & s_,
		trng::yarn2  & R_
) : X(X_), mu(mu_), s(s_), R(R_), NORMAL(mu_, s_) { }

void RV_Normal::jump(double X) { X = X; }

double RV_Normal::draw() {
	X = NORMAL(R);
	return X;
}

double RV_Normal::lpdf(double X) {
	NORMAL.mu(mu);
	NORMAL.sigma(s);
	return log(NORMAL.pdf(X));
}

double RV_Normal::lpdf() { lpdf(X); }

RV_t_walk::RV_t_walk(
		double const & x1_,
		double			 & X,
		double const & os_,
		double const & p1_,
		double const & s1_,
		trng::yarn2  & R_
) : x1(x1_), x2(X), os(os_), p1(p1_), s1(s1_), 
		R(R_), EXPO(1.0), bounds(2) { }

std::map<std::string, double> RV_t_walk::state() const {
	std::map<std::string, double> out;
	out["x1"] = x1;
	out["X"] = x2;
  out["os"] = os,
	out["p1"] = p1;
	out["s1"] = s1;
	return out;
}

void RV_t_walk::jump(double X) { 
	x2 = X; 
}

double RV_t_walk::draw() {
	double S, V, W, T;
	do {
		S = 2.0 * U(R) - 1.0;
		V = 2.0 * U(R) - 1.0;
		W = S*S + V*V;
	} while (W > 1.0);
	T = S * sqrt(p1*(pow(W,-2.0/p1)-1.0)/W);
  x2 = T * s1 + x1 + os;
	return x2; 
}

double RV_t_walk::lpdf(double X) {
	double lpdf;
	lpdf = boost::math::lgamma((p1+1.0)/2.0) - boost::math::lgamma(p1/2.0) -
		0.5 * log(p1*pi*pow(s1,2)) -
		(p1+1.0)/2.0 * log(1.0 + (pow(X-(x1+os),2))/(p1*pow(s1,2)) );
		
//		lgamma((p1+1.0)/2.0) - lgamma(p1/2.0) -
//		0.5 * log(p1*pi*pow(s1,2)) - 
//		(p1+1.0)/2.0 * log(1.0 + (pow(X-x1,2))/(p1*pow(s1,2)) ) +
//					lgamma((p2+1.0)/2.0) - lgamma(p2/2.0) -
//		0.5 * log(p2*pi*pow(s2,2)) - 
//		(p2+1.0)/2.0 * log(1.0 + (pow(x3-X,2))/(p2*pow(s2,2)) );
	return lpdf;
}

double RV_t_walk::lpdf() { return lpdf(x2); }


RV_Missing_t_walk::RV_Missing_t_walk(
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
		companion(3,3), R(R_), eigvalues(3),
		bounds_pk1(2), bounds_pk2(2), EXPO(1.0)//,
//		ecompanion(3,3)
{
	companion.zeros();
	companion(1,0) = 1.0;
	companion(2,1) = 1.0;
// 	ecompanion(0,0) = 0.0; ecompanion(0,1) = 0.0;
//	ecompanion(1,0) = 1.0; ecompanion(1,1) = 0.0;
//	ecompanion(2,0) = 0.0; ecompanion(2,1) = 1.0;
}

std::map<std::string, double> RV_Missing_t_walk::state() const {
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

void RV_Missing_t_walk::jump(double X) { x2 = X; }

double RV_Missing_t_walk::draw() {  // Slice sampler w/ two peaks.
  ly = lpdf() - EXPO(R);
	find_slice();
	double ii = 0;
	while(true) {
		//ii++; // if (ii > 10000) throw std::runtime_error("Too many steps.");
		ii++; if (ii > 10000) { std::cout << "Too many steps." << std::endl; }
		x_new = choose();	
		ly_new = lpdf(x_new);
		if (ly_new >= ly) // && (x_new > min) && (x_new < max))  truncation.
			break; 
		else 
			trim();
	}
	x2 = x_new;
	return x2; 
}

double RV_Missing_t_walk::lpdf(double X) {
	double lpdf;
	lpdf = 	boost::math::lgamma((p1+1.0)/2.0) - boost::math::lgamma(p1/2.0) -
		0.5 * log(p1*pi*pow(s1,2)) - 
		(p1+1.0)/2.0 * log(1.0 + (pow(X -(x1+os1),2))/(p1*pow(s1,2)) ) +
					boost::math::lgamma((p2+1.0)/2.0) - boost::math::lgamma(p2/2.0) -
		0.5 * log(p2*pi*pow(s2,2)) - 
		(p2+1.0)/2.0 * log(1.0 + (pow(x3-(X +os2),2))/(p2*pow(s2,2)) );
	return lpdf;
}

double RV_Missing_t_walk::lpdf() { return lpdf(x2); }

void RV_Missing_t_walk::find_slice() {
	find_peaks();
	if (ly <= lpdf(peak1)) bounds_pk1 = step_out(peak1);
	if (ly <= lpdf(peak2)) bounds_pk2 = step_out(peak2);
	if (ly > lpdf(peak1)) {   // In case peak 1 does not reach slice level.
		bounds_pk1[0] = bounds_pk2[0];
		bounds_pk1[1] = bounds_pk2[0];
	}
	if (ly > lpdf(peak2)) {		// In case peak 2 does not reach slice level.
		bounds_pk2[0] = bounds_pk1[0];
		bounds_pk2[1] = bounds_pk1[0];
	}
	if (  // In case the two parts of the slice overlap, split it:
		(bounds_pk1[0] < std::min(bounds_pk1[1], bounds_pk2[1])) &&
		(bounds_pk2[0] < std::min(bounds_pk1[1], bounds_pk2[1]))
	) {
		double x = (bounds_pk1[0]+bounds_pk2[1])/2.0;
		bounds_pk1[1] = x;
		bounds_pk2[0] = x;
	}
}


std::vector<double> RV_Missing_t_walk::step_out(double peak) {
	std::vector<double> bounds(2);
	int m = 200;
	double w = (s1+s2)/2.0;
	bounds[0] = peak - w * U(R);  // w = use (s1+s2)/2
	bounds[1] = bounds[0] + w;
	int j = std::floor(m * U(R));    // m needed...
	int k = (m-1) - j;
	while ((j>0) && (ly < lpdf(bounds[0])) ) { // trunc. can be added here.
		bounds[0] = bounds[0] - w;
		j = j - 1;
	}
	while ((k>0) && (ly < lpdf(bounds[1])) ) { // truncation can be added here.
		bounds[1] = bounds[1] + w;
		k = k - 1;
	}
	return bounds;
}

void RV_Missing_t_walk::trim() {
	if (x_new <= bounds_pk1[1]) {
		if (x_new <= peak1) {
			bounds_pk1[0] = x_new;
		} else {
			bounds_pk1[1] = x_new;
		}
	} else {
		if (x_new <= peak2) {
			bounds_pk2[0] = x_new;
		} else {
			bounds_pk2[1] = x_new;
		}
	}
}

double RV_Missing_t_walk::choose() {
	double d = (bounds_pk2[1] - bounds_pk2[0]) + (bounds_pk1[1] - bounds_pk1[0]);
	double q = (bounds_pk2[0] - bounds_pk1[1]);
	double l = U(R) * d;
	if ((bounds_pk1[0] + l) > bounds_pk1[1]) {
		return (l + q + bounds_pk1[0]);
	} else {
		return (l + bounds_pk1[0]);
	}
}

void RV_Missing_t_walk::find_peaks() {  
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

  companion(0,2) = -B4/B1;
	companion(1,2) = -B3/B1;
	companion(2,2) = -B2/B1;
	
	if (!arma::eig_gen(cx_eigval, cx_eigvec, companion)) 
		throw std::runtime_error("Failed eigenvalue decomposition.");
	arma::cx_vec::iterator re_eigval_end = 
		std::remove_if(cx_eigval.begin(), cx_eigval.end(), has_imaginary);
	std:sort(cx_eigval.begin(), re_eigval_end);
// ADD remove_if it's too close to the nearest value.
	num_real_eigval = re_eigval_end - cx_eigval.begin() - 1;
	if (num_real_eigval == 0) {
		throw std::runtime_error("No real eigenvalues available to find peaks.");
	} else if (num_real_eigval == 1) {
		peak1 = arma::real(cx_eigval)[0];
		peak2 = arma::real(cx_eigval)[0];
	} else if (num_real_eigval == 3) {
		peak1 = arma::real(cx_eigval)[0];
		peak2 = arma::real(cx_eigval)[2];
	} else {
		throw std::runtime_error("Number of real eigenvalues not equal to 0, 1, or 3.");
	}
//	ecompanion(0,2) = companion(0,2);
//	ecompanion(1,2) = companion(1,2);
//	ecompanion(2,2) = companion(2,2);
//	Eigen::EigenSolver<Eigen::MatrixXd> es(ecompanion);
//	std::cout << "Peaks: \n" << es.eigenvalues() << std::endl;
	
}

// observed_normal after here.
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
//		ecompanion(5,5)
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


// 	ecompanion(0,0) = 0.0; ecompanion(0,1) = 0.0;
//	ecompanion(1,0) = 1.0; ecompanion(1,1) = 0.0;
//	ecompanion(2,0) = 0.0; ecompanion(2,1) = 1.0;
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
		//ii++; // if (ii > 10000) throw std::runtime_error("Too many steps.");
		ii++; if (ii > 10000) { std::cout << "Too many steps." << std::endl; }
		x_new = choose();	
		ly_new = lpdf(x_new);
		if (ly_new >= ly) // && (x_new > min) && (x_new < max))  truncation.
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
		bool fb = (i == peak_bound_lr.begin());
		bool lb = (i == peak_bound_lr.end()  );
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
	bool fp = (peak_iter == peaks.begin());
	bool lp = (peak_iter == peaks.end());
	std::vector<double> bounds(2);
	int m = 200;
	double w = (s1+s2)/2.0;
	bounds[0] = *peak_iter - w * U(R);  // w = use (s1+s2)/2
	bounds[1] = bounds[0] + w;
	int j = std::floor(m * U(R));    // m needed...
	int k = (m-1) - j;
	while ((j>0) && (ly < lpdf(bounds[0])) ) { // trunc. can be added here.
		bounds[0] = bounds[0] - w;
		j = j - 1;
		if (!fp && (bounds[0] < *(peak_iter-1))) break;
	}
	while ((k>0) && (ly < lpdf(bounds[1])) ) { // truncation can be added here.
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

void RV_Missing_t_walk_observed_normal::find_peaks() {
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

	if (!arma::eig_gen(cx_eigval, cx_eigvec, companion)) 
		throw std::runtime_error("Failed eigenvalue decomposition.");
	arma::cx_vec::iterator re_eigval_end = 
		std::remove_if(cx_eigval.begin(), cx_eigval.end(), has_imaginary);
	std:sort(cx_eigval.begin(), re_eigval_end);
	num_real_eigval = re_eigval_end - cx_eigval.begin() - 1;
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
//	ecompanion(0,2) = companion(0,2);
//	ecompanion(1,2) = companion(1,2);
//	ecompanion(2,2) = companion(2,2);
//	Eigen::EigenSolver<Eigen::MatrixXd> es(ecompanion);
//	std::cout << "Peaks: \n" << es.eigenvalues() << std::endl;
	

}	

// Missing_t_walk_observed_interval after here.
RV_Missing_t_walk_observed_interval::RV_Missing_t_walk_observed_interval(
		double const & x1_,
		double 			 & X,
		double const & x3_,
		double const & os1_,
		double const & os2_,
		double const & p1_,   // degrees of freedom 1
		double const & p2_,   // degrees of freedom 2
		double const & s1_,   // scale 1
		double const & s2_,   // scale 2
		double const & Xmin_,
		double const & Xmax_,
		trng::yarn2 & R_
) :	x1(x1_), x2(X), x3(x3_),
		os1(os1_), os2(os2_),
		p1(p1_), p2(p2_), s1(s1_), s2(s2_), 
		Xmin(Xmin_), Xmax(Xmax_),
		companion(3,3), R(R_),
		bounds_pk1(2), bounds_pk2(2), EXPO(1.0)//,
//		ecompanion(3,3)
{
	companion.zeros();
	companion(1,0) = 1.0;
	companion(2,1) = 1.0;
// 	ecompanion(0,0) = 0.0; ecompanion(0,1) = 0.0;
//	ecompanion(1,0) = 1.0; ecompanion(1,1) = 0.0;
//	ecompanion(2,0) = 0.0; ecompanion(2,1) = 1.0;
}

std::map<std::string, double> RV_Missing_t_walk_observed_interval::state() const {
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
	out["Xmin"] = Xmin;
	out["Xmax"] = Xmax;
	return out;
}

void RV_Missing_t_walk_observed_interval::jump(double X) { x2 = X; }

double RV_Missing_t_walk_observed_interval::draw() {  // Slice sampler w/ two peaks.
  ly = lpdf() - EXPO(R);
	find_slice();
	double ii = 0;
	while(true) {
		//ii++; // if (ii > 10000) throw std::runtime_error("Too many steps.");
		ii++; if (ii > 10000) { std::cout << "Too many steps." << std::endl; }
		x_new = choose();	
		ly_new = lpdf(x_new);
		if (ly_new >= ly) // && (x_new > min) && (x_new < max))  truncation.
			break; 
		else 
			trim();
	}
	x2 = x_new;
	return x2; 
}

double RV_Missing_t_walk_observed_interval::lpdf(double X) {
	if ((X < Xmin) || (X > Xmax)) {		
		return -1.0 * std::numeric_limits<double>::infinity();
	}
	double lpdf;
	lpdf = 	boost::math::lgamma((p1+1.0)/2.0) - boost::math::lgamma(p1/2.0) -
		0.5 * log(p1*pi*pow(s1,2)) - 
		(p1+1.0)/2.0 * log(1.0 + (pow(X -(x1+os1),2))/(p1*pow(s1,2)) ) +
					boost::math::lgamma((p2+1.0)/2.0) - boost::math::lgamma(p2/2.0) -
		0.5 * log(p2*pi*pow(s2,2)) - 
		(p2+1.0)/2.0 * log(1.0 + (pow(x3-(X +os2),2))/(p2*pow(s2,2)) );
	return lpdf;
}

double RV_Missing_t_walk_observed_interval::lpdf() { return lpdf(x2); }

void RV_Missing_t_walk_observed_interval::find_slice() {
	find_peaks();
	if (ly <= lpdf(peak1)) bounds_pk1 = step_out(peak1);
	if (ly <= lpdf(peak2)) bounds_pk2 = step_out(peak2);
	if (ly > lpdf(peak1)) {   // In case peak 1 does not reach slice level.
		bounds_pk1[0] = bounds_pk2[0];
		bounds_pk1[1] = bounds_pk2[0];
	}
	if (ly > lpdf(peak2)) {		// In case peak 2 does not reach slice level.
		bounds_pk2[0] = bounds_pk1[0];
		bounds_pk2[1] = bounds_pk1[0];
	}
	if (  // In case the two parts of the slice overlap, split it:
		(bounds_pk1[0] < std::min(bounds_pk1[1], bounds_pk2[1])) &&
		(bounds_pk2[0] < std::min(bounds_pk1[1], bounds_pk2[1]))
	) {
		double x = (bounds_pk1[0]+bounds_pk2[1])/2.0;
		bounds_pk1[1] = x;
		bounds_pk2[0] = x;
	}
}


std::vector<double> RV_Missing_t_walk_observed_interval::step_out(double peak) {
	std::vector<double> bounds(2);
	int m = 200;
	double w = (s1+s2)/2.0;
	bounds[0] = peak - w * U(R);  // w = use (s1+s2)/2
	bounds[1] = bounds[0] + w;
	int j = std::floor(m * U(R));    // m needed...
	int k = (m-1) - j;
	while ((j>0) && (ly < lpdf(bounds[0])) ) { // trunc. can be added here.
		bounds[0] = bounds[0] - w;
		j = j - 1;
	}
	while ((k>0) && (ly < lpdf(bounds[1])) ) { // truncation can be added here.
		bounds[1] = bounds[1] + w;
		k = k - 1;
	}
	return bounds;
}

void RV_Missing_t_walk_observed_interval::trim() {
	if (x_new <= bounds_pk1[1]) {
		if (x_new <= peak1) {
			bounds_pk1[0] = x_new;
		} else {
			bounds_pk1[1] = x_new;
		}
	} else {
		if (x_new <= peak2) {
			bounds_pk2[0] = x_new;
		} else {
			bounds_pk2[1] = x_new;
		}
	}
}

double RV_Missing_t_walk_observed_interval::choose() {
	double d = (bounds_pk2[1] - bounds_pk2[0]) + (bounds_pk1[1] - bounds_pk1[0]);
	double q = (bounds_pk2[0] - bounds_pk1[1]);
	double l = U(R) * d;
	if ((bounds_pk1[0] + l) > bounds_pk1[1]) {
		return (l + q + bounds_pk1[0]);
	} else {
		return (l + bounds_pk1[0]);
	}
}

void RV_Missing_t_walk_observed_interval::find_peaks() {
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

  companion(0,2) = -B4/B1;
	companion(1,2) = -B3/B1;
	companion(2,2) = -B2/B1;
	
	if (!arma::eig_gen(cx_eigval, cx_eigvec, companion)) 
		throw std::runtime_error("Failed eigenvalue decomposition.");
	arma::cx_vec::iterator re_eigval_end = 
		std::remove_if(cx_eigval.begin(), cx_eigval.end(), has_imaginary);
	std:sort(cx_eigval.begin(), re_eigval_end);
	num_real_eigval = re_eigval_end - cx_eigval.begin() - 1;
	if (num_real_eigval == 0) {
		throw std::runtime_error("No real eigenvalues available to find peaks.");
	} else if (num_real_eigval == 1) {
		peak1 = arma::real(cx_eigval)[0];
		peak2 = arma::real(cx_eigval)[0];
	} else if (num_real_eigval == 3) {
		peak1 = arma::real(cx_eigval)[0];
		peak2 = arma::real(cx_eigval)[2];
	} else {
		throw std::runtime_error("Number of real eigenvalues not equal to 0, 1, or 3.");
	}
//	ecompanion(0,2) = companion(0,2);
//	ecompanion(1,2) = companion(1,2);
//	ecompanion(2,2) = companion(2,2);
//	Eigen::EigenSolver<Eigen::MatrixXd> es(ecompanion);
//	std::cout << "Peaks: \n" << es.eigenvalues() << std::endl;
	
	
}
