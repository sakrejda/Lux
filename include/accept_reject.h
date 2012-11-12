#ifndef AR_H
#define AR_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

#include <boost/spirit/home/support/detail/hold_any.hpp>

#include <flipper.h>


template <class P, class C> class AcceptReject {
	P a;
	P b;
	C proposal;
	FlipFlop<P> ff;
	double lpd_current;
	double lpd_propose;
	double asymmetry;
	void calculate_acceptance();
	double acceptance;
	bool a_stale;
	bool accept_move;

public:
	AcceptReject(SEXP x);
	AcceptReject(P init_parameters, C init_proposal);
	AcceptReject();
	void propose();
	double get_a();
	std::map<std::string, boost::spirit::hold_any> accept();
	std::map<std::string, boost::spirit::hold_any> reject();
	std::map<std::string, boost::spirit::hold_any> report();
	std::map<std::string, boost::spirit::hold_any> step(double r);
};



template <class P, class C> 
AcceptReject<P, C>::AcceptReject (SEXP x) : 
	a(x), b(x), proposal(x),
	asymmetry(0.0), acceptance(0.0),
	a_stale(true), accept_move(true)
{
	ff = FlipFlop<P>(&a, &b);
}

template <class P, class C> 
AcceptReject<P, C>::AcceptReject (P init_parameters, C init_proposal) {
	a = init_parameters;
	b = init_parameters;
	proposal = init_proposal;
	ff = FlipFlop<P>(&a, &b);
	asymmetry = 0.0;
	acceptance = 0.0;
	a_stale = true;
	accept_move = true;
}

template <class P, class C> 
AcceptReject<P, C>::AcceptReject () {
	P init_parameters;
	C init_proposal;
	P a = init_parameters;
	P b = init_parameters;
	proposal = init_proposal;
	ff = FlipFlop<P>(&a, &b);
	asymmetry = 0.0;
	acceptance = 0.0;
	a_stale = true;
	accept_move = true;
}


template <class P, class C> 
void AcceptReject<P, C>::calculate_acceptance() {
	lpd_current = ff.focus_pointer()->get_lpd();
	lpd_propose = ff.other_pointer()->get_lpd();
  acceptance = exp(lpd_propose - lpd_current + asymmetry);
}

template <class P, class C> 
double AcceptReject<P, C>::get_a() {
	if (a_stale) {
		calculate_acceptance();
		a_stale = false;
	}
	return acceptance;
}

template <class P, class C> 
void AcceptReject<P, C>::propose() {
	asymmetry = proposal.propose(ff.other_pointer(), ff.focus_pointer());
	a_stale = true;
}

//Stack overflow example, due to Jere.Jones:
//map<int, int> m;
//vector<int> v;
//for(map<int,int>::iterator it = m.begin(); it != m.end(); ++it) {
//  v.push_back(it->first);    // makes a vector of key values...
//  cout << it->first << "\n";
//}

template <class P, class C> 
std::map<std::string, boost::spirit::hold_any> AcceptReject<P, C>::report() {   // init parameter is now meaningless since we don't have to cough up keys separately---sweet!!
	std::map<std::string, boost::spirit::hold_any> parameter_values;
	std::map<std::string, boost::spirit::hold_any> report_map;

//		msg_stream << get_a()  << ", " << ff.focus_pointer()->get_lpd() << ", " << asymmetry  << ", " << accept_move << ", " << ff.focus_pointer()->report(init) << std::endl;
	parameter_values["acceptance"] = get_a();
	parameter_values["lpd"] = ff.focus_pointer()->get_lpd();
	parameter_values["asymmetry"] = asymmetry;
	parameter_values["accept"] = accept_move;
	report_map = ff.focus_pointer()->report();
	parameter_values.insert(report_map.begin(), report_map.end());

	return parameter_values;
}

template <class P, class C> 
std::map<std::string, boost::spirit::hold_any> AcceptReject<P, C>::accept() {
	accept_move = true;
	ff.flip();
	std::map<std::string, boost::spirit::hold_any> parameter_values = report();
	return parameter_values;
}

template <class P, class C> 
std::map<std::string, boost::spirit::hold_any> AcceptReject<P, C>::reject() {
	accept_move = false;
	std::map<std::string, boost::spirit::hold_any> parameter_values = report();
	return parameter_values;
}

template <class P, class C>
std::map<std::string, boost::spirit::hold_any> AcceptReject<P, C>::step(double r) {
	std::map<std::string, boost::spirit::hold_any> parameter_values;
	propose();
	if (r < get_a()) {
		parameter_values = accept();
	} else {
		parameter_values = reject();
	}
	return parameter_values;
}

#endif
