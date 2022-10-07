/*
 * kepler.cc
 *
 *  Created on: Oct 6, 2014
 *      Author: dankas
 */
#include <math.h>
#include <vector>
#include <iostream>

#include "define.hh"
#include "functions.hh"

std::vector<double> kepler_to_qp(std::vector<double> kepler, double m, double M) {
   	//std::vector<double> P;
   	//std::vector<double> S;
	std::vector<double> o;
	std::vector<double> o_dot;
	std::vector<double> x_trans;
	std::vector<double> y_trans;
   	std::vector<double> return_vec;
   	//kepler 0=a, 1=e, 2=i, 3=omega, 4=Omega, 5=nu, 6=E
	double E = kepler[6];
   	//double eta = kepler[0]*(1-kepler[1]*kepler[1]);
   	double r = kepler[0]*(1-kepler[1]*cos(E));

   	o.push_back(r*cos(kepler[5]));
   	o.push_back(r*sin(kepler[5]));

   	o_dot.push_back(((sqrt(G*(M+m)*kepler[0]))/r)*(-sin(E)));
   	o_dot.push_back(((sqrt(G*(M+m)*kepler[0]))/r)*(sqrt(1-kepler[1]*kepler[1])*cos(E)));

   	x_trans.push_back(cos(kepler[3])*cos(kepler[4]+PI/2) - sin(kepler[3])*cos(kepler[2])*sin(kepler[4]+PI/2));
   	x_trans.push_back(-cos(kepler[3])*sin(kepler[4]+PI/2) - sin(kepler[3])*cos(kepler[2])*cos(kepler[4]+PI/2));
   	x_trans.push_back(-sin(kepler[3])*sin(kepler[2]));

   	y_trans.push_back(sin(kepler[3])*cos(kepler[4]+PI/2) + cos(kepler[3])*cos(kepler[2])*sin(kepler[4]+PI/2));
   	y_trans.push_back(-sin(kepler[3])*sin(kepler[4]+PI/2) + cos(kepler[3])*cos(kepler[2])*cos(kepler[4]+PI/2));
   	y_trans.push_back(cos(kepler[3])*sin(kepler[2]));

   	return_vec.push_back(o[0]*x_trans[0] + o[1]*y_trans[0]);
   	return_vec.push_back(o[0]*x_trans[1] + o[1]*y_trans[1]);
   	return_vec.push_back(o[0]*x_trans[2] + o[1]*y_trans[2]);

   	return_vec.push_back(o_dot[0]*x_trans[0] + o_dot[1]*y_trans[0]);
   	return_vec.push_back(o_dot[0]*x_trans[1] + o_dot[1]*y_trans[1]);
   	return_vec.push_back(o_dot[0]*x_trans[2] + o_dot[1]*y_trans[2]);

   	return_vec[3] = return_vec[3]*m;
   	return_vec[4] = return_vec[4]*m;
   	return_vec[5] = return_vec[5]*m;

    return return_vec;
}

std::vector<double> qp_to_kepler(std::vector<double> q, std::vector<double> p, double m, double M) {
	std::vector<double> e;
	std::vector<double> e_temp;
	std::vector<double> h;
	std::vector<double> n;
	std::vector<double> v_v;
	std::vector<double> r_v;

	std::vector<double> return_vec;

	double mu = G*(M+m);

	std::vector<double> k;
	k.push_back(0);k.push_back(0);k.push_back(1);

	v_v.push_back(p[0]/m);
	v_v.push_back(p[1]/m);
	v_v.push_back(p[2]/m);
	r_v.push_back(q[0]);
	r_v.push_back(q[1]);
	r_v.push_back(q[2]);


	h = cross_product(r_v,v_v);
	e_temp = cross_product(v_v,h);
	e.push_back(e_temp[0]/mu - r_v[0]/abs_v(r_v));
	e.push_back(e_temp[1]/mu - r_v[1]/abs_v(r_v));
	e.push_back(e_temp[2]/mu - r_v[2]/abs_v(r_v));
	n = cross_product(k,h);

	double epsilon = dot_product(v_v,v_v)/2 - mu/abs_v(r_v);

	if(abs_v(e) <= 1) { //a
		return_vec.push_back(-mu/(2*epsilon));
	}
	else {
		return_vec.push_back(mu/(2*epsilon));
	}
	return_vec.push_back(abs_v(e)); //e

	return_vec.push_back(acos( h[2]/abs_v(h) )); //i

	if(return_vec[1] < 1e-10) { //omega
		return_vec.push_back(0);
	}
	else if(return_vec[2] > 1e-10) {
		if(e[2] >= 0) {
			return_vec.push_back(acos( dot_product(n,e)/(abs_v(n)*abs_v(e)) ));
		}
		else {
			return_vec.push_back(2*PI - acos( dot_product(n,e)/(abs_v(n)*abs_v(e)) ));
		}
	}
	else {
		if(h[2] >= 0) {
			return_vec.push_back(atan2(e[1],e[0]));
		}
		else {
			return_vec.push_back(2*PI - atan2(e[1],e[0]));
		}
	}

	if(return_vec[1] < 1e-10) { //Omega
			return_vec.push_back(0);
	}
	else if(n[1] >= 0 && return_vec[2] > 1e-10) {
		return_vec.push_back(acos( n[0]/abs_v(n) ));
	}
	else if(n[1] < 0 && return_vec[2] > 1e-10) {
		return_vec.push_back(2*PI - acos( n[0]/abs_v(n) ));
	}
	else {
		return_vec.push_back(0);
	}

	double temp_arg; //Nu
	if(return_vec[1] > 1e-10) {
		temp_arg = dot_product(e,r_v)/(return_vec[1]*abs_v(r_v));
		if(dot_product(r_v,v_v) >= 0) {
			return_vec.push_back(acos(temp_arg));
		}
		else {
			return_vec.push_back(2*PI - acos(temp_arg));
		}

	}
	else if(return_vec[2] > 1e-10) {
		temp_arg = dot_product(n,r_v)/(abs_v(n)*abs_v(r_v));
		if(dot_product(n,r_v) >= 0) {
			return_vec.push_back(acos(temp_arg));
		}
		else {
			return_vec.push_back(2*PI - acos(temp_arg));
		}
	}
	else {
		temp_arg = r_v[0]/abs_v(r_v);
		if(v_v[0] <= 0) {
			return_vec.push_back(acos(temp_arg));
		}
		else {
			return_vec.push_back(2*PI - acos(temp_arg));
		}
	}

	//E
	if(return_vec[5] < PI) {
		return_vec.push_back(acos((return_vec[1] + cos(return_vec[5]))/(1+return_vec[1]*cos(return_vec[5]))));
	}
	else {
		return_vec.push_back(2*PI - acos((return_vec[1] + cos(return_vec[5]))/(1+return_vec[1]*cos(return_vec[5]))));
	}

	return return_vec;

}


std::vector<double> xv_to_kepler(std::vector<double> q, std::vector<double> p, double M) {
	std::vector<double> e;
	std::vector<double> e_temp;
	std::vector<double> h;
	std::vector<double> n;
	std::vector<double> v_v;
	std::vector<double> r_v;

	std::vector<double> return_vec;

	double mu = G*(M);

	std::vector<double> k;
	k.push_back(0);k.push_back(0);k.push_back(1);

	v_v.push_back(p[0]);
	v_v.push_back(p[1]);
	v_v.push_back(p[2]);
	r_v.push_back(q[0]);
	r_v.push_back(q[1]);
	r_v.push_back(q[2]);


	h = cross_product(r_v,v_v);
	e_temp = cross_product(v_v,h);
	e.push_back(e_temp[0]/mu - r_v[0]/abs_v(r_v));
	e.push_back(e_temp[1]/mu - r_v[1]/abs_v(r_v));
	e.push_back(e_temp[2]/mu - r_v[2]/abs_v(r_v));
	n = cross_product(k,h);

	double epsilon = dot_product(v_v,v_v)/2 - mu/abs_v(r_v);

	if(abs_v(e) <= 1) { //a
		return_vec.push_back(-mu/(2*epsilon));
	}
	else {
		return_vec.push_back(mu/(2*epsilon));
	}
	return_vec.push_back(abs_v(e)); //e

	return_vec.push_back(acos( h[2]/abs_v(h) )); //i

	if(return_vec[1] < 1e-10) { //omega
		return_vec.push_back(0);
	}
	else if(return_vec[2] > 1e-10) {
		if(e[2] >= 0) {
			return_vec.push_back(acos( dot_product(n,e)/(abs_v(n)*abs_v(e)) ));
		}
		else {
			return_vec.push_back(2*PI - acos( dot_product(n,e)/(abs_v(n)*abs_v(e)) ));
		}
	}
	else {
		if(h[2] >= 0) {
			return_vec.push_back(atan2(e[1],e[0]));
		}
		else {
			return_vec.push_back(2*PI - atan2(e[1],e[0]));
		}
	}

	if(return_vec[1] < 1e-10) { //Omega
			return_vec.push_back(0);
	}
	else if(n[1] >= 0 && return_vec[2] > 1e-10) {
		return_vec.push_back(acos( n[0]/abs_v(n) ));
	}
	else if(n[1] < 0 && return_vec[2] > 1e-10) {
		return_vec.push_back(2*PI - acos( n[0]/abs_v(n) ));
	}
	else {
		return_vec.push_back(0);
	}

	double temp_arg; //Nu
	if(return_vec[1] > 1e-10) {
		temp_arg = dot_product(e,r_v)/(return_vec[1]*abs_v(r_v));
		if(dot_product(r_v,v_v) >= 0) {
			return_vec.push_back(acos(temp_arg));
		}
		else {
			return_vec.push_back(2*PI - acos(temp_arg));
		}

	}
	else if(return_vec[2] > 1e-10) {
		temp_arg = dot_product(n,r_v)/(abs_v(n)*abs_v(r_v));
		if(dot_product(n,r_v) >= 0) {
			return_vec.push_back(acos(temp_arg));
		}
		else {
			return_vec.push_back(2*PI - acos(temp_arg));
		}
	}
	else {
		temp_arg = r_v[0]/abs_v(r_v);
		if(v_v[0] <= 0) {
			return_vec.push_back(acos(temp_arg));
		}
		else {
			return_vec.push_back(2*PI - acos(temp_arg));
		}
	}

	//E
	if(return_vec[5] < PI) {
		return_vec.push_back(acos((return_vec[1] + cos(return_vec[5]))/(1+return_vec[1]*cos(return_vec[5]))));
	}
	else {
		return_vec.push_back(2*PI - acos((return_vec[1] + cos(return_vec[5]))/(1+return_vec[1]*cos(return_vec[5]))));
	}

	return return_vec;

}

std::vector<double> qp_to_kepler_sun(std::vector<double> qs, std::vector<double> ps, std::vector<double> q, std::vector<double> p, double m, double M) {
	std::vector<double> q_temp;
	std::vector<double> p_temp;
	std::vector<double> return_vec;

	q_temp.push_back(q[0] - qs[0]);
	q_temp.push_back(q[1] - qs[1]);
	q_temp.push_back(q[2] - qs[2]);

	p_temp.push_back(p[0] - m/M*ps[0]);
	p_temp.push_back(p[1] - m/M*ps[1]);
	p_temp.push_back(p[2] - m/M*ps[2]);

	return_vec = qp_to_kepler(q_temp,p_temp,m,M);
	return return_vec;
}
