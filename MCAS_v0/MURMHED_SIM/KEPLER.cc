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

double MOID_bonanno(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<std::vector<double> >* m_vec, unsigned int id_1, unsigned int id_2) {
/*	std::vector<double> kep1;
	std::vector<double> kep2;
	kep1 = qp_to_kepler((*q_vec)[id_1], (*p_vec)[id_1], (*m_vec)[id_1], (*m_vec)[0]);
	kep2 = qp_to_kepler((*q_vec)[id_2], (*p_vec)[id_2], (*m_vec)[id_2], (*m_vec)[0]);
*/
	double Dmin;
	std::vector<double> T;
	T = cross_product((*p_vec)[id_1]/(*m_vec)[id_1],(*p_vec)[id_2]/(*m_vec)[id_2]);

	Dmin = (((*q_vec)[id_1][0])*((*q_vec)[id_1][0])*T[0]*T[0])/dot_product(T,T);

	return Dmin;
}

double Tisserand(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double> *m_vec, unsigned int id_1, unsigned int id_2) {
	std::vector<double> kep1;
	std::vector<double> kep2;
	kep1 = qp_to_kepler((*q_vec)[id_1], (*p_vec)[id_1], (*m_vec)[id_1], (*m_vec)[0]);
	kep2 = qp_to_kepler((*q_vec)[id_2], (*p_vec)[id_2], (*m_vec)[id_2], (*m_vec)[0]);

	double T = kep2[0]/kep1[0] + 2.0*sqrt(((kep1[0])/(kep2[0]))*(1 - kep1[1]*kep1[1]))*cos(kep1[2]);
	
	return T;
}

std::vector<double> total_angular_momentum_qv_heliocentric(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* v_vec, std::vector<double> *m_vec) {
	std::vector<double> ret;
	ret.push_back(0);
	ret.push_back(0);
	ret.push_back(0);

	int i;
	for(i=(*q_vec).size()-1; i >= 1; i--) {
		ret = ret + cross_product(((*q_vec)[i]),((*m_vec)[i])*((*v_vec)[i]));
	}

	return ret;
}

int rot_to_invariable_plane(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* v_vec, std::vector<double>* m_vec) {
	std::cout << "Converting coordinate system to the invariable plane." << std::endl;
	std::vector<double> L;
	L = total_angular_momentum_qv_heliocentric(q_vec, v_vec, m_vec);
	output_stream out(std::cout);
	out.type = 0;

	double L_t = abs_v(L);
	double phi = acos(L[2]/L_t);
	double theta = atan2(L[1],L[0]);
	std::cout << "Rotation around invertable plan of theta=" << theta*(180/PI) << ", phi=" << phi*(180/PI) << " detected." << std::endl;
	std::cout << "Total angular momentum: ";
	print_vector(L/L_t,out);
	std::cout << "Anti rotation result: ";
	print_vector(rot_to_plane(L/L_t,L),out);
	int i;
	for(i=(*q_vec).size()-1; i >= 1; i--) {
		((*q_vec)[i]) = rot_to_plane(((*q_vec)[i]),L);
		((*v_vec)[i]) = rot_to_plane(((*v_vec)[i]),L);
	}
	std::cout << "Checking sun coordinates: q " << abs_v(((*q_vec)[0])) << ", v " << abs_v(((*v_vec)[0])) << std::endl;
	std::cout << "Conversion done." << std::endl << std::endl;
	std::cout << "Checking new angular momentum:" << std::endl;
	L = total_angular_momentum_qv_heliocentric(q_vec, v_vec, m_vec);
	print_vector(L,out);

	return CALC_OK;
}

int rot_all_y(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* v_vec, double theta) {
	std::cout << "Rotating all objects around y axis by " << theta*(180/PI) << " degrees " << std::endl;
	int i;
	for(i=(*q_vec).size()-1; i >= 1; i--) {
		((*q_vec)[i]) = rot_y(((*q_vec)[i]),theta);
		((*v_vec)[i]) = rot_y(((*v_vec)[i]),theta);
	}
	std::cout << "Rotation complete" << std::endl;

	return CALC_OK;
}

int rot_all_z(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* v_vec, double theta) {
	std::cout << "Rotating all objects around z axis by " << theta*(180/PI) << " degrees " << std::endl;
	int i;
	for(i=(*q_vec).size()-1; i >= 1; i--) {
		((*q_vec)[i]) = rot_z(((*q_vec)[i]),theta);
		((*v_vec)[i]) = rot_z(((*v_vec)[i]),theta);
	}
	std::cout << "Rotation complete" << std::endl;

	return CALC_OK;
}

void heliocentric_qp(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vec) {
	int i;
	for(i=(*q_vec).size()-1; i >= 0; i--) {
		((*q_vec)[i])[0] = ((*q_vec)[i])[0] - ((*q_vec)[0])[0];
		((*q_vec)[i])[1] = ((*q_vec)[i])[1] - ((*q_vec)[0])[1];
		((*q_vec)[i])[2] = ((*q_vec)[i])[2] - ((*q_vec)[0])[2];

		((*p_vec)[i])[0] = ((*p_vec)[i])[0] - ((*p_vec)[0])[0]*(((*m_vec)[i])/((*m_vec)[0]));
		((*p_vec)[i])[1] = ((*p_vec)[i])[1] - ((*p_vec)[0])[1]*(((*m_vec)[i])/((*m_vec)[0]));
		((*p_vec)[i])[2] = ((*p_vec)[i])[2] - ((*p_vec)[0])[2]*(((*m_vec)[i])/((*m_vec)[0]));
	}
}

void heliocentric_qv(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* v_vec) {
	int i;
	for(i=(*q_vec).size()-1; i >= 0; i--) {
		((*q_vec)[i])[0] = ((*q_vec)[i])[0] - ((*q_vec)[0])[0];
		((*q_vec)[i])[1] = ((*q_vec)[i])[1] - ((*q_vec)[0])[1];
		((*q_vec)[i])[2] = ((*q_vec)[i])[2] - ((*q_vec)[0])[2];

		((*v_vec)[i])[0] = ((*v_vec)[i])[0] - ((*v_vec)[0])[0];
		((*v_vec)[i])[1] = ((*v_vec)[i])[1] - ((*v_vec)[0])[1];
		((*v_vec)[i])[2] = ((*v_vec)[i])[2] - ((*v_vec)[0])[2];
	}
}

std::vector<double> kepler_to_qp(std::vector<double> kepler, double m, double M) {
   	/*INPUT: a,e,i,omega,Omega,nu,central body mass*/
	unsigned int i;
	std::vector<double> return_vec;
  double r_vec_kep_state[3], v_vec_kep_state[3], mu_kep_state;
  double p_kep_state, ecc_kep_state, incl_kep_state, omega_kep_state, argp_kep_state, nu_kep_state, arglat_kep_state, truelon_kep_state, lonper_kep_state;

	mu_kep_state = G*(m+M)*1e-9;
    p_kep_state = kepler[0]*1e-3*(1-pow(kepler[1],2));
    ecc_kep_state = kepler[1];
        incl_kep_state = kepler[2];
        omega_kep_state = kepler[4];
        argp_kep_state = kepler[3];
        nu_kep_state = kepler[5];
        arglat_kep_state = nu_kep_state + argp_kep_state;
        lonper_kep_state = omega_kep_state + argp_kep_state;
        truelon_kep_state = nu_kep_state + lonper_kep_state;

    coe2rv(p_kep_state, ecc_kep_state, incl_kep_state, omega_kep_state, argp_kep_state, nu_kep_state, arglat_kep_state, truelon_kep_state, lonper_kep_state, r_vec_kep_state, v_vec_kep_state, mu_kep_state);

    for(i=0; i < 3; i++) {
    	return_vec.push_back(r_vec_kep_state[i]*1e3);//pos (m)
    }
	for(i=0; i < 3; i++) {
    	return_vec.push_back(v_vec_kep_state[i]*1e3*m); //momentum kg m /s
    }

	return return_vec;

   	/*//std::vector<double> P;
   	//std::vector<double> S;
	std::vector<double> o;
	std::vector<double> o_dot;
	std::vector<double> x_trans;
	std::vector<double> y_trans;
   	std::vector<double> return_vec;
   	//kepler 0=a, 1=e, 2=i, 3=omega, 4=Omega, 5=nu, 6=E
   	//double eta = kepler[0]*(1-kepler[1]*kepler[1]);
   	double r = kepler[0]*(1-kepler[1]*cos(kepler[6]));

   	o.push_back(r*cos(kepler[5]));
   	o.push_back(r*sin(kepler[5]));
   	o.push_back(0);

   	o_dot.push_back(((sqrt(G*(M)*kepler[0]))/r)*(-sin(kepler[6])));
   	o_dot.push_back(((sqrt(G*(M)*kepler[0]))/r)*(sqrt(1-kepler[1]*kepler[1])*cos(kepler[6])));
   	o_dot.push_back(0);

   	o = rot_z(rot_y(rot_z(o,kepler[3]),kepler[2]),kepler[4]+PI/2);
   	o_dot = rot_z(rot_y(rot_z(o_dot,kepler[3]),kepler[2]),kepler[4]+PI/2);
   	return_vec.insert(return_vec.end(),o.begin(),o.end());
   	return_vec.insert(return_vec.end(),o_dot.begin(),o_dot.end());

   	return_vec[3] = return_vec[3]*m;
   	return_vec[4] = return_vec[4]*m;
   	return_vec[5] = return_vec[5]*m;

    return return_vec;*/
}

double nu_to_E(double e,double nu) {
	return atan(sqrt(1-pow(e,2))*sin(nu)/(e + cos(nu)));
}

std::vector<double> kepler_to_xv(std::vector<double> kepler, double mu) {
	/*INPUT: a,e,i,omega,Omega,nu,central body mass*/
	unsigned int i;
	std::vector<double> return_vec;
  double r_vec_kep_state[3], v_vec_kep_state[3], mu_kep_state;
  double p_kep_state, ecc_kep_state, incl_kep_state, omega_kep_state, argp_kep_state, nu_kep_state, arglat_kep_state, truelon_kep_state, lonper_kep_state;

	mu_kep_state = mu*1e-9; //km^3/s^2
    p_kep_state = kepler[0]*1e-3*(1-pow(kepler[1],2)); //km
    ecc_kep_state = kepler[1];
        incl_kep_state = kepler[2]; //rad
        omega_kep_state = kepler[4];//rad
        argp_kep_state = kepler[3];//rad
        nu_kep_state = kepler[5];//rad
        arglat_kep_state = nu_kep_state + argp_kep_state;
        lonper_kep_state = omega_kep_state + argp_kep_state;
        truelon_kep_state = nu_kep_state + lonper_kep_state;

    coe2rv(p_kep_state, ecc_kep_state, incl_kep_state, omega_kep_state, argp_kep_state, nu_kep_state, arglat_kep_state, truelon_kep_state, lonper_kep_state, r_vec_kep_state, v_vec_kep_state, mu_kep_state);

    for(i=0; i < 3; i++) {
    	return_vec.push_back(r_vec_kep_state[i]*1e3);//m
    }
	for(i=0; i < 3; i++) {
    	return_vec.push_back(v_vec_kep_state[i]*1e3); //m/s
    }

	return return_vec;
   	/*//std::vector<double> P;
   	//std::vector<double> S;
	std::vector<double> o;
	std::vector<double> o_dot;
	std::vector<double> x_trans;
	std::vector<double> y_trans;
   	std::vector<double> return_vec;
   	//kepler 0=a, 1=e, 2=i, 3=omega, 4=Omega, 5=nu, 6=E
   	//double eta = kepler[0]*(1-kepler[1]*kepler[1]);
   	double r = kepler[0]*(1-kepler[1]*cos(kepler[6]));

   	o.push_back(r*cos(kepler[5]));
   	o.push_back(r*sin(kepler[5]));
   	o.push_back(0);

   	o_dot.push_back(((sqrt(G*(M)*kepler[0]))/r)*(-sin(kepler[6])));
   	o_dot.push_back(((sqrt(G*(M)*kepler[0]))/r)*(sqrt(1-kepler[1]*kepler[1])*cos(kepler[6])));
   	o_dot.push_back(0);

   	o = rot_z(rot_y(rot_z(o,kepler[3]),kepler[2]),kepler[4]+PI/2);
   	o_dot = rot_z(rot_y(rot_z(o_dot,kepler[3]),kepler[2]),kepler[4]+PI/2);
   	return_vec.insert(return_vec.end(),o.begin(),o.end());
   	return_vec.insert(return_vec.end(),o_dot.begin(),o_dot.end());
 */
    
}


std::vector<double> xv_to_kepler(std::vector<double> x, std::vector<double> v, double mu) {
	
	/*INPUT: x_vec,v_vec,central body mass*/
	std::vector<double> return_vec;
  double r_vec_kep_state[3], v_vec_kep_state[3], mu_kep_state;
  double p_kep_state, a_kep_state, ecc_kep_state, incl_kep_state, omega_kep_state, argp_kep_state, nu_kep_state, m_kep_state, arglat_kep_state, truelon_kep_state, lonper_kep_state;

r_vec_kep_state[0] = x[0]*1e-3; //km
r_vec_kep_state[1] = x[1]*1e-3;
r_vec_kep_state[2] = x[2]*1e-3;

v_vec_kep_state[0] = v[0]*1e-3;
v_vec_kep_state[1] = v[1]*1e-3;
v_vec_kep_state[2] = v[2]*1e-3; //km/s

	mu_kep_state = mu*1e-9;

    rv2coe(r_vec_kep_state, v_vec_kep_state, mu_kep_state, p_kep_state, a_kep_state, ecc_kep_state, incl_kep_state, omega_kep_state, argp_kep_state, nu_kep_state, m_kep_state, arglat_kep_state, truelon_kep_state, lonper_kep_state);

    return_vec.push_back(a_kep_state*1e3); //m
    return_vec.push_back(ecc_kep_state);
    return_vec.push_back(incl_kep_state); //rad
    return_vec.push_back(argp_kep_state);
    return_vec.push_back(omega_kep_state);
    return_vec.push_back(nu_kep_state);
    return_vec.push_back(m_kep_state);
    return_vec.push_back(arglat_kep_state);
    return_vec.push_back(truelon_kep_state);
    return_vec.push_back(lonper_kep_state);
    return_vec.push_back(p_kep_state*1e3);

	return return_vec;

	/*std::vector<double> e;
	std::vector<double> e_temp;
	std::vector<double> h;
	std::vector<double> n;
	std::vector<double> v_v;
	std::vector<double> r_v;

	std::vector<double> return_vec;

	double mu = G*(M);

	std::vector<double> k;
	k.push_back(0);k.push_back(0);k.push_back(1);

	v_v.push_back(v[0]);
	v_v.push_back(v[1]);
	v_v.push_back(v[2]);
	r_v.push_back(x[0]);
	r_v.push_back(x[1]);
	r_v.push_back(x[2]);


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

	return return_vec;*/

}

std::vector<double> qp_to_kepler(std::vector<double> q, std::vector<double> p, double m, double M) {
		/*INPUT: x_vec,v_vec,central body mass*/
	std::vector<double> return_vec;
  double r_vec_kep_state[3], v_vec_kep_state[3], mu_kep_state;
  double p_kep_state, a_kep_state, ecc_kep_state, incl_kep_state, omega_kep_state, argp_kep_state, nu_kep_state, m_kep_state, arglat_kep_state, truelon_kep_state, lonper_kep_state;

r_vec_kep_state[0] = q[0]*1e-3;
r_vec_kep_state[1] = q[1]*1e-3;
r_vec_kep_state[2] = q[2]*1e-3;

v_vec_kep_state[0] = p[0]*1e-3/m;
v_vec_kep_state[1] = p[1]*1e-3/m;
v_vec_kep_state[2] = p[2]*1e-3/m;

	mu_kep_state = G*(m+M)*1e-9;

    rv2coe(r_vec_kep_state, v_vec_kep_state, mu_kep_state, p_kep_state, a_kep_state, ecc_kep_state, incl_kep_state, omega_kep_state, argp_kep_state, nu_kep_state, m_kep_state, arglat_kep_state, truelon_kep_state, lonper_kep_state);

    return_vec.push_back(a_kep_state*1e3);
    return_vec.push_back(ecc_kep_state);
    return_vec.push_back(incl_kep_state);
    return_vec.push_back(argp_kep_state);
    return_vec.push_back(omega_kep_state);
    return_vec.push_back(nu_kep_state);
    return_vec.push_back(m_kep_state);
    return_vec.push_back(arglat_kep_state);
    return_vec.push_back(truelon_kep_state);
    return_vec.push_back(lonper_kep_state);
    return_vec.push_back(p_kep_state*1e3);

	return return_vec;

	/*std::vector<double> e;
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

	return return_vec;*/

}

/*
std::vector<std::vector<double> > H_kep_qp(std::vector<double> q_vec, std::vector<double> p_vec, double m, double M, double dt) {
	unsigned int j;
	double mu, a, r_0, v_r0, Xi, z, r, alpha;
	double f, g, dotf, dotg, S, C;

	std::vector<double> temp_q, temp_v;
	std::vector<std::vector<double> > ret;
	temp_q.push_back(0);temp_q.push_back(0);temp_q.push_back(0);
	temp_v.push_back(0);temp_v.push_back(0);temp_v.push_back(0);

			mu = G*(m + M);
			r_0 = abs_v(q_vec);
			a = 1/(2/r_0 - abs_v(p_vec)/(mu*pow(m,2)));
			alpha = 1/a;
			v_r0 = dot_product(p_vec, q_vec)/(r_0*m);

			Xi = universal_kepler_equation_newton(dt, 10e-8, a, r_0, v_r0, mu);
			//std::cout << a/AU << " " << v_r0 << std::endl;

			z = alpha*pow(Xi,2);
			S = stumpff(3,6,z);
			C = stumpff(2,6,z);

			f = 1 - pow(Xi,2)/r_0*C;
			g = dt - 1/sqrt(mu)*pow(Xi,3)*S;

			temp_q[0] = f*(q_vec[0]) + g*(p_vec[0]/m);
			temp_q[1] = f*(q_vec[1]) + g*(p_vec[1]/m);
			temp_q[2] = f*(q_vec[2]) + g*(p_vec[2]/m);

			r = abs_v(temp_q);

			dotf = sqrt(mu)/(r*r_0)*(alpha*pow(Xi,3)*S - Xi);
			dotg = 1 - pow(Xi,2)/r*C;

			temp_v[0] = dotf*(q_vec[0]) + dotg*(p_vec[0]/m);
			temp_v[1] = dotf*(q_vec[1]) + dotg*(p_vec[1]/m);
			temp_v[2] = dotf*(q_vec[2]) + dotg*(p_vec[2]/m);

			ret.push_back(temp_q);
			ret.push_back(temp_v);

			return ret;

}
*/
std::vector<std::vector<double> > H_kep_qp(std::vector<double> q_vec, std::vector<double> p_vec, double m, double M, double dt) {

	std::vector<double> temp_vec;
	std::vector<std::vector<double> > ret;

	double ro[3], vo[3], r[3], v[3];

	ro[0] = q_vec[0]*1e-3;
	ro[1] = q_vec[1]*1e-3;
	ro[2] = q_vec[2]*1e-3;

	vo[0] = p_vec[0]*1e-3/m;
	vo[1] = p_vec[1]*1e-3/m;
	vo[2] = p_vec[2]*1e-3/m;

	kepler(ro, vo, dt, G*(M+m)*1e-9, r, v);

	temp_vec.push_back(r[0]*1e3);
	temp_vec.push_back(r[1]*1e3);
	temp_vec.push_back(r[2]*1e3);
	ret.push_back(temp_vec);

	temp_vec[0] = v[0]*m*1e3;
	temp_vec[1] = v[1]*m*1e3;
	temp_vec[2] = v[2]*m*1e3;
	ret.push_back(temp_vec);

			return ret;

}

std::vector<std::vector<double> > H_kep_xv(std::vector<double> q_vec, std::vector<double> v_vec, double M, double dt) {

	std::vector<double> temp_vec;
	std::vector<std::vector<double> > ret;

	double ro[3], vo[3], r[3], v[3];

	ro[0] = q_vec[0]*1e-3;
	ro[1] = q_vec[1]*1e-3;
	ro[2] = q_vec[2]*1e-3;

	vo[0] = v_vec[0]*1e-3;
	vo[1] = v_vec[1]*1e-3;
	vo[2] = v_vec[2]*1e-3;

	kepler(ro, vo, dt, G*M*1e-9, r, v);

	temp_vec.push_back(r[0]*1e3);
	temp_vec.push_back(r[1]*1e3);
	temp_vec.push_back(r[2]*1e3);
	ret.push_back(temp_vec);

	temp_vec[0] = v[0]*1e3;
	temp_vec[1] = v[1]*1e3;
	temp_vec[2] = v[2]*1e3;
	ret.push_back(temp_vec);

			return ret;
}
/*
std::vector<std::vector<double> > H_kep_xv(std::vector<double> q_vec, std::vector<double> v_vec, double M, double dt) {
	unsigned int j;
	double mu, a, r_0, v_r0, Xi, z, r, alpha;
	double f, g, dotf, dotg, S, C;

	std::vector<double> temp_q, temp_v;
	std::vector<std::vector<double> > ret;
	temp_q.push_back(0);temp_q.push_back(0);temp_q.push_back(0);
	temp_v.push_back(0);temp_v.push_back(0);temp_v.push_back(0);

			mu = G*(M);
			r_0 = abs_v(q_vec);
			a = 1/(2/r_0 - abs_v(v_vec)/mu);
			alpha = 1/a;
			v_r0 = dot_product(v_vec, q_vec)/r_0;

			Xi = universal_kepler_equation_newton(dt, 10e-8, a, r_0, v_r0, mu);
			//std::cout << a/AU << " " << v_r0 << std::endl;

			z = alpha*pow(Xi,2);
			S = stumpff(3,6,z);
			C = stumpff(2,6,z);

			f = 1 - pow(Xi,2)/r_0*C;
			g = dt - 1/sqrt(mu)*pow(Xi,3)*S;

			temp_q[0] = f*(q_vec[0]) + g*v_vec[0];
			temp_q[1] = f*(q_vec[1]) + g*v_vec[1];
			temp_q[2] = f*(q_vec[2]) + g*v_vec[2];

			r = abs_v(temp_q);

			dotf = sqrt(mu)/(r*r_0)*(alpha*pow(Xi,3)*S - Xi);
			dotg = 1 - pow(Xi,2)/r*C;

			temp_v[0] = dotf*(q_vec[0]) + dotg*v_vec[0];
			temp_v[1] = dotf*(q_vec[1]) + dotg*v_vec[1];
			temp_v[2] = dotf*(q_vec[2]) + dotg*v_vec[2];

			ret.push_back(temp_q);
			ret.push_back(temp_v);

			return ret;
}*/

double universal_kepler_equation_newton(double dt, double tol, double a, double r_0, double v_r0, double mu) {

	double Xi, Xi_0, dFdXi, F, dXi, S, C, z, alpha, abs_alpha;
	if(alpha < 0) {
		abs_alpha = -alpha;
	}
	else {
		abs_alpha = alpha;
	}

	Xi_0 = sqrt(mu)*abs_alpha*dt;
	alpha = 1/a;
	Xi = Xi_0;

	do {
		z = alpha*pow(Xi,2);
		S = stumpff(3,6,z);
		C = stumpff(2,6,z);

		F = r_0*v_r0/sqrt(mu)*pow(Xi,2)*C + (1- alpha*r_0)*pow(Xi,3)*S + r_0*Xi - sqrt(mu)*dt;
		dFdXi = r_0*v_r0/sqrt(mu)*Xi*(1 - alpha*pow(Xi,2)*S) + (1 - alpha*r_0)*pow(Xi,2)*C + r_0;
		dXi = F/dFdXi;

		Xi = Xi - dXi;
	} while(dXi > tol);
	//std::cout << "Xi0 " << Xi_0 << ", Xi " << Xi << std::endl;
	return Xi;
}


double kepler_equation_newton(double M, double e, double E_0, double tol) {

	double E_temp = E_0;
	double F_temp = E_temp - e*sin(E_temp) - M;

	while(F_temp*F_temp > tol*tol) {
		E_temp = E_temp - (F_temp)/(1-e*cos(E_temp));
		F_temp = E_temp - e*sin(E_temp) - M;
	}
	//E_temp = fmod(E_temp,2*PI);
	return E_temp;
}
