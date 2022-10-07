/*
 * kepler.cc
 *
 *  Created on: Oct 6, 2014
 *      Author: dankas
 */
#include <math.h>
#include <vector>
#include <iostream>

#include "ast2body.h"
#include "define.hh"
#include "functions.hh"


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
