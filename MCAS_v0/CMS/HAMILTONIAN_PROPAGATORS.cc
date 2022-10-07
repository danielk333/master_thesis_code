/*
 * H_symba.cc
 *
 *  Created on: Oct 6, 2014
 *      Author: dankas
 */
#include <iostream>
#include <vector>
#include <math.h>

#include "define.hh"
#include "functions.hh"

void H_pot(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec,double dt, unsigned int include_center) {
	unsigned int i,j;

	std::vector<std::vector<double> > temp_kick_x,temp_kick_y,temp_kick_z;
	double temp_kick_xd = 0,temp_kick_yd = 0,temp_kick_zd = 0;
	double temp_dist;
	std::vector<double> temp_p;

	for(i=0; i < (*q_vec).size()-1; i++) {
		if(((*type_vec)[j] != -1)) {
			temp_kick_x.push_back(temp_p);
			temp_kick_y.push_back(temp_p);
			temp_kick_z.push_back(temp_p);

			for(j=i+1; j < (*q_vec).size(); j++) {
				if(((*type_vec)[j] == 1)) {
					temp_dist = sqrt(( ((*q_vec)[i])[0] - ((*q_vec)[j])[0] )*( ((*q_vec)[i])[0] - ((*q_vec)[j])[0] ) + ( ((*q_vec)[i])[1] - ((*q_vec)[j])[1] )*( ((*q_vec)[i])[1] - ((*q_vec)[j])[1] ) + ( ((*q_vec)[i])[2] - ((*q_vec)[j])[2] )*( ((*q_vec)[i])[2] - ((*q_vec)[j])[2] ));
					temp_dist = pow(temp_dist,3);
					temp_kick_x[i].push_back( G*((*m_vector)[i])*((*m_vector)[j])*( ((*q_vec)[i])[0] - ((*q_vec)[j])[0] )/temp_dist);
					temp_kick_y[i].push_back( G*((*m_vector)[i])*((*m_vector)[j])*( ((*q_vec)[i])[1] - ((*q_vec)[j])[1] )/temp_dist);
					temp_kick_z[i].push_back( G*((*m_vector)[i])*((*m_vector)[j])*( ((*q_vec)[i])[2] - ((*q_vec)[j])[2] )/temp_dist);
				}
				else {
					temp_kick_x[i].push_back(0);
					temp_kick_y[i].push_back(0);
					temp_kick_z[i].push_back(0);
				}
			}
		}
	}
	for(i=(1-include_center); i < (*q_vec).size(); i++) {
		temp_kick_xd = 0; temp_kick_yd = 0; temp_kick_zd = 0;
		for(j=0; j < (*q_vec).size(); j++) {

			if(j < i) {
				//std::cout << "N_ " << temp_kick_x.size() << "M_ " << "N_ " << temp_kick_x[j].size() << " | " << i << " | " << j << std::endl;
				temp_kick_xd = temp_kick_xd + temp_kick_x[j][i-j-1];
				temp_kick_yd = temp_kick_yd + temp_kick_y[j][i-j-1];
				temp_kick_zd = temp_kick_zd + temp_kick_z[j][i-j-1];
			}
			else if(j > i) {
				//std::cout << "N_ " << temp_kick_x.size() << "M_ " << "N_ " << temp_kick_x[i].size() << " | " << i << " | " << j << std::endl;
				temp_kick_xd = temp_kick_xd + temp_kick_x[i][j-i-1];
				temp_kick_yd = temp_kick_yd + temp_kick_y[i][j-i-1];
				temp_kick_zd = temp_kick_zd + temp_kick_z[i][j-i-1];
			}
		}

		((*p_vec)[i])[0] = ((*p_vec)[i])[0] + dt*temp_kick_xd;
		((*p_vec)[i])[1] = ((*p_vec)[i])[1] + dt*temp_kick_yd;
		((*p_vec)[i])[2] = ((*p_vec)[i])[2] + dt*temp_kick_zd;
	}

}

void H_kin(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec, double dt) {
	unsigned int i;

	for(i=0; i < (*q_vec).size(); i++) {
		if((*type_vec)[i] != -1 && (*type_vec)[i] == 1) {
			((*q_vec)[i])[0] = ((*q_vec)[i])[0] + dt*((*p_vec)[i])[0]/((*m_vector)[i]);
			((*q_vec)[i])[1] = ((*q_vec)[i])[1] + dt*((*p_vec)[i])[1]/((*m_vector)[i]);
			((*q_vec)[i])[2] = ((*q_vec)[i])[2] + dt*((*p_vec)[i])[2]/((*m_vector)[i]);
		}
		else if((*type_vec)[i] != -1 && (*type_vec)[i] == 0) {
			((*q_vec)[i])[0] = ((*q_vec)[i])[0] + dt*((*p_vec)[i])[0];
			((*q_vec)[i])[1] = ((*q_vec)[i])[1] + dt*((*p_vec)[i])[1];
			((*q_vec)[i])[2] = ((*q_vec)[i])[2] + dt*((*p_vec)[i])[2];
		}
	}
}


void H_pot_tp(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec,double dt, unsigned int include_center) {
	unsigned int i,j;

	std::vector<std::vector<double> > temp_kick_x,temp_kick_y,temp_kick_z;
	double temp_kick_xd = 0,temp_kick_yd = 0,temp_kick_zd = 0;
	double temp_dist;
	std::vector<double> temp_p;

	for(i=0; i < (*q_vec).size()-1; i++) {
		if((*type_vec)[i] == 1) {
			temp_kick_x.push_back(temp_p);
			temp_kick_y.push_back(temp_p);
			temp_kick_z.push_back(temp_p);

			for(j=i+1; j < (*q_vec).size(); j++) {
				temp_dist = sqrt(( ((*q_vec)[i])[0] - ((*q_vec)[j])[0] )*( ((*q_vec)[i])[0] - ((*q_vec)[j])[0] ) + ( ((*q_vec)[i])[1] - ((*q_vec)[j])[1] )*( ((*q_vec)[i])[1] - ((*q_vec)[j])[1] ) + ( ((*q_vec)[i])[2] - ((*q_vec)[j])[2] )*( ((*q_vec)[i])[2] - ((*q_vec)[j])[2] ));
				temp_dist = pow(temp_dist,3);
				temp_kick_x[i].push_back( G*((*m_vector)[i])*((*m_vector)[j])*( ((*q_vec)[i])[0] - ((*q_vec)[j])[0] )/temp_dist);
				temp_kick_y[i].push_back( G*((*m_vector)[i])*((*m_vector)[j])*( ((*q_vec)[i])[1] - ((*q_vec)[j])[1] )/temp_dist);
				temp_kick_z[i].push_back( G*((*m_vector)[i])*((*m_vector)[j])*( ((*q_vec)[i])[2] - ((*q_vec)[j])[2] )/temp_dist);
			}

		}
		//REMOVED FOR SPEED
		/*else if((*type_vec)[i] != -1) { 
			temp_kick_x.push_back(temp_p);
			temp_kick_y.push_back(temp_p);
			temp_kick_z.push_back(temp_p);

			for(j=i+1; j < (*q_vec).size(); j++) {
				temp_kick_x[i].push_back( 0 );
				temp_kick_y[i].push_back( 0 );
				temp_kick_z[i].push_back( 0 );
			}
		}*/
	}

	for(i=(1-include_center); i < (*q_vec).size(); i++) {
		temp_kick_xd = 0; temp_kick_yd = 0; temp_kick_zd = 0;
		for(j=0; j < (*q_vec).size(); j++) {
			if((*type_vec)[j] == 1) { //ADDED FOR SPEED
			if(j < i) {
				temp_kick_xd = temp_kick_xd + temp_kick_x[j][i-j-1];
				temp_kick_yd = temp_kick_yd + temp_kick_y[j][i-j-1];
				temp_kick_zd = temp_kick_zd + temp_kick_z[j][i-j-1];
			}
			else if(j > i) {
				temp_kick_xd = temp_kick_xd + temp_kick_x[i][j-i-1];
				temp_kick_yd = temp_kick_yd + temp_kick_y[i][j-i-1];
				temp_kick_zd = temp_kick_zd + temp_kick_z[i][j-i-1];
			}
			}
		}

		((*p_vec)[i])[0] = ((*p_vec)[i])[0] + dt*temp_kick_xd;
		((*p_vec)[i])[1] = ((*p_vec)[i])[1] + dt*temp_kick_yd;
		((*p_vec)[i])[2] = ((*p_vec)[i])[2] + dt*temp_kick_zd;
	}


}


void H_pot_tp_DIS(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec,double dt, unsigned int include_center) {
	unsigned int i,j;

	std::vector<std::vector<double> > temp_kick_x,temp_kick_y,temp_kick_z;
	double temp_kick_xd = 0,temp_kick_yd = 0,temp_kick_zd = 0;
	double temp_dist;
	double EM_mag;
	double rho;
	double v_dot;
	double S;
	std::vector<double> temp_x;
	std::vector<double> temp_v;
	std::vector<double> temp_p;
	std::vector<double> debug_pr;double beta_1,beta_2;

	for(i=0; i < (*q_vec).size()-1; i++) {
		if((*type_vec)[i] == 1) {
			temp_kick_x.push_back(temp_p);
			temp_kick_y.push_back(temp_p);
			temp_kick_z.push_back(temp_p);

			for(j=i+1; j < (*q_vec).size(); j++) {
				temp_dist = sqrt(( ((*q_vec)[i])[0] - ((*q_vec)[j])[0] )*( ((*q_vec)[i])[0] - ((*q_vec)[j])[0] ) + ( ((*q_vec)[i])[1] - ((*q_vec)[j])[1] )*( ((*q_vec)[i])[1] - ((*q_vec)[j])[1] ) + ( ((*q_vec)[i])[2] - ((*q_vec)[j])[2] )*( ((*q_vec)[i])[2] - ((*q_vec)[j])[2] ));
				temp_dist = pow(temp_dist,3);
				temp_kick_x[i].push_back( G*((*m_vector)[i])*((*m_vector)[j])*( ((*q_vec)[i])[0] - ((*q_vec)[j])[0] )/temp_dist);
				temp_kick_y[i].push_back( G*((*m_vector)[i])*((*m_vector)[j])*( ((*q_vec)[i])[1] - ((*q_vec)[j])[1] )/temp_dist);
				temp_kick_z[i].push_back( G*((*m_vector)[i])*((*m_vector)[j])*( ((*q_vec)[i])[2] - ((*q_vec)[j])[2] )/temp_dist);
			}

		}
	}

	for(i=(1-include_center); i < (*q_vec).size(); i++) {
		temp_kick_xd = 0; temp_kick_yd = 0; temp_kick_zd = 0;
		for(j=0; j < (*q_vec).size(); j++) {
			if((*type_vec)[j] == 1) { //ADDED FOR SPEED
			if(j < i) {
				temp_kick_xd = temp_kick_xd + temp_kick_x[j][i-j-1];
				temp_kick_yd = temp_kick_yd + temp_kick_y[j][i-j-1];
				temp_kick_zd = temp_kick_zd + temp_kick_z[j][i-j-1];
			}
			else if(j > i) {
				temp_kick_xd = temp_kick_xd + temp_kick_x[i][j-i-1];
				temp_kick_yd = temp_kick_yd + temp_kick_y[i][j-i-1];
				temp_kick_zd = temp_kick_zd + temp_kick_z[i][j-i-1];
			}
			}
		}

		if((*type_vec)[i] == 0) { //EM effect
			temp_dist = abs_v((*q_vec)[i] - (*q_vec)[0]);
			S = pow(((*m_vector)[i]/rho)*0.75/PI,1.0/3.0);
			EM_mag = (sol_L*(S*S*PI))/(c0*4.0*PI*temp_dist*temp_dist);
			temp_v = (*p_vec)[i]/(*m_vector)[i];
			temp_x = (*q_vec)[i]/temp_dist;
			v_dot = dot_product(temp_v,temp_x);

			temp_kick_xd = temp_kick_xd + EM_mag*((1 - v_dot/c0)*temp_x[0] - temp_v[0]/c0);
			temp_kick_yd = temp_kick_yd + EM_mag*((1 - v_dot/c0)*temp_x[1] - temp_v[1]/c0);
			temp_kick_zd = temp_kick_zd + EM_mag*((1 - v_dot/c0)*temp_x[2] - temp_v[2]/c0);
		}

		((*p_vec)[i])[0] = ((*p_vec)[i])[0] + dt*temp_kick_xd;
		((*p_vec)[i])[1] = ((*p_vec)[i])[1] + dt*temp_kick_yd;
		((*p_vec)[i])[2] = ((*p_vec)[i])[2] + dt*temp_kick_zd;
	}


}



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

	//kepler(ro, vo, dt, G*(M+m)*1e-9, r, v);

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

	//kepler(ro, vo, dt, G*M*1e-9, r, v);

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
void H_dri(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec, double dt) {
	unsigned int i;
	double x_temp,y_temp,z_temp;
	x_temp = 0;
	y_temp = 0;
	z_temp = 0;

	for(i=1; i < (*p_vec).size(); i++) {
		if((*type_vec)[i] != -1) {
		x_temp = x_temp + ((*p_vec)[i])[0];
		y_temp = y_temp + ((*p_vec)[i])[1];
		z_temp = z_temp + ((*p_vec)[i])[2];
		}
	}

	x_temp = x_temp/(*m_vector)[0];
	y_temp = y_temp/(*m_vector)[0];
	z_temp = z_temp/(*m_vector)[0];

	for(i=1; i < (*q_vec).size(); i++) {
		if((*type_vec)[i] != -1) {
		((*q_vec)[i])[0] = ((*q_vec)[i])[0] + dt*x_temp;
		((*q_vec)[i])[1] = ((*q_vec)[i])[1] + dt*y_temp;
		((*q_vec)[i])[2] = ((*q_vec)[i])[2] + dt*z_temp;
		}
	}

}

void H_kep(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec, double dt, double sun_collision, double ejection_criteria) {
	unsigned int j;
	double mu, a, r_0, v_r0, Xi, z, r, alpha;
	double f, g, dotf, dotg, S, C;

	std::vector<double> temp_q, temp_v;
	temp_q.push_back(0);temp_q.push_back(0);temp_q.push_back(0);
	temp_v.push_back(0);temp_v.push_back(0);temp_v.push_back(0);

	for(j=1; j < (*q_vec).size(); j++) {
		if((*type_vec)[j] != -1) {
			mu = G*((*m_vector)[j] + (*m_vector)[0]);
			r_0 = abs_v((*q_vec)[j]);
			a = 1/(2/r_0 - abs_v((*p_vec)[j])/(mu*pow((*m_vector)[j],2)));
			alpha = 1/a;
			v_r0 = dot_product((*p_vec)[j], (*q_vec)[j])/(r_0*(*m_vector)[j]);

			Xi = universal_kepler_equation_newton(dt, 10e-8, a, r_0, v_r0, mu);
			//std::cout << a/AU << " " << v_r0 << std::endl;

			z = alpha*pow(Xi,2);
			S = stumpff(3,6,z);
			C = stumpff(2,6,z);

			f = 1 - pow(Xi,2)/r_0*C;
			g = dt - 1/sqrt(mu)*pow(Xi,3)*S;

			temp_q[0] = f*((*q_vec)[j])[0] + g*((*p_vec)[j])[0]/(*m_vector)[j];
			temp_q[1] = f*((*q_vec)[j])[1] + g*((*p_vec)[j])[1]/(*m_vector)[j];
			temp_q[2] = f*((*q_vec)[j])[2] + g*((*p_vec)[j])[2]/(*m_vector)[j];

			r = abs_v(temp_q);

			dotf = sqrt(mu)/(r*r_0)*(alpha*pow(Xi,3)*S - Xi);
			dotg = 1 - pow(Xi,2)/r*C;

			temp_v[0] = dotf*((*q_vec)[j])[0] + dotg*((*p_vec)[j])[0]/(*m_vector)[j];
			temp_v[1] = dotf*((*q_vec)[j])[1] + dotg*((*p_vec)[j])[1]/(*m_vector)[j];
			temp_v[2] = dotf*((*q_vec)[j])[2] + dotg*((*p_vec)[j])[2]/(*m_vector)[j];

			((*q_vec)[j])[0] = temp_q[0];
			((*q_vec)[j])[1] = temp_q[1];
			((*q_vec)[j])[2] = temp_q[2];

			((*p_vec)[j])[0] = temp_v[0]*(*m_vector)[j];
			((*p_vec)[j])[1] = temp_v[1]*(*m_vector)[j];
			((*p_vec)[j])[2] = temp_v[2]*(*m_vector)[j];
		}
	}

}



void H_mod_mid(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec, double dt, double n, double sun_collision, double ejection_criteria) {
	unsigned int i,j;
	double h = dt/n;
	double temp_dist;
	double temp_kick_x,temp_kick_y,temp_kick_z;
	double temp_drift_x,temp_drift_y,temp_drift_z;
	std::vector<double> temp;
	temp.push_back(0);
	temp.push_back(0);
	temp.push_back(0);

	std::vector<std::vector<double> > q_temp;
	std::vector<std::vector<double> > p_temp;

	std::vector<std::vector<double> > q_temp_2;
	std::vector<std::vector<double> > p_temp_2;

	std::vector<std::vector<double> > q_temp_3;
	std::vector<std::vector<double> > p_temp_3;


	for(i=0; i < (*p_vec).size(); i++) {
		q_temp.push_back(temp);
		p_temp.push_back(temp);
		q_temp_2.push_back(temp);
		p_temp_2.push_back(temp);
		q_temp_3.push_back(temp);
		p_temp_3.push_back(temp);

		(q_temp[i])[0] = ((*q_vec)[i])[0];
		(q_temp[i])[1] = ((*q_vec)[i])[1];
		(q_temp[i])[2] = ((*q_vec)[i])[2];

		(p_temp[i])[0] = ((*p_vec)[i])[0];
		(p_temp[i])[1] = ((*p_vec)[i])[1];
		(p_temp[i])[2] = ((*p_vec)[i])[2];
	}

	for(i=1; i < (*p_vec).size(); i++) {
		if((*type_vec)[i] != -1) {
				temp_drift_x = ((*p_vec)[i])[0]/(*m_vector)[i];
				temp_drift_y = ((*p_vec)[i])[1]/(*m_vector)[i];
				temp_drift_z = ((*p_vec)[i])[2]/(*m_vector)[i];

				(q_temp_2[i])[0] = ((*q_vec)[i])[0] + h*temp_drift_x;
				(q_temp_2[i])[1] = ((*q_vec)[i])[1] + h*temp_drift_y;
				(q_temp_2[i])[2] = ((*q_vec)[i])[2] + h*temp_drift_z;

				temp_dist = pow(sqrt(( (q_temp_2[i])[0] - ((*q_vec)[0])[0] )*( (q_temp_2[i])[0] - ((*q_vec)[0])[0] ) + ( (q_temp_2[i])[1] - ((*q_vec)[0])[1] )*( (q_temp_2[i])[1] - ((*q_vec)[0])[1] ) + ( (q_temp_2[i])[2] - ((*q_vec)[0])[2] )*( (q_temp_2[i])[2] - ((*q_vec)[0])[2] )),3);
				temp_kick_x = -G*((*m_vector)[i])*((*m_vector)[0])*( (q_temp_2[i])[0] - ((*q_vec)[0])[0] )/temp_dist;
				temp_kick_y = -G*((*m_vector)[i])*((*m_vector)[0])*( (q_temp_2[i])[1] - ((*q_vec)[0])[1] )/temp_dist;
				temp_kick_z = -G*((*m_vector)[i])*((*m_vector)[0])*( (q_temp_2[i])[2] - ((*q_vec)[0])[2] )/temp_dist;

				(p_temp_2[i])[0] = ((*p_vec)[i])[0] + h*temp_kick_x;
				(p_temp_2[i])[1] = ((*p_vec)[i])[1] + h*temp_kick_y;
				(p_temp_2[i])[2] = ((*p_vec)[i])[2] + h*temp_kick_z;
		}
	}

	for(j=0; j < n-1; j++) {
		for(i=1; i < (*p_vec).size(); i++) {
			if((*type_vec)[i] != -1) {
			temp_drift_x = (p_temp_2[i])[0]/(*m_vector)[i];
			temp_drift_y = (p_temp_2[i])[1]/(*m_vector)[i];
			temp_drift_z = (p_temp_2[i])[2]/(*m_vector)[i];

			(q_temp_3[i])[0] = (q_temp[i])[0] + 2*h*temp_drift_x;
			(q_temp_3[i])[1] = (q_temp[i])[1] + 2*h*temp_drift_y;
			(q_temp_3[i])[2] = (q_temp[i])[2] + 2*h*temp_drift_z;

			temp_dist = pow(sqrt(( (q_temp_3[i])[0] - ((*q_vec)[0])[0] )*( (q_temp_3[i])[0] - ((*q_vec)[0])[0] ) + ( (q_temp_3[i])[1] - ((*q_vec)[0])[1] )*( (q_temp_3[i])[1] - ((*q_vec)[0])[1] ) + ( (q_temp_3[i])[2] - ((*q_vec)[0])[2] )*( (q_temp_3[i])[2] - ((*q_vec)[0])[2] )),3);
			temp_kick_x = -G*((*m_vector)[i])*((*m_vector)[0])*( (q_temp_3[i])[0] - ((*q_vec)[0])[0] )/temp_dist;
			temp_kick_y = -G*((*m_vector)[i])*((*m_vector)[0])*( (q_temp_3[i])[1] - ((*q_vec)[0])[1] )/temp_dist;
			temp_kick_z = -G*((*m_vector)[i])*((*m_vector)[0])*( (q_temp_3[i])[2] - ((*q_vec)[0])[2] )/temp_dist;

			(p_temp_3[i])[0] = (p_temp[i])[0] + 2*h*temp_kick_x;
			(p_temp_3[i])[1] = (p_temp[i])[1] + 2*h*temp_kick_y;
			(p_temp_3[i])[2] = (p_temp[i])[2] + 2*h*temp_kick_z;

			(q_temp[i])[0] = (q_temp_2[i])[0];
			(q_temp[i])[1] = (q_temp_2[i])[1];
			(q_temp[i])[2] = (q_temp_2[i])[2];

			(p_temp[i])[0] = (p_temp_2[i])[0];
			(p_temp[i])[1] = (p_temp_2[i])[1];
			(p_temp[i])[2] = (p_temp_2[i])[2];

			(q_temp_2[i])[0] = (q_temp_3[i])[0];
			(q_temp_2[i])[1] = (q_temp_3[i])[1];
			(q_temp_2[i])[2] = (q_temp_3[i])[2];

			(p_temp_2[i])[0] = (p_temp_3[i])[0];
			(p_temp_2[i])[1] = (p_temp_3[i])[1];
			(p_temp_2[i])[2] = (p_temp_3[i])[2];
			}
		}
	}

	for(i=1; i < (*p_vec).size(); i++) {
		if((*type_vec)[i] != -1) {
		temp_drift_x = (p_temp_2[i])[0]/(*m_vector)[i];
		temp_drift_y = (p_temp_2[i])[1]/(*m_vector)[i];
		temp_drift_z = (p_temp_2[i])[2]/(*m_vector)[i];

		((*q_vec)[i])[0] = 0.5*((q_temp_2[i])[0] + (q_temp[i])[0] + h*temp_drift_x);
		((*q_vec)[i])[1] = 0.5*((q_temp_2[i])[1] + (q_temp[i])[1] + h*temp_drift_y);
		((*q_vec)[i])[2] = 0.5*((q_temp_2[i])[2] + (q_temp[i])[2] + h*temp_drift_z);

		temp_dist = sqrt(( (q_temp_2[i])[0] - ((*q_vec)[0])[0] )*( (q_temp_2[i])[0] - ((*q_vec)[0])[0] ) + ( (q_temp_2[i])[1] - ((*q_vec)[0])[1] )*( (q_temp_2[i])[1] - ((*q_vec)[0])[1] ) + ( (q_temp_2[i])[2] - ((*q_vec)[0])[2] )*( (q_temp_2[i])[2] - ((*q_vec)[0])[2] ));
		if((temp_dist < sun_collision) | (temp_dist > ejection_criteria)) {
			(*type_vec)[i] = -1;
		}
		temp_dist = pow(temp_dist,3);
		temp_kick_x = -G*((*m_vector)[i])*((*m_vector)[0])*( (q_temp_2[i])[0] - ((*q_vec)[0])[0] )/temp_dist;
		temp_kick_y = -G*((*m_vector)[i])*((*m_vector)[0])*( (q_temp_2[i])[1] - ((*q_vec)[0])[1] )/temp_dist;
		temp_kick_z = -G*((*m_vector)[i])*((*m_vector)[0])*( (q_temp_2[i])[2] - ((*q_vec)[0])[2] )/temp_dist;

		((*p_vec)[i])[0] = 0.5*((p_temp_2[i])[0] + (p_temp[i])[0] + h*temp_kick_x);
		((*p_vec)[i])[1] = 0.5*((p_temp_2[i])[1] + (p_temp[i])[1] + h*temp_kick_y);
		((*p_vec)[i])[2] = 0.5*((p_temp_2[i])[2] + (p_temp[i])[2] + h*temp_kick_z);
		}
	}

}

void H_euler(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec, double dt, double sun_collision, double ejection_criteria) {
	unsigned int i,j;
	double temp_dist;
	double temp_kick_x,temp_kick_y,temp_kick_z;
	double temp_drift_x,temp_drift_y,temp_drift_z;

	for(i=0; i < (*p_vec).size(); i++) {
		if((*type_vec)[i] != -1) {
		temp_drift_x = ((*p_vec)[i])[0]/(*m_vector)[i];
		temp_drift_y = ((*p_vec)[i])[1]/(*m_vector)[i];
		temp_drift_z = ((*p_vec)[i])[2]/(*m_vector)[i];

		((*q_vec)[i])[0] = ((*q_vec)[i])[0] + dt*temp_drift_x;
		((*q_vec)[i])[1] = ((*q_vec)[i])[1] + dt*temp_drift_y;
		((*q_vec)[i])[2] = ((*q_vec)[i])[2] + dt*temp_drift_z;
		temp_kick_x = 0;
		temp_kick_y = 0;
		temp_kick_z = 0;
		for(j=0; j < (*p_vec).size(); j++) {
			if(((*type_vec)[j] != -1) & (j != i)) {
				temp_dist = sqrt(( ((*q_vec)[j])[0] - ((*q_vec)[i])[0] )*( ((*q_vec)[j])[0] - ((*q_vec)[i])[0] ) + ( ((*q_vec)[j])[1] - ((*q_vec)[i])[1] )*( ((*q_vec)[j])[1] - ((*q_vec)[i])[1] ) + ( ((*q_vec)[j])[2] - ((*q_vec)[i])[2] )*( ((*q_vec)[j])[2] - ((*q_vec)[i])[2] ));
				if((j == 0) & ((temp_dist < sun_collision) | (temp_dist > ejection_criteria)) ) {
					(*type_vec)[i] = -1;
				}
				temp_dist = pow(temp_dist,3);
				temp_kick_x = temp_kick_x -G*((*m_vector)[i])*((*m_vector)[j])*( ((*q_vec)[i])[0] - ((*q_vec)[j])[0] )/temp_dist;
				temp_kick_y = temp_kick_y -G*((*m_vector)[i])*((*m_vector)[j])*( ((*q_vec)[i])[1] - ((*q_vec)[j])[1] )/temp_dist;
				temp_kick_z = temp_kick_z -G*((*m_vector)[i])*((*m_vector)[j])*( ((*q_vec)[i])[2] - ((*q_vec)[j])[2] )/temp_dist;
			}
		}
		((*p_vec)[i])[0] = ((*p_vec)[i])[0] + dt*temp_kick_x;
		((*p_vec)[i])[1] = ((*p_vec)[i])[1] + dt*temp_kick_y;
		((*p_vec)[i])[2] = ((*p_vec)[i])[2] + dt*temp_kick_z;
		}
	}

}
