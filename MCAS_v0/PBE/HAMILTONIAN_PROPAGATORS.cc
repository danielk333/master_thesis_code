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
/*
 std::vector<double> barycentre(std::vector<std::vector<double> > *q_vec, std::vector<double> *m_vector);
 std::vector<double> barycentre_v(std::vector<std::vector<double> > *p_vec, std::vector<double> *m_vector);

std::vector<double> barycentre(std::vector<std::vector<double> > *q_vec, std::vector<double> *m_vector) {
	std::vector<double> C;
	C.push_back(0);C.push_back(0);C.push_back(0);
	unsigned int i;
	double m_tot = 0;
	for(i=0; i < (*q_vec).size(); i++) {
		C[0] = C[0] + (*m_vector)[i]*(*q_vec[i])[0];
		C[1] = C[1] + (*m_vector)[i]*(*q_vec[i])[1];
		C[2] = C[2] + (*m_vector)[i]*(*q_vec[i])[2];
		m_tot = m_tot + (*m_vector)[i];
	}
	C[0] = C[0]/m_tot;
	C[1] = C[1]/m_tot;
	C[2] = C[2]/m_tot;

	return C;
}

std::vector<double> barycentre_v(std::vector<std::vector<double> > *p_vec, std::vector<double> *m_vector) {
	std::vector<double> C;
	double m_temp = 0;
	C.push_back(0);C.push_back(0);C.push_back(0);
	unsigned int i;
	for(i=0; i < (*p_vec).size(); i++) {
		C[0] = C[0] + (*p_vec)[i][0];
		C[1] = C[1] + (*p_vec)[i][1];
		C[2] = C[2] + (*p_vec)[i][2];
		m_temp = m_temp + (*m_vector)[i];
	}

	C[0] = C[0]/m_temp;
	C[1] = C[1]/m_temp;
	C[2] = C[2]/m_temp;
	return C;
}



void barycentric_p(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vec) {
	int i;
	std::vector<double> bary;
	std::vector<double> bary_v;

	bary_v = barycentre(p_vec,m_vec);

	for(i=(*q_vec).size()-1; i >= 0; i--) {
		((*q_vec)[i])[0] = ((*q_vec)[i])[0] - ((*q_vec)[0])[0];
		((*q_vec)[i])[1] = ((*q_vec)[i])[1] - ((*q_vec)[0])[1];
		((*q_vec)[i])[2] = ((*q_vec)[i])[2] - ((*q_vec)[0])[2];

		((*p_vec)[i])[0] = ((*p_vec)[i])[0] - ((*p_vec)[0])[0]*(((*m_vec)[i])/((*m_vec)[0]));
		((*p_vec)[i])[1] = ((*p_vec)[i])[1] - ((*p_vec)[0])[1]*(((*m_vec)[i])/((*m_vec)[0]));
		((*p_vec)[i])[2] = ((*p_vec)[i])[2] - ((*p_vec)[0])[2]*(((*m_vec)[i])/((*m_vec)[0]));
	}
}
*/

void restore_coord_qp(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vec, std::vector<double> coord) {
	int i;

	for(i=(*q_vec).size()-1; i >= 0; i--) {
		((*q_vec)[i])[0] = ((*q_vec)[i])[0] + coord[0];
		((*q_vec)[i])[1] = ((*q_vec)[i])[1] + coord[1];
		((*q_vec)[i])[2] = ((*q_vec)[i])[2] + coord[2];

		((*p_vec)[i])[0] = ((*p_vec)[i])[0] + ((*m_vec)[i])*coord[3];
		((*p_vec)[i])[1] = ((*p_vec)[i])[1] + ((*m_vec)[i])*coord[4];
		((*p_vec)[i])[2] = ((*p_vec)[i])[2] + ((*m_vec)[i])*coord[5];
	}

}
std::vector<double> temp_heliocentric_qp(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vec) {
	int i;
	std::vector<double> ret;
	ret.push_back(((*q_vec)[0])[0]);
	ret.push_back(((*q_vec)[0])[1]);
	ret.push_back(((*q_vec)[0])[2]);
	ret.push_back(((*p_vec)[0])[0]/((*m_vec)[0]));
	ret.push_back(((*p_vec)[0])[1]/((*m_vec)[0]));
	ret.push_back(((*p_vec)[0])[2]/((*m_vec)[0]));
	
	for(i=(*q_vec).size()-1; i >= 0; i--) {
		((*q_vec)[i])[0] = ((*q_vec)[i])[0] - ((*q_vec)[0])[0];
		((*q_vec)[i])[1] = ((*q_vec)[i])[1] - ((*q_vec)[0])[1];
		((*q_vec)[i])[2] = ((*q_vec)[i])[2] - ((*q_vec)[0])[2];

		((*p_vec)[i])[0] = ((*p_vec)[i])[0] - ((*p_vec)[0])[0]*(((*m_vec)[i])/((*m_vec)[0]));
		((*p_vec)[i])[1] = ((*p_vec)[i])[1] - ((*p_vec)[0])[1]*(((*m_vec)[i])/((*m_vec)[0]));
		((*p_vec)[i])[2] = ((*p_vec)[i])[2] - ((*p_vec)[0])[2]*(((*m_vec)[i])/((*m_vec)[0]));
	}
	return ret;
}

void heliocentric_q(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vec) {
	int i;
	for(i=(*q_vec).size()-1; i >= 0; i--) {
		((*q_vec)[i])[0] = ((*q_vec)[i])[0] - ((*q_vec)[0])[0];
		((*q_vec)[i])[1] = ((*q_vec)[i])[1] - ((*q_vec)[0])[1];
		((*q_vec)[i])[2] = ((*q_vec)[i])[2] - ((*q_vec)[0])[2];
	}
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

void H_pot(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec,double dt, unsigned int include_center) {
	unsigned int i,j;

	std::vector<std::vector<double> > temp_kick_x,temp_kick_y,temp_kick_z;
	double temp_kick_xd = 0,temp_kick_yd = 0,temp_kick_zd = 0;
	double temp_dist;
	std::vector<double> temp_p;

	for(i=0; i < (*q_vec).size()-1; i++) {
		if((*type_vec)[i] != -1) {
		temp_kick_x.push_back(temp_p);
		temp_kick_y.push_back(temp_p);
		temp_kick_z.push_back(temp_p);

		for(j=i+1; j < (*q_vec).size(); j++) {
				temp_dist = sqrt(( ((*q_vec)[i])[0] - ((*q_vec)[j])[0] )*( ((*q_vec)[i])[0] - ((*q_vec)[j])[0] ) + ( ((*q_vec)[i])[1] - ((*q_vec)[j])[1] )*( ((*q_vec)[i])[1] - ((*q_vec)[j])[1] ) + ( ((*q_vec)[i])[2] - ((*q_vec)[j])[2] )*( ((*q_vec)[i])[2] - ((*q_vec)[j])[2] ));
				temp_dist = pow(temp_dist,3);
				temp_kick_x[i].push_back( Ggrav*((*m_vector)[i])*((*m_vector)[j])*( ((*q_vec)[i])[0] - ((*q_vec)[j])[0] )/temp_dist);
				temp_kick_y[i].push_back( Ggrav*((*m_vector)[i])*((*m_vector)[j])*( ((*q_vec)[i])[1] - ((*q_vec)[j])[1] )/temp_dist);
				temp_kick_z[i].push_back( Ggrav*((*m_vector)[i])*((*m_vector)[j])*( ((*q_vec)[i])[2] - ((*q_vec)[j])[2] )/temp_dist);
		}

		}
	}

	for(i=(1-include_center); i < (*q_vec).size(); i++) {
		temp_kick_xd = 0; temp_kick_yd = 0; temp_kick_zd = 0;
		for(j=0; j < (*q_vec).size(); j++) {
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

		((*p_vec)[i])[0] = ((*p_vec)[i])[0] + dt*temp_kick_xd;
		((*p_vec)[i])[1] = ((*p_vec)[i])[1] + dt*temp_kick_yd;
		((*p_vec)[i])[2] = ((*p_vec)[i])[2] + dt*temp_kick_zd;
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
				temp_kick_x[i].push_back( Ggrav*((*m_vector)[i])*((*m_vector)[j])*( ((*q_vec)[i])[0] - ((*q_vec)[j])[0] )/temp_dist);
				temp_kick_y[i].push_back( Ggrav*((*m_vector)[i])*((*m_vector)[j])*( ((*q_vec)[i])[1] - ((*q_vec)[j])[1] )/temp_dist);
				temp_kick_z[i].push_back( Ggrav*((*m_vector)[i])*((*m_vector)[j])*( ((*q_vec)[i])[2] - ((*q_vec)[j])[2] )/temp_dist);
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

void H_kin(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec, double dt) {
	unsigned int i;

	for(i=0; i < (*q_vec).size(); i++) {
		if((*type_vec)[i] != -1) {
		((*q_vec)[i])[0] = ((*q_vec)[i])[0] + dt*((*p_vec)[i])[0]/((*m_vector)[i]);
		((*q_vec)[i])[1] = ((*q_vec)[i])[1] + dt*((*p_vec)[i])[1]/((*m_vector)[i]);
		((*q_vec)[i])[2] = ((*q_vec)[i])[2] + dt*((*p_vec)[i])[2]/((*m_vector)[i]);
		}
	}
}


void H_pot_tp_DIS(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec,double dt, unsigned int include_center, double rho) {
	unsigned int i,j;

	std::vector<std::vector<double> > temp_kick_x,temp_kick_y,temp_kick_z;
	double temp_kick_xd = 0,temp_kick_yd = 0,temp_kick_zd = 0;
	double temp_dist;
	double EM_mag;
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
				temp_kick_x[i].push_back( Ggrav*((*m_vector)[i])*((*m_vector)[j])*( ((*q_vec)[i])[0] - ((*q_vec)[j])[0] )/temp_dist);
				temp_kick_y[i].push_back( Ggrav*((*m_vector)[i])*((*m_vector)[j])*( ((*q_vec)[i])[1] - ((*q_vec)[j])[1] )/temp_dist);
				temp_kick_z[i].push_back( Ggrav*((*m_vector)[i])*((*m_vector)[j])*( ((*q_vec)[i])[2] - ((*q_vec)[j])[2] )/temp_dist);
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
			temp_dist = abs_v(add_v_v((*q_vec)[i],multiply_v(-1.0,(*q_vec)[0])));
			S = pow(((*m_vector)[i]/rho)*0.75/PI,1.0/3.0);
			EM_mag = (sol_L*(S*S*PI))/(c0*4.0*PI*temp_dist*temp_dist);
			temp_v = multiply_v(1.0/(*m_vector)[i],(*p_vec)[i]);
			temp_x = multiply_v(1.0/temp_dist,(*q_vec)[i]);
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

