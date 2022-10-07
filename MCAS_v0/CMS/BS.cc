/*
 * BS.cc
 *
 *  Created on: Jun 17, 2015
 *      Author: dankas
 */

#include <vector>
#include <math.h>
#include <algorithm>    // std::random_shuffle

#include "define.hh"
#include "functions.hh"



void BS_step(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec, std::vector<unsigned int> N, double dt) {
	unsigned int i,j,k;double rho = 0;
	std::vector<std::vector<double> > Z,P,X;
	std::vector<std::vector<double> > Zinit;
	std::vector<std::vector<double> > all_kicks_init;
	std::vector<std::vector<double> > all_kicks;
	std::vector<double> H;
	int TE = 0;
	unsigned int OBJS = (*q_vec).size();
	unsigned int VARS = OBJS*3;
	unsigned int SIZE = VARS*2;

	P = fill_mat(SIZE,N.size(),0.0);
	Zinit.push_back(H);
	Zinit[0].resize(SIZE);
	for(i=0; i < OBJS; i++) {
		Zinit[0][i*3]   = (*q_vec)[i][0];
		Zinit[0][i*3+1] = (*q_vec)[i][1];
		Zinit[0][i*3+2] = (*q_vec)[i][2];

		Zinit[0][VARS + i*3]   = (*p_vec)[i][0];
		Zinit[0][VARS + i*3+1] = (*p_vec)[i][1];
		Zinit[0][VARS + i*3+2] = (*p_vec)[i][2];
	}

	all_kicks_init = dZ(&Zinit, m_vector, type_vec, 0, rho);

	for(k=0; k < N.size(); k++) {
		H.push_back(dt/N[k]);
		Z = fill_mat(N[k]+1,SIZE,0.0);
		Z[0] = Zinit[0];

		for(i=0; i < OBJS; i++) {
			Z[1][i*3]   = (*q_vec)[i][0] + H[k]*(Z[0][VARS + i*3  ]/((*m_vector)[i]));
			Z[1][i*3+1] = (*q_vec)[i][1] + H[k]*(Z[0][VARS + i*3+1]/((*m_vector)[i]));
			Z[1][i*3+2] = (*q_vec)[i][2] + H[k]*(Z[0][VARS + i*3+2]/((*m_vector)[i]));

			Z[1][VARS + i*3]   = (*p_vec)[i][0] + H[k]*all_kicks_init[i][0];
			Z[1][VARS + i*3+1] = (*p_vec)[i][1] + H[k]*all_kicks_init[i][1];
			Z[1][VARS + i*3+2] = (*p_vec)[i][2] + H[k]*all_kicks_init[i][2];
		}

		for(j=1; j < N[k]; j++) {

			all_kicks = dZ(&Z, m_vector, type_vec, j, rho);

			for(i=0; i < OBJS; i++) {
				Z[j+1][i*3  ] = Z[j-1][i*3  ] + 2.0*H[k]*( Z[j][VARS + i*3  ] /(*m_vector)[i]);
				Z[j+1][i*3+1] = Z[j-1][i*3+1] + 2.0*H[k]*( Z[j][VARS + i*3+1] /(*m_vector)[i]);
				Z[j+1][i*3+2] = Z[j-1][i*3+2] + 2.0*H[k]*( Z[j][VARS + i*3+2] /(*m_vector)[i]);

				Z[j+1][VARS + i*3  ] = Z[j-1][VARS + i*3  ] + 2.0*H[k]*all_kicks[i][0];
				Z[j+1][VARS + i*3+1] = Z[j-1][VARS + i*3+1] + 2.0*H[k]*all_kicks[i][1];
				Z[j+1][VARS + i*3+2] = Z[j-1][VARS + i*3+2] + 2.0*H[k]*all_kicks[i][2];
			}
		}

		all_kicks = dZ(&Z, m_vector, type_vec, N[k], rho);
		for(i=0; i < OBJS; i++) {
			P[i*3  ][k] = 0.5*(Z[N[k]][i*3  ] + Z[N[k]-2][i*3  ] + H[k]*( Z[N[k]][VARS + i*3  ] /(*m_vector)[i]) );
			P[i*3+1][k] = 0.5*(Z[N[k]][i*3+1] + Z[N[k]-2][i*3+1] + H[k]*( Z[N[k]][VARS + i*3+1] /(*m_vector)[i]) );
			P[i*3+2][k] = 0.5*(Z[N[k]][i*3+2] + Z[N[k]-2][i*3+2] + H[k]*( Z[N[k]][VARS + i*3+2] /(*m_vector)[i]) );

			P[VARS + i*3  ][k] = 0.5*(Z[N[k]][VARS + i*3  ] + Z[N[k]-2][VARS + i*3  ] + H[k]*( all_kicks[i][0] ) );
			P[VARS + i*3+1][k] = 0.5*(Z[N[k]][VARS + i*3+1] + Z[N[k]-2][VARS + i*3+1] + H[k]*( all_kicks[i][1] ) );
			P[VARS + i*3+2][k] = 0.5*(Z[N[k]][VARS + i*3+2] + Z[N[k]-2][VARS + i*3+2] + H[k]*( all_kicks[i][2] ) );
		}

	}

	for(i=0; i < OBJS; i++) {
		((*q_vec)[i])[0] = neville(&H,&P[i*3  ],0.0);
		((*q_vec)[i])[1] = neville(&H,&P[i*3+1],0.0);
		((*q_vec)[i])[2] = neville(&H,&P[i*3+2],0.0);

		((*p_vec)[i])[0] = neville(&H,&P[VARS + i*3  ],0.0);
		((*p_vec)[i])[1] = neville(&H,&P[VARS + i*3+1],0.0);
		((*p_vec)[i])[2] = neville(&H,&P[VARS + i*3+2],0.0);
	}

}


std::vector<std::vector<double> > dZ(std::vector<std::vector<double> >* Z, std::vector<double>* m_vector, std::vector<int>* type_vec, unsigned int ind, double rho) {
	std::vector<std::vector<double> > temp_kick_x,temp_kick_y,temp_kick_z;
	std::vector<std::vector<double> > all_kicks;
	double temp_kick_xd = 0,temp_kick_yd = 0,temp_kick_zd = 0;
	double temp_dist;
	double EM_mag;
	double v_dot;
	double S;
	std::vector<double> temp_x;
	std::vector<double> temp_v;
	std::vector<double> temp_p;

	unsigned int i,j;


	temp_x.resize(3);
	temp_v.resize(3);
	
	unsigned int OBJS = (*type_vec).size();

	for(i=0; i < OBJS-1; i++) {
		if((*type_vec)[i] == 1) {
			temp_kick_x.push_back(temp_p);
			temp_kick_y.push_back(temp_p);
			temp_kick_z.push_back(temp_p);

			for(j=i+1; j < OBJS; j++) {
				temp_dist = sqrt(( (*Z)[ind][i*3] - (*Z)[ind][j*3] )*( (*Z)[ind][i*3] - (*Z)[ind][j*3] ) + ( (*Z)[ind][i*3+1] - (*Z)[ind][j*3+1] )*( (*Z)[ind][i*3+1] - (*Z)[ind][j*3+1] ) + ( (*Z)[ind][i*3+2] - (*Z)[ind][j*3+2] )*( (*Z)[ind][i*3+2] - (*Z)[ind][j*3+2] ));
				temp_dist = pow(temp_dist,3);
				temp_kick_x[i].push_back( G*((*m_vector)[i])*((*m_vector)[j])*( (*Z)[ind][i*3  ] - (*Z)[ind][j*3  ] )/temp_dist);
				temp_kick_y[i].push_back( G*((*m_vector)[i])*((*m_vector)[j])*( (*Z)[ind][i*3+1] - (*Z)[ind][j*3+1] )/temp_dist);
				temp_kick_z[i].push_back( G*((*m_vector)[i])*((*m_vector)[j])*( (*Z)[ind][i*3+2] - (*Z)[ind][j*3+2] )/temp_dist);
			}

		}
	}

	for(i=0; i < OBJS; i++) {
		temp_kick_xd = 0; temp_kick_yd = 0; temp_kick_zd = 0;
		for(j=0; j < OBJS; j++) {
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
			
			temp_x[0] = (*Z)[ind][i*3  ]-(*Z)[ind][0];
			temp_x[1] = (*Z)[ind][i*3+1]-(*Z)[ind][1];
			temp_x[2] = (*Z)[ind][i*3+2]-(*Z)[ind][2];


			temp_dist = abs_v(temp_x);
			S = pow(((*m_vector)[i]/rho)*0.75/PI,1.0/3.0);
			EM_mag = (sol_L*(S*S*PI))/(c0*4.0*PI*temp_dist*temp_dist);

			temp_v[0] = ((*Z)[ind][(OBJS+i)*3  ]/(*m_vector)[i]);
			temp_v[1] = ((*Z)[ind][(OBJS+i)*3+1]/(*m_vector)[i]);
			temp_v[2] = ((*Z)[ind][(OBJS+i)*3+2]/(*m_vector)[i]);

			temp_x[0] = (*Z)[ind][i*3  ]/temp_dist;
			temp_x[1] = (*Z)[ind][i*3+1]/temp_dist;
			temp_x[2] = (*Z)[ind][i*3+2]/temp_dist;

			v_dot = dot_product(temp_v,temp_x);

			temp_kick_xd = temp_kick_xd + EM_mag*((1 - v_dot/c0)*temp_x[0] - temp_v[0]/c0);
			temp_kick_yd = temp_kick_yd + EM_mag*((1 - v_dot/c0)*temp_x[1] - temp_v[1]/c0);
			temp_kick_zd = temp_kick_zd + EM_mag*((1 - v_dot/c0)*temp_x[2] - temp_v[2]/c0);
		}

		temp_p.clear();
		temp_p.push_back(temp_kick_xd);
		temp_p.push_back(temp_kick_yd);
		temp_p.push_back(temp_kick_zd);
		all_kicks.push_back(temp_p);
	}
	return all_kicks;
}
