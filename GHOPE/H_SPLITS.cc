
#include <iostream>
#include <vector>
#include <math.h>

#include "resources.hh"
#include "functions.hh"

//Avrage step execution             : 6.13856e-05 s
//Pre-allocate memory
//Avrage step execution             : 3.01266e-05 s
//Changed to inline cube
//Avrage step execution             : 2.93164e-05 s
//Changed to inline sqrt
//Avrage step execution             : 2.87688e-05 s
//Changed to newly designed data type
//Avrage step execution             : 1.91815e-05 s
//Introduce iterators
//
void H_pot(	tensor<double> & q_vec, 
			tensor<double> & p_vec, 
			tensor<double>& m_vector, 
			tensor<int> & type_vec, 
			double dt, 
			tensor<double>& interaction_data) {

	unsigned int i,j,im,jm;
	double temp_kick_xd = 0,temp_kick_yd = 0,temp_kick_zd = 0;
	double temp_dist;

	for(i=0; i < q_vec.size(0)-1; i++) {
		if(type_vec(i) == 1) {
			for(j=i+1; j < q_vec.size(0); j++) {
				temp_dist = inline_CUB(sqrt(
					inline_SQR( q_vec(i,0) - q_vec(j,0) ) + 
					inline_SQR( q_vec(i,1) - q_vec(j,1) ) + 
					inline_SQR( q_vec(i,2) - q_vec(j,2) )
					));

				interaction_data(i,j,0) = ( G*(m_vector(i))*(m_vector(j))*( q_vec(i,0) - q_vec(j,0) )/temp_dist);
				interaction_data(i,j,1) = ( G*(m_vector(i))*(m_vector(j))*( q_vec(i,1) - q_vec(j,1) )/temp_dist);
				interaction_data(i,j,2) = ( G*(m_vector(i))*(m_vector(j))*( q_vec(i,2) - q_vec(j,2) )/temp_dist);
			}

		}
	}

	for(i=0; i < q_vec.size(0); i++) {
		temp_kick_xd = 0; 
		temp_kick_yd = 0; 
		temp_kick_zd = 0;
		for(j=0; j < q_vec.size(0); j++) {
			if(type_vec(j) == 1 && j != i) {
				im = inline_MIN(i,j);
				jm = inline_MAX(i,j);
				temp_kick_xd += interaction_data(im,jm,0);
				temp_kick_yd += interaction_data(im,jm,1);
				temp_kick_zd += interaction_data(im,jm,2);
			}
		}

		p_vec(i,0) += dt*temp_kick_xd;
		p_vec(i,1) += dt*temp_kick_yd;
		p_vec(i,2) += dt*temp_kick_zd;
		
	}
	

}

void H_kin(	tensor<double> & q_vec, 
			tensor<double> & p_vec, 
			tensor<double>& m_vector, 
			tensor<int> & type_vec, 
			double dt) {
	unsigned int i;

	for(i=0; i < q_vec.size(0); i++) {
		if(type_vec(i) != -1) {
			q_vec(i,0) += dt*p_vec(i,0)/m_vector(i);
			q_vec(i,1) += dt*p_vec(i,1)/m_vector(i);
			q_vec(i,2) += dt*p_vec(i,2)/m_vector(i);
		}
	}
}
