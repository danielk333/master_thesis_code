

#include <vector>
#include <math.h>

#include "resources.hh"
#include "functions.hh"


tensor<double> barycentre(tensor<double> &q_vec, tensor<double> &m_vector) {
	int D[] = {3};
	tensor<double> C(D,1,0.0);
	double m_tot = 0;
	for(int i=0; i < q_vec.size(0); i++) {
		C(0) += m_vector(i)*q_vec(i,0);
		C(1) += m_vector(i)*q_vec(i,1);
		C(2) += m_vector(i)*q_vec(i,2);
		m_tot += m_vector(i);
	}
	C(0) /= m_tot;
	C(1) /= m_tot;
	C(2) /= m_tot;

	return C;
}

tensor<double> barycentre_v(tensor<double> &p_vec, tensor<double> &m_vector) {
	int D[] = {3};
	tensor<double> C(D,1,0.0);
	double m_tot = 0;
	for(int i=0; i < p_vec.size(0); i++) {
		C(0) += p_vec(i,0);
		C(1) += p_vec(i,1);
		C(2) += p_vec(i,2);
		m_tot += m_vector(i);
	}
	C(0) /= m_tot;
	C(1) /= m_tot;
	C(2) /= m_tot;

	return C;
}

double Total_energy(tensor<double> &q_vec, 
					tensor<double> &p_vec, 
					tensor<double> &m_vector, 
					tensor<int> &type_vec) {
	int D[] = {3};
	tensor<double> C;
	tensor<double> temp_v(D,1,0.0);
	int i,j;
	double E_tot = 0;

	C = barycentre_v(p_vec,m_vector);

	// Kinetic
	for(i=0; i < p_vec.size(0); i++) {
		if(type_vec(i) == 1) {
			temp_v(0) = p_vec(i,0)/m_vector(i) - C(0);
			temp_v(1) = p_vec(i,1)/m_vector(i) - C(1);
			temp_v(2) = p_vec(i,2)/m_vector(i) - C(2);

			E_tot += 0.5*m_vector(i)*dot_product(temp_v,temp_v)/inline_SQR(AU); //why au^-2
		}
	}
	//Potential
	for(i=0; i < (p_vec.size(0)-1); i++) {
		if(type_vec(i) == 1) {
			for(j=i+1; j < p_vec.size(0); j++) {
				if(type_vec(j) == 1) {
					E_tot -= G*m_vector(i)*m_vector(j)/(sqrt(inline_SQR(q_vec(i,0) - q_vec(j,0)) + inline_SQR(q_vec(i,1) - q_vec(j,1)) + inline_SQR(q_vec(i,2) - q_vec(j,2)) ) /inline_SQR(AU)); //why au^-2
				}
			}
		}
	}

	return E_tot;
}

tensor<double> Total_angular_momentum(tensor<double> &q_vec, 
					tensor<double> &p_vec, 
					tensor<double> &m_vector, 
					tensor<int> &type_vec) {
	int D[] = {3};
	tensor<double> Cv,Cx;
	tensor<double> ret(D,1,0.0),tmp(D,1),tmpq(D,1),tmpp(D,1);
	int i,j;
	double E_tot = 0;

	Cx = barycentre(q_vec,m_vector);
	Cv = barycentre_v(p_vec,m_vector);

	for(i=0; i < q_vec.size(0); i++) {
		if(type_vec(i) == 1) {
			tmpq(0) = q_vec(i,0) - Cx(0);
			tmpq(1) = q_vec(i,1) - Cx(1);
			tmpq(2) = q_vec(i,2) - Cx(2);

			tmpp(0) = p_vec(i,0) - Cv(0);
			tmpp(1) = p_vec(i,1) - Cv(1);
			tmpp(2) = p_vec(i,2) - Cv(2);

			tmp = cross_product(tmpq,tmpp);
			ret(0) += tmp(0);
			ret(1) += tmp(1);
			ret(2) += tmp(2);
		}
	}

	return ret;
}

