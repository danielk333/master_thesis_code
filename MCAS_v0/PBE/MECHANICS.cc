

#include <vector>
#include <math.h>

#include "define.hh"
#include "functions.hh"

 
std::vector<double> barycentre(std::vector<std::vector<double> > q_vec, std::vector<double> m_vector) {
	std::vector<double> C;
	C.push_back(0);C.push_back(0);C.push_back(0);
	unsigned int i;
	double m_tot = 0;
	for(i=0; i < q_vec.size(); i++) {
		C[0] = C[0] + m_vector[i]*(q_vec[i])[0];
		C[1] = C[1] + m_vector[i]*(q_vec[i])[1];
		C[2] = C[2] + m_vector[i]*(q_vec[i])[2];
		m_tot = m_tot + m_vector[i];
	}
	C[0] = C[0]/m_tot;
	C[1] = C[1]/m_tot;
	C[2] = C[2]/m_tot;

	return C;
}

std::vector<double> barycentre_v(std::vector<std::vector<double> > p_vec, std::vector<double> m_vector) {
	std::vector<double> C;
	double m_temp = 0;
	C.push_back(0);C.push_back(0);C.push_back(0);
	unsigned int i;
	for(i=0; i < p_vec.size(); i++) {
		C[0] = C[0] + (p_vec[i])[0];
		C[1] = C[1] + (p_vec[i])[1];
		C[2] = C[2] + (p_vec[i])[2];
		m_temp = m_temp + m_vector[i];
	}

	C[0] = C[0]/m_temp;
	C[1] = C[1]/m_temp;
	C[2] = C[2]/m_temp;
	return C;
}

double E_t(std::vector<std::vector<double> > q_vec, std::vector<std::vector<double> > p_vec, std::vector<double> m_vector, std::vector<int> type_vec) {
	std::vector<double> C;
	std::vector<double> temp_v;
	C.push_back(0);C.push_back(0);C.push_back(0);
	unsigned int i,j;
	double m_tot = 0;
	double E_tot = 0;
	for(i=0; i < p_vec.size(); i++) {
		if(type_vec[i] == 1) {
			C[0] = C[0] + (p_vec[i])[0];
			C[1] = C[1] + (p_vec[i])[1];
			C[2] = C[2] + (p_vec[i])[2];
			m_tot = m_tot + m_vector[i];
		}
	}
	C[0] = C[0]/m_tot;
	C[1] = C[1]/m_tot;
	C[2] = C[2]/m_tot;
	temp_v.push_back(0);
	temp_v.push_back(0);
	temp_v.push_back(0);

	// Kinetic
	for(i=0; i < p_vec.size(); i++) {
		if(type_vec[i] == 1) {
			temp_v[0] = (p_vec[i])[0]/m_vector[i] - C[0];
			temp_v[1] = (p_vec[i])[1]/m_vector[i] - C[1];
			temp_v[2] = (p_vec[i])[2]/m_vector[i] - C[2];

			E_tot = E_tot + 0.5*m_vector[i]*dot_product(temp_v,temp_v)/pow(AU,2);
		}
	}
	//Potential
	for(i=0; i < (p_vec.size()-1); i++) {
		if(type_vec[i] == 1) {
			for(j=i+1; j < p_vec.size(); j++) {
				if(type_vec[j] == 1) {
					E_tot = E_tot - Ggrav*m_vector[i]*m_vector[j]/sqrt(pow((q_vec[i])[0] - (q_vec[j])[0],2) + pow((q_vec[i])[1] - (q_vec[j])[1],2) + pow((q_vec[i])[2] - (q_vec[j])[2],2))/pow(AU,2);
				}
			}
		}
	}

	return E_tot;
}

std::vector<double> total_angular_momentum_qp_heliocentric(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<int>* type_vec) {
	std::vector<double> ret;
	ret.push_back(0);
	ret.push_back(0);
	ret.push_back(0);

	int i;
	for(i=(*q_vec).size()-1; i >= 1; i--) {
		if((*type_vec)[i] == 1) {
			ret = add_v_v(ret,cross_product(((*q_vec)[i]),((*p_vec)[i])));
		}
	}

	return ret;
}

