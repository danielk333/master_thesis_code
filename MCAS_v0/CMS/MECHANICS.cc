

#include <vector>
#include <math.h>

#include "define.hh"
#include "functions.hh"


double JDtoJ(double JD) {
	return 2000.0 + (JD - 2451545.0)/365.25;
}
 
std::vector<double> barycentre(std::vector<std::vector<double> >* q_vec, std::vector<double> *m_vector, std::vector<int>* type_vec) {
	std::vector<double> C;
	C.push_back(0);C.push_back(0);C.push_back(0);
	unsigned int i;
	double m_tot = 0;
	for(i=0; i < (*q_vec).size(); i++) {
		if((*type_vec)[i] == 1) {
			C[0] = C[0] + (*m_vector)[i]*((*q_vec)[i])[0];
			C[1] = C[1] + (*m_vector)[i]*((*q_vec)[i])[1];
			C[2] = C[2] + (*m_vector)[i]*((*q_vec)[i])[2];
			m_tot = m_tot + (*m_vector)[i];
		}
	}
	C[0] = C[0]/m_tot;
	C[1] = C[1]/m_tot;
	C[2] = C[2]/m_tot;

	return C;
}

std::vector<double> barycentre_v(std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec) {
	std::vector<double> C;
	double m_temp = 0;
	C.push_back(0);C.push_back(0);C.push_back(0);
	unsigned int i;
	for(i=0; i < (*p_vec).size(); i++) {
		if((*type_vec)[i] == 1) {
			C[0] = C[0] + ((*p_vec)[i])[0];
			C[1] = C[1] + ((*p_vec)[i])[1];
			C[2] = C[2] + ((*p_vec)[i])[2];
			m_temp = m_temp + (*m_vector)[i];
		}
	}

	C[0] = C[0]/m_temp;
	C[1] = C[1]/m_temp;
	C[2] = C[2]/m_temp;
	return C;
}

std::vector<double> barycentre_p(std::vector<std::vector<double> >* p_vec, std::vector<int>* type_vec) {
	std::vector<double> C;
	double m_temp = 0;
	C.push_back(0);C.push_back(0);C.push_back(0);
	unsigned int i;
	for(i=0; i < (*p_vec).size(); i++) {
		if((*type_vec)[i] == 1) {
			C[0] = C[0] + ((*p_vec)[i])[0];
			C[1] = C[1] + ((*p_vec)[i])[1];
			C[2] = C[2] + ((*p_vec)[i])[2];
		}
	}
	return C;
}

double E_t(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec) {
	std::vector<double> C;
	std::vector<double> temp_v;
	C.push_back(0);C.push_back(0);C.push_back(0);
	unsigned int i,j;
	double m_tot = 0;
	double E_tot = 0;
	for(i=0; i < (*p_vec).size(); i++) {
		if((*type_vec)[i] == 1) {
			C[0] = C[0] + ((*p_vec)[i])[0];
			C[1] = C[1] + ((*p_vec)[i])[1];
			C[2] = C[2] + ((*p_vec)[i])[2];
			m_tot = m_tot + (*m_vector)[i];
		}
	}
	C[0] = C[0]/m_tot;
	C[1] = C[1]/m_tot;
	C[2] = C[2]/m_tot;
	temp_v.push_back(0);
	temp_v.push_back(0);
	temp_v.push_back(0);

	// Kinetic
	for(i=0; i < (*p_vec).size(); i++) {
		if((*type_vec)[i] == 1) {
			temp_v[0] = ((*p_vec)[i])[0]/(*m_vector)[i] - C[0];
			temp_v[1] = ((*p_vec)[i])[1]/(*m_vector)[i] - C[1];
			temp_v[2] = ((*p_vec)[i])[2]/(*m_vector)[i] - C[2];

			E_tot = E_tot + 0.5*(*m_vector)[i]*dot_product(temp_v,temp_v);
		}
	}
	//Potential
	for(i=0; i < ((*p_vec).size()-1); i++) {
		if((*type_vec)[i] == 1) {
			for(j=i+1; j < (*p_vec).size(); j++) {
				if((*type_vec)[i] == 1) {
					E_tot = E_tot - G*(*m_vector)[i]*(*m_vector)[j]/abs_v( (*q_vec)[i] - (*q_vec)[j] );
				}
			}
		}
	}

	return E_tot;
}

std::vector<double> total_angular_momentum_qp_barycenter(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec) {
	std::vector<double> ret;
	ret.push_back(0);
	ret.push_back(0);
	ret.push_back(0);
	std::vector<double> bary,bary_p;
	bary = barycentre(q_vec,m_vector,type_vec);
	bary_p = barycentre_p(p_vec,type_vec);

	int i;
	for(i=(*q_vec).size()-1; i >= 1; i--) {
		if((*type_vec)[i] == 1) {
			ret = ret + cross_product(((*q_vec)[i]) - bary,((*p_vec)[i]) - bary_p);
		}
	}

	return ret;
}

std::vector<double> total_angular_momentum_qp(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<int>* type_vec) {
	std::vector<double> ret;
	ret.push_back(0);
	ret.push_back(0);
	ret.push_back(0);

	int i;
	for(i=(*q_vec).size()-1; i >= 1; i--) {
		if((*type_vec)[i] == 1) {
			ret = ret + cross_product(((*q_vec)[i]),((*p_vec)[i]));
		}
	}

	return ret;
}



std::vector<double> hill_sphere(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec) {
	unsigned int i;
	std::vector<double> kepler_elements;
	std::vector<double> hills;
	for(i=1; i < (*p_vec).size(); i++) {
		if((*type_vec)[i] == 0) {
			kepler_elements = qp_to_kepler((*q_vec)[i],(*p_vec)[i], (*m_vector)[i], (*m_vector)[0]);
			hills.push_back(kepler_elements[0]*(1-kepler_elements[1])*cbrt((*m_vector)[i]/(3*(*m_vector)[0])));
		}
	}
	return hills;

}

void merge(unsigned int ind_A, unsigned int ind_B, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec) {

}



void check_coll(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec,std::vector<std::vector<double> >* q_vec_old, std::vector<std::vector<double> >* p_vec_old, std::vector<double>* m_vector, std::vector<int>* type_vec) {

}

//coll data

std::vector<std::vector<double> > distance_matrix(std::vector<std::vector<double> >* q_vec) {
	std::vector<std::vector<double> > r_MAT;
	std::vector<double> r_MAT_row;
	std::vector<double> temp_q;
	temp_q.push_back(0);temp_q.push_back(0);temp_q.push_back(0);
	unsigned int i,j;

	for(i=0; i < (*q_vec).size()-1; i++) {
		r_MAT_row.clear();
		for(j=i+1; j < (*q_vec).size(); j++) {
			temp_q[0] = (*q_vec)[i][0] - (*q_vec)[j][0];
			temp_q[1] = (*q_vec)[i][1] - (*q_vec)[j][1];
			temp_q[2] = (*q_vec)[i][2] - (*q_vec)[j][2];
			r_MAT_row.push_back(abs_v(temp_q));
		}
		r_MAT.push_back(r_MAT_row);
	}

	return r_MAT;
	
}

