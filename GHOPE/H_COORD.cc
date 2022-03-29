#include "resources.hh"
#include "functions.hh"

void map2centric(	tensor<double> & q_vec, 
					tensor<double> & p_vec, 
					tensor<double>& m_vector,
					int ID) {
	tensor<double> central_body_q;
	tensor<double> central_body_v;
	central_body_q = q_vec.col(ID);
	central_body_v = p_vec.col(ID)/m_vector(ID);

	for(int i=0; i < q_vec.size(0); i++) {
		q_vec(i,0) -= central_body_q(0);
		q_vec(i,1) -= central_body_q(1);
		q_vec(i,2) -= central_body_q(2);

		p_vec(i,0) -= central_body_v(0)*m_vector(i);
		p_vec(i,1) -= central_body_v(1)*m_vector(i);
		p_vec(i,2) -= central_body_v(2)*m_vector(i);
	}
}

//eccliptic and vernal equinox
void map2ecliptic(	tensor<double> & q_vec, 
					tensor<double> & p_vec, 
					int ID) {
	tensor<double> body_q;
	tensor<double> body_p;
	tensor<double> h;

	body_q = q_vec.col(ID);
	body_p = p_vec.col(ID);

	h = cross_product(body_q,body_p);

	rot_cols_to_plane(q_vec,h);
	rot_cols_to_plane(p_vec,h);
	rot_cols_z(q_vec,-PI*0.5);
	rot_cols_z(p_vec,-PI*0.5);
}
