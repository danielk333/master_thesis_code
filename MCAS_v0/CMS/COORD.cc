/*
 * COORD.cc
 *
 *  Created on: Mar 14, 2016
 *      Author: dankas
 */

#include "define.hh"
#include "functions.hh"




void jacobi(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector) {

	std::vector<double> Q_temp;
	std::vector<double> P_temp;
	double m_temp = 0;
	unsigned int i;
	Q_temp.push_back(0);Q_temp.push_back(0);Q_temp.push_back(0);
	P_temp.push_back(0);P_temp.push_back(0);P_temp.push_back(0);

	for(i=1; i < (*q_vec).size(); i++) {
		//Q_temp[0] = 0;Q_temp[1] = 0;Q_temp[2] = 0;
		//P_temp[0] = 0;P_temp[1] = 0;P_temp[2] = 0;
		//for(j=0; j < i; i++) {

			Q_temp[0] = Q_temp[0] + (*m_vector)[i-1]*((*q_vec)[i-1])[0];
			Q_temp[1] = Q_temp[1] + (*m_vector)[i-1]*((*q_vec)[i-1])[1];
			Q_temp[2] = Q_temp[2] + (*m_vector)[i-1]*((*q_vec)[i-1])[2];

			P_temp[0] = P_temp[0] + ((*p_vec)[i-1])[0];
			P_temp[1] = P_temp[1] + ((*p_vec)[i-1])[1];
			P_temp[2] = P_temp[2] + ((*p_vec)[i-1])[2];

			m_temp = m_temp + (*m_vector)[i-1];
		//}

		((*q_vec)[i])[0] = ((*q_vec)[i])[0] - Q_temp[0]/m_temp;
		((*q_vec)[i])[1] = ((*q_vec)[i])[1] - Q_temp[1]/m_temp;
		((*q_vec)[i])[2] = ((*q_vec)[i])[2] - Q_temp[2]/m_temp;

		((*p_vec)[i])[0] = m_temp/(m_temp + (*m_vector)[i])*((*p_vec)[i])[0] - (*m_vector)[i]/(m_temp + (*m_vector)[i])*P_temp[0];
		((*p_vec)[i])[1] = m_temp/(m_temp + (*m_vector)[i])*((*p_vec)[i])[1] - (*m_vector)[i]/(m_temp + (*m_vector)[i])*P_temp[1];
		((*p_vec)[i])[2] = m_temp/(m_temp + (*m_vector)[i])*((*p_vec)[i])[2] - (*m_vector)[i]/(m_temp + (*m_vector)[i])*P_temp[2];
		//(m_temp + (*m_vector)[i])

	}
	Q_temp[0] = Q_temp[0] + (*m_vector)[0]*((*q_vec)[0])[0];
	Q_temp[1] = Q_temp[1] + (*m_vector)[0]*((*q_vec)[0])[1];
	Q_temp[2] = Q_temp[2] + (*m_vector)[0]*((*q_vec)[0])[2];
	P_temp[0] = P_temp[0] + ((*p_vec)[0])[0];
	P_temp[1] = P_temp[1] + ((*p_vec)[0])[1];
	P_temp[2] = P_temp[2] + ((*p_vec)[0])[2];
	m_temp = m_temp + (*m_vector)[0];

	((*q_vec)[0])[0] = Q_temp[0]/m_temp;
	((*q_vec)[0])[1] = Q_temp[1]/m_temp;
	((*q_vec)[0])[2] = Q_temp[2]/m_temp;

	((*p_vec)[0])[0] = P_temp[0];
	((*p_vec)[0])[1] = P_temp[1];
	((*p_vec)[0])[2] = P_temp[2];
	//std::cout << "Q: " << ((*q_vec)[0])[0]/AU << " " << ((*q_vec)[0])[1]/AU << " " << ((*q_vec)[0])[2]/AU << std::endl;
	//std::cout << " : " << ((*q_vec)[1])[0] << " " << ((*q_vec)[1])[1] << " " << ((*q_vec)[1])[2] << std::endl;
	//std::cout << " : " << ((*q_vec)[2])[0] << " " << ((*q_vec)[2])[1] << " " << ((*q_vec)[2])[2] << std::endl;
}

void jacobi_inv(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector) {

	std::vector<double> Q_temp,Q_temp2;
	std::vector<double> P_temp,P_temp2;
	double m_temp = 0;
	unsigned int i;
	Q_temp.push_back(((*q_vec)[0])[0]);Q_temp.push_back(((*q_vec)[0])[1]);Q_temp.push_back(((*q_vec)[0])[2]);
	P_temp.push_back(((*p_vec)[0])[0]);P_temp.push_back(((*p_vec)[0])[1]);P_temp.push_back(((*p_vec)[0])[2]);

	Q_temp2.push_back(0);Q_temp2.push_back(0);Q_temp2.push_back(0);
	P_temp2.push_back(0);P_temp2.push_back(0);P_temp2.push_back(0);

	for(i=1; i < (*q_vec).size(); i++) {
		Q_temp2[0] = Q_temp2[0] + (*m_vector)[i]*((*q_vec)[i])[0];
		Q_temp2[1] = Q_temp2[1] + (*m_vector)[i]*((*q_vec)[i])[1];
		Q_temp2[2] = Q_temp2[2] + (*m_vector)[i]*((*q_vec)[i])[2];
		P_temp2[0] = P_temp2[0] + ((*p_vec)[i])[0];
		P_temp2[1] = P_temp2[1] + ((*p_vec)[i])[1];
		P_temp2[2] = P_temp2[2] + ((*p_vec)[i])[2];
		m_temp = m_temp + (*m_vector)[i];
	}
	m_temp = m_temp + (*m_vector)[0];

	((*q_vec)[0])[0] = ((*q_vec)[0])[0] - Q_temp2[0]/m_temp;
	((*q_vec)[0])[1] = ((*q_vec)[0])[1] - Q_temp2[1]/m_temp;
	((*q_vec)[0])[2] = ((*q_vec)[0])[2] - Q_temp2[2]/m_temp;

	((*p_vec)[0])[0] = (*m_vector)[0]/m_temp*((*p_vec)[0])[0] - P_temp2[0];
	((*p_vec)[0])[1] = (*m_vector)[0]/m_temp*((*p_vec)[0])[0] - P_temp2[1];
	((*p_vec)[0])[2] = (*m_vector)[0]/m_temp*((*p_vec)[0])[0] - P_temp2[2];

	m_temp = 0;
	for(i=1; i < (*q_vec).size(); i++) {
		//Q_temp[0] = 0;Q_temp[1] = 0;Q_temp[2] = 0;
		//P_temp[0] = 0;P_temp[1] = 0;P_temp[2] = 0;
		//for(j=0; j < i; i++) {
		m_temp = m_temp + (*m_vector)[i-1]; //(m_temp + (*m_vector)[i])
		Q_temp[0] = Q_temp[0] + (*m_vector)[i]/(m_temp + (*m_vector)[i])*((*q_vec)[i])[0];
		Q_temp[1] = Q_temp[1] + (*m_vector)[i]/(m_temp + (*m_vector)[i])*((*q_vec)[i])[1];
		Q_temp[2] = Q_temp[2] + (*m_vector)[i]/(m_temp + (*m_vector)[i])*((*q_vec)[i])[2];

		P_temp[0] = (P_temp[0] + ((*p_vec)[i])[0])*(m_temp + (*m_vector)[i])/m_temp;
		P_temp[1] = (P_temp[1] + ((*p_vec)[i])[1])*(m_temp + (*m_vector)[i])/m_temp;
		P_temp[2] = (P_temp[2] + ((*p_vec)[i])[2])*(m_temp + (*m_vector)[i])/m_temp;

		//}

		((*q_vec)[i])[0] = m_temp/(m_temp + (*m_vector)[i])*((*q_vec)[i])[0] - Q_temp[0];
		((*q_vec)[i])[1] = m_temp/(m_temp + (*m_vector)[i])*((*q_vec)[i])[1] - Q_temp[1];
		((*q_vec)[i])[2] = m_temp/(m_temp + (*m_vector)[i])*((*q_vec)[i])[2] - Q_temp[2];

		((*p_vec)[i])[0] = ((*p_vec)[i])[0] + (*m_vector)[i]/(m_temp + (*m_vector)[i])*P_temp[0];
		((*p_vec)[i])[1] = ((*p_vec)[i])[1] + (*m_vector)[i]/(m_temp + (*m_vector)[i])*P_temp[1];
		((*p_vec)[i])[2] = ((*p_vec)[i])[2] + (*m_vector)[i]/(m_temp + (*m_vector)[i])*P_temp[2];
		//(m_temp + (*m_vector)[i])

	}

}
void canonical_heliocentric(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector) {

	std::vector<double> Q_temp;
	std::vector<double> P_temp;
	double m_temp = 0;
	unsigned int i;
	Q_temp.push_back(0);Q_temp.push_back(0);Q_temp.push_back(0);
	P_temp.push_back(0);P_temp.push_back(0);P_temp.push_back(0);

	for(i=0; i < (*q_vec).size(); i++) {
		Q_temp[0] = Q_temp[0] + (*m_vector)[i]*((*q_vec)[i])[0];
		Q_temp[1] = Q_temp[1] + (*m_vector)[i]*((*q_vec)[i])[1];
		Q_temp[2] = Q_temp[2] + (*m_vector)[i]*((*q_vec)[i])[2];

		P_temp[0] = P_temp[0] + ((*p_vec)[i])[0];
		P_temp[1] = P_temp[1] + ((*p_vec)[i])[1];
		P_temp[2] = P_temp[2] + ((*p_vec)[i])[2];

		m_temp = m_temp + (*m_vector)[i];
	}

	for(i=1; i < (*q_vec).size(); i++) {
		((*q_vec)[i])[0] = ((*q_vec)[i])[0] - ((*q_vec)[0])[0];
		((*q_vec)[i])[1] = ((*q_vec)[i])[1] - ((*q_vec)[0])[1];
		((*q_vec)[i])[2] = ((*q_vec)[i])[2] - ((*q_vec)[0])[2];

		((*p_vec)[i])[0] = ((*p_vec)[i])[0] - (*m_vector)[i]/m_temp*P_temp[0];
		((*p_vec)[i])[1] = ((*p_vec)[i])[1] - (*m_vector)[i]/m_temp*P_temp[1];
		((*p_vec)[i])[2] = ((*p_vec)[i])[2] - (*m_vector)[i]/m_temp*P_temp[2];
	}

	((*q_vec)[0])[0] = Q_temp[0]/m_temp;
	((*q_vec)[0])[1] = Q_temp[1]/m_temp;
	((*q_vec)[0])[2] = Q_temp[2]/m_temp;

	((*p_vec)[0])[0] = P_temp[0];
	((*p_vec)[0])[1] = P_temp[1];
	((*p_vec)[0])[2] = P_temp[2];

}

void canonical_heliocentric_inv(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector) {

	std::vector<double> Q_temp;
	std::vector<double> P_temp;
	double m_temp;
	unsigned int i;
	Q_temp.push_back(0);Q_temp.push_back(0);Q_temp.push_back(0);
	P_temp.push_back(0);P_temp.push_back(0);P_temp.push_back(0);

	m_temp = (*m_vector)[0];
	for(i=1; i < (*q_vec).size(); i++) {
		Q_temp[0] = Q_temp[0] + (*m_vector)[i]*((*q_vec)[i])[0];
		Q_temp[1] = Q_temp[1] + (*m_vector)[i]*((*q_vec)[i])[1];
		Q_temp[2] = Q_temp[2] + (*m_vector)[i]*((*q_vec)[i])[2];

		P_temp[0] = P_temp[0] + ((*p_vec)[i])[0];
		P_temp[1] = P_temp[1] + ((*p_vec)[i])[1];
		P_temp[2] = P_temp[2] + ((*p_vec)[i])[2];

		m_temp = m_temp + (*m_vector)[i];
	}

	((*q_vec)[0])[0] = ((*q_vec)[0])[0] - Q_temp[0]/m_temp;
	((*q_vec)[0])[1] = ((*q_vec)[0])[1] - Q_temp[1]/m_temp;
	((*q_vec)[0])[2] = ((*q_vec)[0])[2] - Q_temp[2]/m_temp;

	for(i=1; i < (*q_vec).size(); i++) {
		((*q_vec)[i])[0] = ((*q_vec)[i])[0] + ((*q_vec)[0])[0];
		((*q_vec)[i])[1] = ((*q_vec)[i])[1] + ((*q_vec)[0])[1];
		((*q_vec)[i])[2] = ((*q_vec)[i])[2] + ((*q_vec)[0])[2];

		((*p_vec)[i])[0] = ((*p_vec)[i])[0] + (*m_vector)[i]/m_temp*((*p_vec)[0])[0];
		((*p_vec)[i])[1] = ((*p_vec)[i])[1] + (*m_vector)[i]/m_temp*((*p_vec)[0])[1];
		((*p_vec)[i])[2] = ((*p_vec)[i])[2] + (*m_vector)[i]/m_temp*((*p_vec)[0])[2];
	}

	((*p_vec)[0])[0] = ((*m_vector)[0]/m_temp)*((*p_vec)[0])[0] - P_temp[0];
	((*p_vec)[0])[1] = ((*m_vector)[0]/m_temp)*((*p_vec)[0])[1] - P_temp[1];
	((*p_vec)[0])[2] = ((*m_vector)[0]/m_temp)*((*p_vec)[0])[2] - P_temp[2];
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

int rot_all_y(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, double theta) {
	std::cout << "Rotating all objects around y axis by " << theta*(180/PI) << " degrees " << std::endl;
	int i;
	for(i=(*q_vec).size()-1; i >= 0; i--) {
		((*q_vec)[i]) = rot_y(((*q_vec)[i]),theta);
		((*p_vec)[i]) = rot_y(((*p_vec)[i]),theta);
	}
	std::cout << "Rotation complete" << std::endl;

	return SIM_OK;
}

int rot_all_z(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, double theta) {
	std::cout << "Rotating all objects around z axis by " << theta*(180/PI) << " degrees " << std::endl;
	int i;
	for(i=(*q_vec).size()-1; i >= 1; i--) {
		((*q_vec)[i]) = rot_z(((*q_vec)[i]),theta);
		((*p_vec)[i]) = rot_z(((*p_vec)[i]),theta);
	}
	std::cout << "Rotation complete" << std::endl;

	return SIM_OK;
}

int rot_to_invariable_plane(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<int>* type_vec) {
	std::cout << "Converting coordinate system to the invariable plane." << std::endl;
	std::vector<double> L;
	L = total_angular_momentum_qp(q_vec, p_vec, type_vec);

	double L_t = abs_v(L);
	double phi = acos(L[2]/L_t);
	double theta = atan2(L[1],L[0]);
	std::cout << "Rotation around invertable plan of theta=" << theta*(180/PI) << ", phi=" << phi*(180/PI) << " detected." << std::endl;
	std::cout << "Total angular momentum: " << L_t << std::endl;

	std::cout << "Anti rotation result: ";

	int i;
	for(i=(*q_vec).size()-1; i >= 1; i--) {
		((*q_vec)[i]) = rot_to_plane(((*q_vec)[i]),L);
		((*p_vec)[i]) = rot_to_plane(((*p_vec)[i]),L);
	}
	std::cout << "Conversion done." << std::endl << std::endl;

	return SIM_OK;
}