#include "resources.hh"

//FUNCTIONS RELATING TO DATA READ
//
//*****************************************************************************************
#ifndef DATA_READ_HH
#define DATA_READ_HH

void load_config(std::string file, OPTIONS *opt);
void load_data_matrix(std::vector<std::vector<double> > *data,std::string file, std::string delim);

#endif /* DATA_READ_HH */



//FUNCTIONS RELATING TO DATA WRITE
//
//*****************************************************************************************
#ifndef DATA_WRITE_HH
#define DATA_WRITE_HH

void save_mat(std::string out_file_name, std::vector<std::vector<double> > *M, unsigned int precision, std::string delim);

void save_tensor(std::string out_file_name, tensor<double> *M, unsigned int precision, std::string delim);

void merge_data(tensor<double> *collect_tens,
                          int save_cnt,
                          tensor<double> *q_vec, 
                          tensor<double> *p_vec, 
                          tensor<double> *m_vector, 
                          tensor<int> *type_vector, 
                          tensor<int> *body_id_vector,
                          double T);

#endif /* DATA_WRITE_HH */


//FUNCTIONS RELATING TO DATA TYPE
// These functions are declared in the DATA_TYPE.hh file and are only given 
// here as comments for compleatness of the functions header
//
// SEE ALSO OVERLOADING.hh 
// for a complete record of all overloading of operators (+,-,*,/)
//*****************************************************************************************
/*

template <class T>
void tensor<T>::print(std::ostream& out)

template <class T>
void vector2tensor(tensor<T> &X,const std::vector<T> &Y);
template <class T>
void vvector2tensor(tensor<T> &X,const std::vector<std::vector<T> > &Y);
template <class T>
void tensor2vector(std::vector<T> &X, const tensor<T> &Y);
template <class T>
void tensor2vvector(std::vector<std::vector<T> > &X,const tensor<T> &Y);
template <class T>
tensor<T> sort_tensor_ret(const tensor<T> X);
template <class T>
void sort_tensor(const tensor<T> &X);
*/


//FUNCTIONS RELATING TO MATH
//
//*****************************************************************************************
#ifndef MATH_FUNCTIONS_HH_
#define MATH_FUNCTIONS_HH_

//Vector manipulation/function
tensor<double> cross_product(tensor<double> &u, tensor<double> &v);
tensor<double> cross_product(const tensor<double> &u, const tensor<double> &v);
double dot_product(tensor<double> &u, tensor<double> &v);
double dot_product(const tensor<double> &u, const tensor<double> &v);
double norm(const tensor<double> &u);
double abs_v(tensor<double> &u);
double abs_v(const tensor<double> &u);
inline void rot_z(tensor<double> &u ,double theta);
inline void rot_y(tensor<double> &u ,double theta);
void rot_cols_z(tensor<double> &u ,double theta);
void rot_cols_y(tensor<double> &u ,double theta);
void rot_rows_z(tensor<double> &u ,double theta);
void rot_rows_y(tensor<double> &u ,double theta);
void rot_to_plane(tensor<double> &u ,const tensor<double> &n);
void rot_cols_to_plane(tensor<double> &u ,const tensor<double> &n);
void rot_rows_to_plane(tensor<double> &u ,const tensor<double> &n);

//Number checks/manipulation
inline bool is_finite(double x);
inline double abs2(double x);

//Basic math functions
inline int factorial(int n);

#endif /* MATH_FUNCTIONS_HH_ */

//Inline functions for speed

template<class T>
inline T inline_SQR(const T a) {return a*a;}
template<class T>
inline T inline_CUB(const T a) {return a*a*a;}

template<class T>
inline const T &inline_MAX(const T &a, const T &b)
        {return b > a ? (b) : (a);}

template<class T>
inline const T &inline_MIN(const T &a, const T &b)
        {return b < a ? (b) : (a);}

//FUNCTIONS RELATING TO DATES AND TIME
//
//*****************************************************************************************
#ifndef DATETIME_HH_
#define DATETIME_HH_

double JDtoJ(double JD);

#endif /* DATETIME_HH_ */

//FUNCTIONS RELATING TO MECHANICS (Physics)
//
//*****************************************************************************************
#ifndef MECHANICS_HH_
#define MECHANICS_HH_

//CLASSIC MECHANICS
tensor<double> barycentre(tensor<double> &q_vec, tensor<double> &m_vector);
tensor<double> barycentre_v(tensor<double> &p_vec, tensor<double> &m_vector);
double Total_energy(tensor<double> &q_vec, 
					tensor<double> &p_vec, 
					tensor<double> &m_vector, 
					tensor<int> &type_vec);
tensor<double> Total_angular_momentum(tensor<double> &q_vec, 
					tensor<double> &p_vec, 
					tensor<double> &m_vector, 
					tensor<int> &type_vec);

//Celestial mechnics 

//Helpers

#endif /* MECHANICS_HH_ */



//FUNCTIONS RELATING TO SYSTEM INFRASTRUCTURE
//
//*****************************************************************************************
#ifndef SYSTEM_FUNCTIONS_HH_
#define SYSTEM_FUNCTIONS_HH_

std::string exec(const char* cmd);
bool file_exists(const std::string& name);
std::string get_selfpath();

#endif /* SYSTEM_FUNCTIONS_HH_ */


//FUNCTIONS RELATING TO HAMILTONIAN SPLITS
//
//*****************************************************************************************
#ifndef H_SPLITS_HH_
#define H_SPLITS_HH_
#endif /* H_SPLITS_HH_ */

void H_kin(	tensor<double> & q_vec, 
			tensor<double> & p_vec, 
			tensor<double>& m_vector, 
			tensor<int> & type_vec, 
			double dt);
void H_pot(	tensor<double> & q_vec, 
			tensor<double> & p_vec, 
			tensor<double>& m_vector, 
			tensor<int> & type_vec, 
			double dt, 
			tensor<double>& interaction_data);


//FUNCTIONS RELATING TO COORDINATE TRANSFORMATIONS
//
//*****************************************************************************************
#ifndef H_COORD_HH_
#define H_COORD_HH_

void map2centric(	tensor<double> & q_vec, 
					tensor<double> & p_vec, 
					tensor<double>& m_vector,
					int ID);

void map2ecliptic(	tensor<double> & q_vec, 
					tensor<double> & p_vec, 
					int ID);


#endif /* H_COORD_HH_ */


//FUNCTIONS RELATING TO ODE INTEGRATIONS
//
//*****************************************************************************************
#ifndef ODE_INT_HH_
#define ODE_INT_HH_

//COMPLETE THIS WITH TENSOR

void BS_step(	std::vector<std::vector<double> >* q_vec, 
				std::vector<std::vector<double> >* p_vec, 
				std::vector<double>* m_vector, 
				std::vector<int>* type_vec, 
				std::vector<unsigned int> N, 
				double dt);

std::vector<std::vector<double> > dZ(	std::vector<std::vector<double> >* Z, 
										std::vector<double>* m_vector, 
										std::vector<int>* type_vec, 
										unsigned int ind);


#endif /* ODE_INT_HH_ */

