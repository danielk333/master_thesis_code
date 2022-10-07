/*
 * functions.hh
 *
 *  Created on: May 6, 2015
 *      Author: dankas
 */
#include "define.hh"

#ifndef DATA_READ_HH
#define DATA_READ_HH

int load_config(std::string file, OPTIONS *opt);
int load_data_vector(std::vector<double> *data,std::string file);
int load_data_matrix(std::vector<std::vector<double> > *data,std::string file);
int load_previus_settings(char *in, OPTIONS *opt, SIMULATION *sim);

#endif /* DATA_READ_HH */

#ifndef DATA_WRITE_HH
#define DATA_WRITE_HH

int save_mat(std::string out_file_name, std::vector<std::vector<double> > *M);
int save(const char *out, std::vector<std::vector<std::vector<double> > > *q, std::vector<std::vector<std::vector<double> > > *p, std::vector<std::vector<double> > *m_vector, std::vector<std::vector<int> > *type_vec, std::vector<double>* time_vec, OPTIONS *opt, SIMULATION *sim);
//int save(const char *out, std::vector<std::vector<std::vector<double> > > *q, std::vector<std::vector<std::vector<double> > > *p, std::vector<double> *m_vector, std::vector<int> *type_vec, OPTIONS *opt, SIMULATION *sim);
int save_diagnostics(const char *out, std::vector<double> E_vec, std::vector<double> T_vec);
int save_settings(const char *out, OPTIONS *opt, SIMULATION *sim);
int add_save_data_to_memory(OPTIONS *opt, std::vector<std::vector<std::vector<double> > > *q_vec_save, std::vector<std::vector<std::vector<double> > > *p_vec_save, std::vector<std::vector<double> > *q_vec, std::vector<std::vector<double> > *p_vec, std::vector<double> *m_vector);

#endif /* DATA_WRITE_HH */

#ifndef MATH_FUNCTIONS_HH_
#define MATH_FUNCTIONS_HH_

// SORTING
std::vector<size_t> sort_indexes_2(std::vector<double> x);

//Function fitting
double neville(std::vector<double>* x,std::vector<double>* y,double x0);

//Matrix manipulation/functions
std::vector<std::vector<double> > fill_mat(unsigned int a, unsigned int b, double X);

//Vector manipulation/function
std::vector<double> cross_product(std::vector<double> u,std::vector<double> v);
double dot_product(std::vector<double> u,std::vector<double> v);
double abs_v(std::vector<double> u);
std::vector<double> rot_z(std::vector<double> u ,double theta);
std::vector<double> rot_y(std::vector<double> u ,double theta);
std::vector<double> rot_to_plane(std::vector<double> u ,std::vector<double> n);

//Number checks/manipulation
bool is_finite(double x);
double abs_working(double x);

//Basic math functions
int factorial(int n);
int stumpff(int k, int n, double x);

#endif /* MATH_FUNCTIONS_HH_ */


#ifndef MECHANICS_HH_
#define MECHANICS_HH_

//TIME
double JDtoJ(double JD);

//CLASSIC MECHANICS
std::vector<double> barycentre(std::vector<std::vector<double> >* q_vec, std::vector<double> *m_vector, std::vector<int>* type_vec);
std::vector<double> barycentre_v(std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec);
std::vector<double> barycentre_p(std::vector<std::vector<double> >* p_vec, std::vector<int>* type_vec);
double E_t(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec);
std::vector<double> total_angular_momentum_qp_barycenter(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec);
std::vector<double> total_angular_momentum_qp(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<int>* type_vec);

//Celestial mechnics 
std::vector<double> hill_sphere(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec);
void merge(unsigned int ind_A, unsigned int ind_B, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec);
void check_coll(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec,std::vector<std::vector<double> >* q_vec_old, std::vector<std::vector<double> >* p_vec_old, std::vector<double>* m_vector, std::vector<int>* type_vec);

//Helpers
std::vector<std::vector<double> > distance_matrix(std::vector<std::vector<double> >* q_vec);

#endif /* MECHANICS_HH_ */


#ifndef SYSTEM_FUNCTIONS_HH_
#define SYSTEM_FUNCTIONS_HH_

bool file_exists(const std::string& name);
int calculate_memory_management(OPTIONS *opt);
int console_output(SIMULATION *sim, OPTIONS *opt);
std::string get_selfpath();

#endif /* SYSTEM_FUNCTIONS_HH_ */

#ifndef INTEGRATOR_SETTINGS_HH_
#define INTEGRATOR_SETTINGS_HH_

int load_integrator_scheme(std::vector<std::vector<double> > *integrator_coef, OPTIONS *opt);

#endif /* INTEGRATOR_SETTINGS_HH_ */


#ifndef H_INTEGRATORS_HH_
#define H_INTEGRATORS_HH_

void BS_step(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec, std::vector<unsigned int> N, double dt);
std::vector<std::vector<double> > dZ(std::vector<std::vector<double> >* Z, std::vector<double>* m_vector, std::vector<int>* type_vec, unsigned int ind, double rho);

/*
void H_sun_DH(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec,double dt);
void H_int_DH(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec,double dt);
void H_int_hills_DH(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec, double dt);
*/
void H_kep(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec, double dt, double sun_collision, double ejection_criteria);
void H_dri(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec, double dt);
void H_kin(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec, double dt);
void H_pot(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec, double dt, unsigned int include_center);
void H_pot_tp(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec,double dt, unsigned int include_center);
void H_pot_tp_DIS(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec,double dt, unsigned int include_center);
std::vector<std::vector<double> > H_kep_qp(std::vector<double> q_vec, std::vector<double> p_vec, double m, double M, double dt);
std::vector<std::vector<double> > H_kep_xv(std::vector<double> q_vec, std::vector<double> v_vec, double M, double dt);

/*
void H_int_leapfrog(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec,double dt);
void H_int_hills_leapfrog(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec, double dt);
*/
/*
void H_mod_mid(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec, double dt, double n, double sun_collision, double ejection_criteria);
void H_euler(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec, double dt, double sun_collision, double ejection_criteria);
*/
#endif /* H_INTEGRATORS_HH_ */

#ifndef H_COORD_HH_
#define H_COORD_HH_

void jacobi(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector);
void jacobi_inv(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector);

void canonical_heliocentric(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector);
void canonical_heliocentric_inv(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector);

void heliocentric_qp(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vec);
void heliocentric_qv(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* v_vec);

int rot_to_invariable_plane(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<int>* type_vec);
int rot_all_y(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* v_vec, double theta);
int rot_all_z(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* v_vec, double theta);

#endif /* H_COORD_HH_ */

#ifndef KEPLER_HH_
#define KEPLER_HH_


double universal_kepler_equation_newton(double dt, double tol, double a, double r_0, double v_r0, double mu);
double kepler_equation_newton(double M, double e, double E_0, double tol);

double nu_to_E(double e,double nu);

std::vector<double> qp_to_kepler(std::vector<double> q, std::vector<double> p, double m, double M);
std::vector<double> kepler_to_qp(std::vector<double> kepler, double m, double M);
std::vector<double> xv_to_kepler(std::vector<double> x, std::vector<double> v, double mu);
std::vector<double> kepler_to_xv(std::vector<double> kepler, double mu);

#endif /* KEPLER_HH_ */

#ifndef FUNCTIONS_HH_
#define FUNCTIONS_HH_

#endif /* FUNCTIONS_HH_ */
