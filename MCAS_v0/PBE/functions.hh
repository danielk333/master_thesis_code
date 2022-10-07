/*
 * functions.hh
 *
 *  Created on: May 6, 2015
 *      Author: dankas
 */


#include "define.hh"

  #include "ast2body.h"

#ifndef DATA_READ_HH
#define DATA_READ_HH

int load_config(std::string file, OPTIONS *opt);
int load_data_vector(std::vector<double> *data,std::string file);
int load_data_matrix(std::vector<std::vector<double> > *data,std::string file);
int load_previus_settings(char *in, OPTIONS *opt, SIMULATION *sim);

#endif /* DATA_READ_HH */

#ifndef DATA_WRITE_HH
#define DATA_WRITE_HH

int save_mat(const char *out, std::vector<std::vector<double> > *M, OPTIONS *opt);

#endif /* DATA_WRITE_HH */

#ifndef MATH_FUNCTIONS_HH_
#define MATH_FUNCTIONS_HH_

double neville(std::vector<double>* x,std::vector<double>* y,double x0);
std::vector<std::vector<double> > fill_mat(unsigned int a, unsigned int b, double X);

double abs_working(double x);
bool is_finite(double x);
std::vector<double> cross_product(std::vector<double> u,std::vector<double> v);
double dot_product(std::vector<double> u,std::vector<double> v);
double abs_v(std::vector<double> u);
std::vector<std::vector<double> > generate_random_sphere_normals(unsigned int N);
std::vector<std::vector<double> > multiply_mat(double a, std::vector<std::vector<double> > M);
std::vector<std::vector<double> > add_v_mat(std::vector<double> a, std::vector<std::vector<double> > M);
double sum_v(std::vector<double> u);
std::vector<double> multiply_v(double a, std::vector<double> u);
std::vector<double> add_v(double a, std::vector<double> u);
unsigned int draw_from_dist(std::vector<double> dist);
std::vector<double> add_v_v(std::vector<double> a, std::vector<double> M);
std::vector<std::vector<double> > matrix_transpose(const std::vector<std::vector<double> > M);

#endif /* MATH_FUNCTIONS_HH_ */

#ifndef SYSTEM_FUNCTIONS_HH_
#define SYSTEM_FUNCTIONS_HH_

bool file_exists(const std::string& name);
int console_output(SIMULATION *sim, OPTIONS *opt);
std::string get_selfpath();

#endif /* SYSTEM_FUNCTIONS_HH_ */

#ifndef MECHANICS_HH_
#define MECHANICS_HH_

std::vector<double> total_angular_momentum_qv_heliocentric(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* v_vec, std::vector<double>* m_vec);
std::vector<double> total_angular_momentum_qp_heliocentric(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<int>* type_vec);
double E_t(std::vector<std::vector<double> > q_vec, std::vector<std::vector<double> > p_vec, std::vector<double> m_vector, std::vector<int> type_vec);
std::vector<double> barycentre(std::vector<std::vector<double> > q_vec, std::vector<double> m_vector);
std::vector<double> barycentre_v(std::vector<std::vector<double> > p_vec, std::vector<double> m_vector);

#endif /* MECHANICS_HH_ */

#ifndef INTEGRATOR_SETTINGS_HH_
#define INTEGRATOR_SETTINGS_HH_

int load_integrator_scheme(std::vector<std::vector<double> > *integrator_coef, OPTIONS *opt);

#endif /* INTEGRATOR_SETTINGS_HH_ */

#ifndef H_INTEGRATORS_HH_
#define H_INTEGRATORS_HH_

void restore_coord_qp(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vec, std::vector<double> coord);
std::vector<double> temp_heliocentric_qp(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vec);

void heliocentric_q(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vec);
void heliocentric_qp(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vec);
void heliocentric_qv(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* v_vec);

void H_kin(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec, double dt);
void H_pot(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec, double dt, unsigned int include_center);

void H_pot_tp(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec,double dt, unsigned int include_center);
void H_pot_tp_DIS(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec,double dt, unsigned int include_center,double rho);

void BS_step(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vector, std::vector<int>* type_vec, std::vector<unsigned int> N, double dt, double rho);
std::vector<std::vector<double> > dZ(std::vector<std::vector<double> >* Z, std::vector<double>* m_vector, std::vector<int>* type_vec, unsigned int ind, double rho);

#endif /* H_INTEGRATORS_HH_ */

#ifndef KEPLER_HH_
#define KEPLER_HH_

/*std::vector<double> qp_to_kepler(std::vector<double> q, std::vector<double> p, double m, double M);
std::vector<double> qp_to_kepler_sun(std::vector<double> qs, std::vector<double> ps, std::vector<double> q, std::vector<double> p, double m, double M);
std::vector<double> xv_to_kepler(std::vector<double> x, std::vector<double> v, double mu);
std::vector<double> kepler_to_qp(std::vector<double> kepler, double m, double M);*/

std::vector<double> qp_to_kepler(std::vector<double> q, std::vector<double> p, double m, double M);
std::vector<double> kepler_to_qp(std::vector<double> kepler, double m, double M);
std::vector<double> xv_to_kepler(std::vector<double> x, std::vector<double> v, double mu);
std::vector<double> kepler_to_xv(std::vector<double> kepler, double mu);

#endif /* KEPLER_HH_ */


#ifndef FUNCTIONS_HH_
#define FUNCTIONS_HH_

#endif /* FUNCTIONS_HH_ */
