/*
 * functions.hh
 *
 *  Created on: Jun 17, 2015
 *      Author: dankas
 */


#include "define.hh"
#include "class.hh"

#ifndef DATA_READ_HH
#define DATA_READ_HH

int load_config(OPTIONS *opt,std::string file);
int load_config_vec(OPTIONS *opt,std::string file);
int load_data_vector(std::vector<double> *data,std::string file);
int load_data_matrix(std::vector<std::vector<double> > *data,std::string file);
int fetch_rows_from_matrix_file(std::vector<std::vector<double> > *data, std::string matrix_file, unsigned int n, unsigned int m);

int load_file(const void *data, std::string file, int type);

#endif /* DATA_READ_HH */

#ifndef DATA_WRITE_HH
#define DATA_WRITE_HH

int save_mat(std::string out_file_name, std::vector<std::vector<double> > *M);

int save_batch(std::string file_name, std::string temp_file_name, std::vector<std::vector<double> > *M, unsigned int row_skip);

#endif /* DATA_WRITE_HH */

#ifndef SYSTEM_FUNCTIONS_HH_
#define SYSTEM_FUNCTIONS_HH_

std::vector<std::string> list_files(const std::string& folder);
bool file_exists(const std::string& name);
bool folder_exists(const std::string& name);
std::string get_selfpath();
bool create_file(std::string name);

#endif /* SYSTEM_FUNCTIONS_HH_ */

#ifndef MATH_FUNCTIONS_HH_
#define MATH_FUNCTIONS_HH_

double max_v(std::vector<double> u);
double min_v(std::vector<double> u);

std::vector<std::vector<double> > uint_mat_to_double(std::vector<std::vector<unsigned int> > M);

double abs_double(double a);

int factorial(int n);
int stumpff(int k, int n, double x);

 std::vector<double> rot_z(std::vector<double> u ,double theta);
 std::vector<double> rot_y(std::vector<double> u ,double theta);
 std::vector<double> rot_to_plane(std::vector<double> u ,std::vector<double> n);

std::vector<double> cross_product(std::vector<double> u,std::vector<double> v);
double dot_product(std::vector<double> u,std::vector<double> v);
double abs_v(std::vector<double> u);
std::vector<std::vector<double> > multiply_mat(double a, std::vector<std::vector<double> > M);
std::vector<std::vector<double> > add_v_mat(std::vector<double> a, std::vector<std::vector<double> > M);
std::vector<std::vector<double> > matrix_transpose(const std::vector<std::vector<double> > M);
std::vector<std::vector<double> > swap_cols(const std::vector<std::vector<double> > M, unsigned int a, unsigned int b);

std::vector<double> abs_element_v(std::vector<double> u);
std::vector<double> extract_col(std::vector<std::vector<double> > M, unsigned int id);
std::string IntToString(unsigned int a);

unsigned int draw_from_dist(std::vector<double> dist);
double draw_from_uniform(double a, double b);

std::vector<unsigned int> histogram(const std::vector<double> *data, unsigned int n);
std::vector<double> histogram_bins(const std::vector<double> *data, unsigned int n);

bool is_finite(double x);

#endif /* MATH_FUNCTIONS_HH_ */

#ifndef D_CRITERION_HH_
#define D_CRITERION_HH_

std::vector<std::vector<double> > calculate_metric_matrix(const std::vector<OBSERVATION_record> *list, std::string function);
std::vector<std::vector<double> > augment_metric_matrix(const std::vector<OBSERVATION_record> *list, const std::vector<OBSERVATION_record> *list_addition, std::string function);

double d_SH(const OBSERVATION_record *orb1, const OBSERVATION_record *orb2);
double d_D(const OBSERVATION_record *orb1, const OBSERVATION_record *orb2);
double d_V(const OBSERVATION_record *orb1, const OBSERVATION_record *orb2);
double d_N(const OBSERVATION_record *orb1, const OBSERVATION_record *orb2);
double rho2(const OBSERVATION_record *orb1, const OBSERVATION_record *orb2);
double varrho1(const OBSERVATION_record *orb1, const OBSERVATION_record *orb2);

#endif /* D_CRITERION_HH_ */

#ifndef ASSOCIATION_HH_
#define ASSOCIATION_HH_


int node_associate_file(std::vector<node> *nodes, unsigned int *cluster_counter, double d_crit, unsigned int node_id, int mother_cluster, std::string matrix_folder, std::vector<std::string> file_list, const std::vector<std::vector<unsigned int> > &index_matrix);
std::vector<std::vector<unsigned int> > single_linkage_clustering_files(std::string matrix_folder, std::vector<std::string> file_list, std::vector<std::vector<unsigned int> > index_matrix, double d_crit);

std::vector<std::vector<unsigned int> > single_linkage_clustering_matrix(std::vector<std::vector<double> > *d_matrix, double d_crit);
int node_associate(std::vector<node> *nodes, unsigned int *cluster_counter, double d_crit, unsigned int node_id, int mother_cluster);



#endif /* ASSOCIATION_HH_ */



#ifndef FUNCTIONS_HH_
#define FUNCTIONS_HH_


#endif /* FUNCTIONS_HH_ */


