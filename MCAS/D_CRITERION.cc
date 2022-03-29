#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "define.hh"
#include "functions.hh"

/*


float d_N(data_t *data, unsigned int a_ind, unsigned int b_ind) {
   float lambda_a = (*data).at(a_ind).at(36)*PI/180;
   float lambda_b = (*data).at(b_ind).at(36)*PI/180;
   float v_g_a = (*data).at(a_ind).at(17);
   float v_g_b = (*data).at(b_ind).at(17);
   float Ra_a = (*data).at(a_ind).at(9)*PI/180;
   float Ra_b = (*data).at(b_ind).at(9)*PI/180;
   float Dec_a = (*data).at(a_ind).at(10)*PI/180;
   float Dec_b = (*data).at(b_ind).at(10)*PI/180;
   float epsilon = 23.5*PI/180;

   float U_a_1 = v_g_a/29.7*(-cos(lambda_a)*cos(Ra_a)*cos(Dec_a) + sin(lambda_a)*cos(epsilon)*cos(Dec_a)*sin(Ra_a)+sin(lambda_a)*sin(epsilon)*sin(Dec_a));
   float U_a_2 = v_g_a/29.7*(-sin(lambda_a)*cos(Ra_a)*cos(Dec_a) - cos(lambda_a)*cos(epsilon)*cos(Dec_a)*sin(Ra_a) - cos(lambda_a)*sin(epsilon)*sin(Dec_a));
   float U_a_3 = v_g_a/29.7*(sin(epsilon)*cos(Dec_a)*sin(Ra_a) - cos(epsilon)*sin(Dec_a));

   float U_b_1 = v_g_b/29.7*(-cos(lambda_b)*cos(Ra_b)*cos(Dec_b) + sin(lambda_b)*cos(epsilon)*cos(Dec_b)*sin(Ra_b)+sin(lambda_b)*sin(epsilon)*sin(Dec_b));
   float U_b_2 = v_g_b/29.7*(-sin(lambda_b)*cos(Ra_b)*cos(Dec_b) - cos(lambda_b)*cos(epsilon)*cos(Dec_b)*sin(Ra_b) - cos(lambda_b)*sin(epsilon)*sin(Dec_b));
   float U_b_3 = v_g_b/29.7*(sin(epsilon)*cos(Dec_b)*sin(Ra_b) - cos(epsilon)*sin(Dec_b));

   float phi_a = atan(U_a_1/U_a_3);
   float phi_b = atan(U_b_1/U_b_3);

   float U_a_norm = sqrt(pow(U_a_1,2)+pow(U_a_2,2)+pow(U_a_3,2));
   float U_b_norm = sqrt(pow(U_b_1,2)+pow(U_b_2,2)+pow(U_b_3,2));

   float c_theta_a = U_a_2/(U_a_norm);
   float c_theta_b = U_b_2/(U_b_norm);

   float w_1 = 1;
   float w_2 = 1;
   float w_3 = 1;

   float delta_phi_1 = 2*sin((phi_b - phi_a)/2);
   float delta_phi_2 = 2*sin((PI + phi_b - phi_a)/2);

   float delta_lambda_1 = 2*sin((lambda_b - lambda_a)/2);
   float delta_lambda_2 = 2*sin((PI + lambda_b - lambda_a)/2);

   float delta_xi2 = fmin(w_2*pow(delta_phi_1,2) + w_3*pow(delta_lambda_1,2),w_2*pow(delta_phi_2,2) + w_3*pow(delta_lambda_2,2));

   float d = pow(U_b_norm - U_a_norm,2) + w_1*pow(c_theta_b - c_theta_a,2) + delta_xi2;

   return d;
}
*/

std::vector<std::vector<double> > calculate_metric_matrix(const std::vector<OBSERVATION_record> *list, std::string function) {
   unsigned int i,j;
   std::vector<double> v;
   std::vector<std::vector<double> > M;
   double temp_val;
   std::cout << "-- Metric matrix " << function << " row " << 0 << " of " << ((*list).size()-1) << " done." << std::flush;
   if(function.compare("SH") == 0) {
      for(i=0; i < ((*list).size()-1); i++) {
         v.clear();
         std::cout << "\r-- Metric matrix " << function << " row " << i << " of " << ((*list).size()-1) << " done." << std::flush;
         for(j=i+1; j < (*list).size(); j++) {
            temp_val = d_SH( &((*list)[i]), &((*list)[j]));
            if( is_finite(temp_val) ) {
               v.push_back( temp_val );
            }
            else {
               std::cout << std::endl;
               std::cout << "Calculated D value " << temp_val << " is not finite, please examine data and function." << std::endl;
               std::cout << "Object " << i << std::endl;

   /*std::cout << "-- Orbital elements: " << std::endl;
   std::cout << "a=" <<(*list)[i].a << std::endl;
   std::cout << "e=" <<(*list)[i].e << std::endl;
   std::cout << "i=" <<(*list)[i].i << std::endl;
   std::cout << "omega=" <<(*list)[i].omega << std::endl;
   std::cout << "Omega=" <<(*list)[i].Omega << std::endl;
   std::cout << "-- Observational elements: " << std::endl;
   std::cout << "ra=" <<(*list)[i].ra << std::endl;
   std::cout << "dec=" <<(*list)[i].dec << std::endl;
   std::cout << "v_g=" <<(*list)[i].v_g << std::endl;
   std::cout << "lambda=" <<(*list)[i].lambda << std::endl;
   std::cout << "-- Record end" << std::endl;*/

               std::cout << "Object " << j << std::endl;
/*
   std::cout << "a=" <<(*list)[j].a << std::endl;
   std::cout << "e=" <<(*list)[j].e << std::endl;
   std::cout << "i=" <<(*list)[j].i << std::endl;
   std::cout << "omega=" <<(*list)[j].omega << std::endl;
   std::cout << "Omega=" <<(*list)[j].Omega << std::endl;
   std::cout << "-- Observational elements: " << std::endl;
   std::cout << "ra=" <<(*list)[j].ra << std::endl;
   std::cout << "dec=" <<(*list)[j].dec << std::endl;
   std::cout << "v_g=" <<(*list)[j].v_g << std::endl;
   std::cout << "lambda=" <<(*list)[j].lambda << std::endl;
   std::cout << "-- Record end" << std::endl;*/

            }
          
         }
            
         M.push_back(v);
      }
   }
   else if(function.compare("D") == 0) {
      for(i=0; i < ((*list).size()-1); i++) {
         v.clear();
         std::cout << "\r-- Metric matrix " << function << " row " << i << " of " << ((*list).size()-1) << " done." << std::flush;
         for(j=i+1; j < (*list).size(); j++) {
            v.push_back( d_D( &((*list)[i]), &((*list)[j])) );
         }
         M.push_back(v);
      }
   }
   else if(function.compare("V") == 0) {
      for(i=0; i < ((*list).size()-1); i++) {
         v.clear();
         std::cout << "\r-- Metric matrix " << function << " row " << i << " of " << ((*list).size()-1) << " done." << std::flush;
         for(j=i+1; j < (*list).size(); j++) {
            v.push_back( d_V( &((*list)[i]), &((*list)[j])) );
         }
         M.push_back(v);
      }
   }
   else if(function.compare("N") == 0) {
      for(i=0; i < ((*list).size()-1); i++) {
         v.clear();
         std::cout << "\r-- Metric matrix " << function << " row " << i << " of " << ((*list).size()-1) << " done." << std::flush;
         for(j=i+1; j < (*list).size(); j++) {
            v.push_back( d_N( &((*list)[i]), &((*list)[j])) );
         }
         M.push_back(v);
      }
   }
   else if(function.compare("rho2") == 0) {
      for(i=0; i < ((*list).size()-1); i++) {
         v.clear();
         std::cout << "\r-- Metric matrix " << function << " row " << i << " of " << ((*list).size()-1) << " done." << std::flush;
         for(j=i+1; j < (*list).size(); j++) {
            v.push_back( rho2( &((*list)[i]), &((*list)[j])) );
         }
         M.push_back(v);
      }
   }
   else if(function.compare("varrho1") == 0) {
      for(i=0; i < ((*list).size()-1); i++) {
         v.clear();
         std::cout << "\r-- Metric matrix " << function << " row " << i << " of " << ((*list).size()-1) << " done." << std::flush;
         for(j=i+1; j < (*list).size(); j++) {
            v.push_back( varrho1( &((*list)[i]), &((*list)[j])) );
         }
         M.push_back(v);
      }
   }
   
   std::cout << "\r-- Metric matrix " << function << " of size " << ((*list).size()) << " done.                       " << std::endl;

   return M;
}

std::vector<std::vector<double> > augment_metric_matrix(const std::vector<OBSERVATION_record> *list, const std::vector<OBSERVATION_record> *list_addition, std::string function) {
   unsigned int i,j;
   std::vector<double> v;
   std::vector<std::vector<double> > M;
   std::cout << "-- Augment Metric matrix " << function << " row " << 0 << " of " << ((*list).size()-1) << " done, augment list: " << (*list_addition).size() << "." << std::flush;
   if(function.compare("SH") == 0) {
      for(i=0; i < (*list).size(); i++) {
         v.clear();
         std::cout << "\r-- Augment Metric matrix " << function << " row " << i << " of " << ((*list).size()-1) << " done, augment list: " << (*list_addition).size() << "." << std::flush;
         for(j=0; j < (*list_addition).size(); j++) {
            v.push_back( d_SH( &((*list)[i]), &((*list_addition)[j])) );
         }
         M.push_back(v);
      }
   }
   else if(function.compare("D") == 0) {
      for(i=0; i < (*list).size(); i++) {
         v.clear();
         std::cout << "\r-- Augment Metric matrix " << function << " row " << i << " of " << ((*list).size()-1) << " done, augment list: " << (*list_addition).size() << "." << std::flush;
         for(j=0; j < (*list_addition).size(); j++) {
            v.push_back( d_D( &((*list)[i]), &((*list_addition)[j])) );
         }
         M.push_back(v);
      }
   }
   else if(function.compare("V") == 0) {
      for(i=0; i < (*list).size(); i++) {
         v.clear();
         std::cout << "\r-- Augment Metric matrix " << function << " row " << i << " of " << ((*list).size()-1) << " done, augment list: " << (*list_addition).size() << "." << std::flush;
         for(j=0; j < (*list_addition).size(); j++) {
            v.push_back( d_V( &((*list)[i]), &((*list_addition)[j])) );
         }
         M.push_back(v);
      }
   }
   else if(function.compare("N") == 0) {
      for(i=0; i < (*list).size(); i++) {
         v.clear();
         std::cout << "\r-- Augment Metric matrix " << function << " row " << i << " of " << ((*list).size()-1) << " done, augment list: " << (*list_addition).size() << "." << std::flush;
         for(j=0; j < (*list_addition).size(); j++) {
            v.push_back( d_N( &((*list)[i]), &((*list_addition)[j])) );
         }
         M.push_back(v);
      }
   }
   else if(function.compare("rho2") == 0) {
      for(i=0; i < (*list).size(); i++) {
         v.clear();
         std::cout << "\r-- Augment Metric matrix " << function << " row " << i << " of " << ((*list).size()-1) << " done, augment list: " << (*list_addition).size() << "." << std::flush;
         for(j=0; j < (*list_addition).size(); j++) {
            v.push_back( rho2( &((*list)[i]), &((*list_addition)[j])) );
         }
         M.push_back(v);
      }
   }
   else if(function.compare("varrho1") == 0) {
      for(i=0; i < (*list).size(); i++) {
         v.clear();
         std::cout << "\r-- Augment Metric matrix " << function << " row " << i << " of " << ((*list).size()-1) << " done, augment list: " << (*list_addition).size() << "." << std::flush;
         for(j=0; j < (*list_addition).size(); j++) {
            v.push_back( varrho1( &((*list)[i]), &((*list_addition)[j])) );
         }
         M.push_back(v);
      }
   }
   
 std::cout << "\r-- Augment Metric matrix " << function << " " <<  ((*list).size()) << " rows done, augment list size  " << (*list_addition).size() << ".                       " << std::endl;
   return M;
}


/*
std::vector<std::vector<unsigned int> > cluster_analysis_on_file(std::string matrix_file, double d_crit) {
   unsigned int i,j;
   //unsigned int objs = ((*d_matrix).size()) + 1;

   std::ifstream data_matrix;
   data_matrix.open(matrix_file.c_str(), std::ios::in);

   std::vector<std::vector<unsigned int> > clusters;
   std::vector<unsigned int> empty_vec;

   std::vector<node> nodes;
   unsigned int cluster_counter;
   node temp_node;
   for(i=0; i < objs; i++) {
      temp_node.set_node(d_matrix,i);
      nodes.push_back(temp_node);
   }

   cluster_counter = 0;
   for(i=0; i < objs; i++) {
      node_associate(&nodes, &cluster_counter, d_crit, i, -1);
   }
   for(i=0; i < cluster_counter; i++) {
      clusters.push_back(empty_vec);
   }
   for(i=0; i < objs; i++) {
      clusters[nodes[i].cluster_id-1].push_back(i);
   }

   return clusters;
}

int node_associate_on_file(std::vector<node> *nodes, unsigned int *cluster_counter, double d_crit, unsigned int node_id, int mother_cluster, std::ifstream data_matrix) {
   if((*nodes)[node_id].member_of_cluster) {
      return (*nodes)[node_id].cluster_id;
   }
   else if(mother_cluster == -1) {
      (*cluster_counter)++;
      (*nodes)[node_id].member_of_cluster = TRUE;
      (*nodes)[node_id].cluster_id = (*cluster_counter);
   }
   else if(mother_cluster >= 0) {
      (*nodes)[node_id].member_of_cluster = TRUE;
      (*nodes)[node_id].cluster_id = (unsigned int)mother_cluster;
   }

   unsigned int i;
   int ret;

   for(i=0; i < (*nodes)[node_id].d_vec.size(); i++) {
      if((*nodes)[node_id].d_vec[i] < d_crit && (*nodes)[node_id].d_vec[i] < d_crit != 0) {
         ret = node_associate(nodes, cluster_counter, d_crit, i, (*nodes)[node_id].cluster_id);
         if(ret > (*cluster_counter)) {
            std::cout << "Something went wrong in cluster association" << std::endl;
            return -1;
         }
      }
   }

   return (*nodes)[node_id].cluster_id;
}
*/

double d_SH(const OBSERVATION_record *orb1, const OBSERVATION_record *orb2) {
   double e_a = (*orb1).e;
   double e_b = (*orb2).e;
   double q_a = (*orb1).a*(1 - e_a);
   double q_b = (*orb2).a*(1 - e_b);
   double i_a = (*orb1).i*PI/180;
   double i_b = (*orb2).i*PI/180;
   double omega_a = (*orb1).omega*PI/180;
   double omega_b = (*orb2).omega*PI/180;
   double Omega_a = (*orb1).Omega*PI/180;
   double Omega_b = (*orb2).Omega*PI/180;

   //i_a = i_a - fmin(i_a,i_b);
   //i_b = i_b - fmin(i_a,i_b);

   double I_21 = 2*asin(0.5*sqrt(pow(2*sin((i_b-i_a)*0.5),2)+sin(i_a)*sin(i_b)*pow(2*sin((Omega_b-Omega_a)*0.5),2)));
   double arcsin_arg = cos((i_b+i_a)*0.5)*sin((Omega_b-Omega_a)*0.5)/cos(I_21*0.5);
   if(arcsin_arg > 1) {
      arcsin_arg = 1;
   }
   else if(arcsin_arg < -1) {
      arcsin_arg = -1;
   }
   if((Omega_b-Omega_a) > PI || (Omega_b-Omega_a) < -PI) { arcsin_arg = (-1)*arcsin_arg; }
   double pi_21 = omega_b - omega_a + 2*asin(arcsin_arg);
   //double pi_21 = Omega_b + omega_b - (Omega_a + omega_a);

   double d = sqrt(pow(e_b-e_a,2) + pow(q_b-q_a,2) + pow(2*sin(I_21*0.5),2) + pow(((e_b+e_a)/2)*(2*sin(pi_21*0.5)),2));
   return d;
}

double d_D(const OBSERVATION_record *orb1, const OBSERVATION_record *orb2) {
   double e_a = (*orb1).e;
   double e_b = (*orb2).e;
   double q_a = (*orb1).a*(1 - e_a);
   double q_b = (*orb2).a*(1 - e_b);
   double i_a = (*orb1).i*PI/180;
   double i_b = (*orb2).i*PI/180;
   double omega_a = (*orb1).omega*PI/180;
   double omega_b = (*orb2).omega*PI/180;
   double Omega_a = (*orb1).Omega*PI/180;
   double Omega_b = (*orb2).Omega*PI/180;

   double lambda_1 = Omega_a + atan(cos(i_a)*tan(omega_a));
   if(cos(omega_a) < 0) { lambda_1 += PI; }
   double lambda_2 = Omega_b + atan(cos(i_b)*tan(omega_b));
   if(cos(omega_b) < 0) { lambda_2 += PI; }
   double beta_1 = asin(sin(i_a)*sin(omega_a));
   double beta_2 = asin(sin(i_b)*sin(omega_b));

   double I_21 = 2*asin(0.5*sqrt(pow(2*sin((i_b-i_a)/2),2)+sin(i_a)*sin(i_b)*pow(2*sin((Omega_b-Omega_a)/2),2)));
   double Theta_21 = acos(sin(beta_1)*sin(beta_2)+cos(beta_1)*cos(beta_2)*cos(lambda_2-lambda_1));

   double d = sqrt(pow((e_b-e_a)/(e_b+e_a),2) + pow((q_b-q_a)/(q_b+q_a),2) + pow(I_21/PI,2) + pow((e_b+e_a)/2*(Theta_21/PI),2));
   return d;
}

double rho2(const OBSERVATION_record *orb1, const OBSERVATION_record *orb2) {
   double e_a = (*orb1).e;
   double e_b = (*orb2).e;
   double q_a = (*orb1).a*(1 - e_a);
   double q_b = (*orb2).a*(1 - e_b);
   double i_a = (*orb1).i*PI/180;
   double i_b = (*orb2).i*PI/180;
   double omega_a = (*orb1).omega*PI/180;
   double omega_b = (*orb2).omega*PI/180;
   double Omega_a = (*orb1).Omega*PI/180;
   double Omega_b = (*orb2).Omega*PI/180;

   double a_a = q_a/(1-e_a);
   double a_b = q_b/(1-e_b);

   double eta_a = sqrt(1-e_a*e_a);
   double eta_b = sqrt(1-e_b*e_b);

   std::vector<double> P1;
   P1.push_back(cos(omega_a)*cos(Omega_a) - cos(i_a)*sin(omega_a)*sin(Omega_a));
   P1.push_back(cos(omega_a)*sin(Omega_a) + cos(i_a)*sin(omega_a)*cos(Omega_a));
   P1.push_back(sin(i_a)*sin(omega_a));

   std::vector<double> P2;
   P2.push_back(cos(omega_b)*cos(Omega_b) - cos(i_b)*sin(omega_b)*sin(Omega_b));
   P2.push_back(cos(omega_b)*sin(Omega_b) + cos(i_b)*sin(omega_b)*cos(Omega_b));
   P2.push_back(sin(i_b)*sin(omega_b));

   std::vector<double> S1;
   S1.push_back(eta_a*(-sin(omega_a)*cos(Omega_a) - cos(i_a)*cos(omega_a)*sin(Omega_a)));
   S1.push_back(eta_a*(-sin(omega_a)*sin(Omega_a) + cos(i_a)*cos(omega_a)*cos(Omega_a)));
   S1.push_back(eta_a*(sin(i_a)*cos(omega_a)));

   std::vector<double> S2;
   S2.push_back(eta_b*(-sin(omega_b)*cos(Omega_b) - cos(i_b)*cos(omega_b)*sin(Omega_b)));
   S2.push_back(eta_b*(-sin(omega_b)*sin(Omega_b) + cos(i_b)*cos(omega_b)*cos(Omega_b)));
   S2.push_back(eta_b*(sin(i_b)*cos(omega_b)));

   double alpha1 = a_a/a_b;
   double alpha2 = a_b/a_a;

   double W_0 = 0.25*(2*(alpha1 + alpha2) + alpha1*e_a*e_a + alpha2*e_b*e_b - 4*dot_product(P1,P2)*e_a*e_b);

   double W_5 = 0.5*(-dot_product(P1,P2));
   double W_6 = 0.5*(-dot_product(P1,S2));
   double W_7 = 0.5*(-dot_product(P2,S1));
   double W_8 = 0.5*(-dot_product(S1,S2));

   double d = sqrt(2*a_a*a_b*(W_0 - sqrt((W_5 + W_8)*(W_5 + W_8) + (W_6 - W_7)*(W_6 - W_7))));

   if(!is_finite(d)) {
      d = 0;
   }

   return d;
}

double d_V(const OBSERVATION_record *orb1, const OBSERVATION_record *orb2) {
   double d = 0;
   return d;
}

double d_N(const OBSERVATION_record *orb1, const OBSERVATION_record *orb2) {
   double d = 0;
   return d;
}

double varrho1(const OBSERVATION_record *orb1, const OBSERVATION_record *orb2) {
   double e_a = (*orb1).e;
   double e_b = (*orb2).e;
   double q_a = (*orb1).a*(1 - e_a);
   double q_b = (*orb2).a*(1 - e_b);
   double i_a = (*orb1).i*PI/180;
   double i_b = (*orb2).i*PI/180;
   double omega_a = (*orb1).omega*PI/180;
   double omega_b = (*orb2).omega*PI/180;
   double Omega_a = (*orb1).Omega*PI/180;
   double Omega_b = (*orb2).Omega*PI/180;

   double L1 = 1;
   double L2 = 1;

   double a_a = q_a/(1-e_a);
   double a_b = q_b/(1-e_b);
   double b_a = a_a*(1-e_a*e_a);
   double b_b = a_b*(1-e_b*e_b);
   double p_a = 2*b_a*b_a/a_a;
   double p_b = 2*b_b*b_b/a_b;

   double cos_xi = cos(i_a)*cos(i_b) + sin(i_a)*sin(i_b)*cos(Omega_a - Omega_b);
   double cos_eta = (cos(omega_a)*cos(omega_b) + cos(i_a)*cos(i_b)*sin(omega_a)*sin(omega_b))*cos(Omega_a-Omega_b) + (cos(i_b)*cos(omega_a)*sin(omega_b) - cos(i_a)*sin(omega_a)*cos(omega_b))*sin(Omega_a-Omega_b) + sin(i_a)*sin(i_b)*sin(omega_a)*sin(omega_b);

   double d = sqrt(1/L1*(p_a + p_b - 2*sqrt(p_a*p_b)*cos_xi) + (e_a*e_a + e_b*e_b - 2*e_a*e_b*cos_eta) + L2*L2/4*(1/a_a - 1/a_b)*(1/a_a - 1/a_b));
   if(!is_finite(d)) {
      d = 0;
   }
   return d;
}

