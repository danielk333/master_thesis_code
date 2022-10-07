#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "define.hh"
#include "functions.hh"


int save_mat(std::string out_file_name, std::vector<std::vector<double> > *M) {

        unsigned int i;
        unsigned int j;

        std::ofstream out_pos(out_file_name.c_str(), std::ios::app);
        out_pos << std::setprecision(20) << std::scientific;
        if(out_pos.fail()) {
            return DATA_SAVE_FAILED;
        }

        for(j = 0; j < (*M).size(); j++) {
            for(i = 0; i < (*M)[j].size(); i++) {
                out_pos << (((*M)[j])[i]);
                if(i < ((*M)[j].size() - 1)) {
                    out_pos << " ";
                }
                else {
                    out_pos << "\r\n";
                }
            }
        }

        out_pos.close();

        return SIM_OK;
}


int save(const char *out, std::vector<std::vector<std::vector<double> > > *q, std::vector<std::vector<std::vector<double> > > *p, std::vector<std::vector<double> > *m_vector, std::vector<std::vector<int> > *type_vec, std::vector<double>* time_vec, OPTIONS *opt, SIMULATION *sim) {
        unsigned int i;
        unsigned int j;
        std::string out_file_name = out;
        std::ofstream out_pos(out_file_name.c_str(), std::ios::app);
        out_pos << std::setprecision((*opt).precision) << std::scientific;
        if(out_pos.fail()) {
            return DATA_SAVE_FAILED;
        }

        for(j = 0; j < (*q).size(); j++) {
            for(i = 0; i < (*q)[j].size(); i++) {
                out_pos << (*time_vec)[j] << " ";
                out_pos << (*q)[j].size() << "\r\n";

                out_pos << i << " ";
                out_pos << (*type_vec)[j][i] << " ";
                out_pos << (*m_vector)[j][i] << "\r\n";
                
                out_pos << ((*q)[j])[i][0] << " ";
                out_pos << ((*q)[j])[i][1] << " ";
                out_pos << ((*q)[j])[i][2] << "\r\n";

                if((*type_vec)[j][i] == 1) {
                    out_pos << ((*p)[j])[i][0]/(*m_vector)[j][i] << " ";
                    out_pos << ((*p)[j])[i][1]/(*m_vector)[j][i] << " ";
                    out_pos << ((*p)[j])[i][2]/(*m_vector)[j][i];
                }
                else {
                    out_pos << ((*p)[j])[i][0]/(*m_vector)[j][i] << " ";
                    out_pos << ((*p)[j])[i][1]/(*m_vector)[j][i] << " ";
                    out_pos << ((*p)[j])[i][2]/(*m_vector)[j][i];
                }
                out_pos << "\r\n";
            }
        }

        out_pos.close();

        return SIM_OK;

/*
        unsigned int i;
        unsigned int j;
        std::vector<double> kepler_elements;
        std::string out_file_name = out;

        std::ofstream out_pos(out_file_name.c_str(), std::ios::app);
        out_pos << std::setprecision((*opt).precision);
        if(out_pos.fail()) {
            return DATA_SAVE_FAILED;
        }

        if((*opt).output_type == 0) {
        	for(j = 0; j < (*q).size(); j++) {
        			for(i = 0; i < (*q)[j].size(); i++) {
        				if((*type_vec)[i] != -1) {
        					out_pos << i << " " << ((double)((*sim).loops_performed-((*q).size()-j-1)*(*opt).save_interval))*(*opt).dt/((double)(3600*24*365.25)) << " " << (((*q)[j])[i])[0] << " " << (((*q)[j])[i])[1] << " " << (((*q)[j])[i])[2] << " " << (((*p)[j])[i])[0] << " " << (((*p)[j])[i])[1] << " " << (((*p)[j])[i])[2];
        					out_pos << "\r\n";
        				}
        			}
        	}
        }
        else if((*opt).output_type == 1) {
        	for(j = 0; j < (*q).size(); j++) {
        			for(i = 0; i < (*q)[j].size(); i++) {
        				if((*type_vec)[i] != -1) {
        					kepler_elements.clear();
        					if(i==0) {
        						kepler_elements.push_back(0);
        						kepler_elements.push_back(0);
        						kepler_elements.push_back(0);

        						kepler_elements.push_back(0);
        						kepler_elements.push_back(0);
        						kepler_elements.push_back(0);
        					}
        					else {
        						//kepler_elements = qp_to_kepler_sun(((*q)[j])[0], ((*p)[j])[0], ((*q)[j])[i], ((*p)[j])[i], (*m_vector)[i], (*m_vector)[0]);
        					}
        					out_pos << i << " " << ((double)((*sim).loops_performed-((*q).size()-j-1)*(*opt).save_interval))*(*opt).dt/((double)(3600*24*365.25)) << " " << kepler_elements[0]/AU << " " << kepler_elements[1] << " " << kepler_elements[2] << " " << kepler_elements[3] << " " << kepler_elements[4] << " " << kepler_elements[5];
        					out_pos << "\r\n";
        				}
        			}
        	}
        }
        else if((*opt).output_type == 2) {
        	for(j = 0; j < (*q).size(); j++) {
        			for(i = 0; i < (*q)[j].size(); i++) {
        				if((*type_vec)[i] != -1) {
        					kepler_elements.clear();
        					if(i==0) {
        					    kepler_elements.push_back(0);
        					    kepler_elements.push_back(0);
        					    kepler_elements.push_back(0);

        					    kepler_elements.push_back(0);
        					    kepler_elements.push_back(0);
       					       	kepler_elements.push_back(0);
           					}
         					else {
         						//kepler_elements = qp_to_kepler_sun(((*q)[j])[0], ((*p)[j])[0], ((*q)[j])[i], ((*p)[j])[i], (*m_vector)[i], (*m_vector)[0]);
        					}

        					out_pos << i << " " << ((double)((*sim).loops_performed-((*q).size()-j-1)*(*opt).save_interval))*(*opt).dt/((double)(3600*24*365.25)) << " " << kepler_elements[0]/AU << " " << kepler_elements[1] << " " << kepler_elements[2] << " " << kepler_elements[3] << " " << kepler_elements[4] << " " << kepler_elements[5];
        					out_pos << " " << (((*q)[j])[i])[0] << " " << (((*q)[j])[i])[1] << " " << (((*q)[j])[i])[2] << " " << (((*p)[j])[i])[0] << " " << (((*p)[j])[i])[1] << " " << (((*p)[j])[i])[2];
        					out_pos << "\r\n";
        				}
        			}
        	}
        }
        out_pos.close();
*/
        return SIM_OK;

}

 int save_diagnostics(const char *out, std::vector<double> E_vec, std::vector<double> T_vec) {

        unsigned int i;

        std::string out_file_name = out;

        std::ofstream out_pos(out_file_name.c_str(), std::ios::app);
        if(out_pos.fail()) {
            return DIAG_SAVE_FAILED;
        }

        for(i = 0; i < E_vec.size(); i++) {
        	out_pos << T_vec[i] << " " << E_vec[i] << "\r\n";
        }

        out_pos.close();

        return SIM_OK;

}

 int save_settings(const char *out, OPTIONS *opt, SIMULATION *sim) {

        unsigned int i;

        std::string out_file_name = out;

        std::ofstream out_pos(out_file_name.c_str(), std::ios::app);
        if(out_pos.fail()) {
            return SETTINGS_SAVE_FAILED;
        }

    out_pos << (*opt).loops << "\r\n";
    out_pos << (*opt).precision << "\r\n";
    out_pos << (*opt).memory_alloc << "\r\n";
    out_pos << (*opt).output_type << "\r\n";
    out_pos << (*opt).save_interval << "\r\n";
    out_pos << (*opt).coord << "\r\n";
    out_pos << (*opt).energy_diag << "\r\n";
    out_pos << (*opt).sun_collision << "\r\n";
    out_pos << (*opt).hill_sphere << "\r\n";
    out_pos << (*opt).ejection_criteria << "\r\n";
    out_pos << (*opt).integrator << "\r\n";
    out_pos << (*opt).dt << "\r\n";

    out_pos << (*opt).body_n << "\r\n";
    out_pos << (*opt).massive_body_n << "\r\n";
    out_pos << (*opt).memory_clear_interval << "\r\n";

    out_pos << (*sim).loops_performed << "\r\n";

    out_pos << (*sim).save_number << "\r\n";
    out_pos << (*sim).next_clear << "\r\n";
    out_pos << (*sim).diag_counter << "\r\n";
    out_pos << (*sim).save_counter << "\r\n";
    
    out_pos << (*sim).t_elapsed << "\r\n";

    out_pos << (*sim).removed_bodies;

        out_pos.close();

        return SIM_OK;

}

int add_save_data_to_memory(OPTIONS *opt, std::vector<std::vector<std::vector<double> > > *q_vec_save, std::vector<std::vector<std::vector<double> > > *p_vec_save, std::vector<std::vector<double> > *q_vec, std::vector<std::vector<double> > *p_vec, std::vector<double> *m_vector) {
   /* std::vector<std::vector<double> > q_vec_temp;
    std::vector<std::vector<double> > p_vec_temp;
    std::vector<double> C_v;
    std::vector<double> C;

        unsigned int j;

            q_vec_temp = (*q_vec);
            p_vec_temp = (*p_vec);
            if((*opt).coord == 0) {
                for(j=0; j < (*q_vec).size(); j++) {
                    (q_vec_temp[j])[0] = ((*q_vec)[j])[0] - ((*q_vec)[0])[0];
                    (q_vec_temp[j])[1] = ((*q_vec)[j])[1] - ((*q_vec)[0])[1];
                    (q_vec_temp[j])[2] = ((*q_vec)[j])[2] - ((*q_vec)[0])[2];
                    (p_vec_temp[j])[0] = ((*p_vec)[j])[0] - (*m_vector)[j]/(*m_vector)[0]*((*p_vec)[0])[0];
                    (p_vec_temp[j])[1] = ((*p_vec)[j])[1] - (*m_vector)[j]/(*m_vector)[0]*((*p_vec)[0])[1];
                    (p_vec_temp[j])[2] = ((*p_vec)[j])[2] - (*m_vector)[j]/(*m_vector)[0]*((*p_vec)[0])[2];
                }
            }
            else { 
                C = barycentre((*q_vec),*m_vector);
                C_v = barycentre_v((*p_vec),*m_vector);
                for(j=0; j < (*q_vec).size(); j++) {
                    (q_vec_temp[j])[0] = ((*q_vec)[j])[0] - C[0];
                    (q_vec_temp[j])[1] = ((*q_vec)[j])[1] - C[1];
                    (q_vec_temp[j])[2] = ((*q_vec)[j])[2] - C[2];
                    (p_vec_temp[j])[0] = ((*p_vec)[j])[0] - C_v[0]*(*m_vector)[j];
                    (p_vec_temp[j])[1] = ((*p_vec)[j])[1] - C_v[1]*(*m_vector)[j];
                    (p_vec_temp[j])[2] = ((*p_vec)[j])[2] - C_v[2]*(*m_vector)[j];
                }
            }
            (*q_vec_save).push_back(q_vec_temp);
            (*p_vec_save).push_back(p_vec_temp);
*/
            return SIM_OK;
}
