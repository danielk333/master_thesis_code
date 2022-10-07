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

        return CALC_OK;
}

int save_batch(std::string file_name, std::string temp_file_name, std::vector<std::vector<double> > *M, unsigned int row_skip) {
        unsigned int i;
        unsigned int j;
    std::string buffer;
    std::string buffer_temp;

    std::ifstream file;
    std::ofstream temp_file;

    std::ofstream file_write;
    std::ifstream temp_file_read;

    file.open(file_name.c_str(), std::ios::in);
    if(file.fail()) {
        return MATRIX_LOAD_FAILED;
    }

    temp_file.open(temp_file_name.c_str(), std::ofstream::out | std::ofstream::trunc);
    temp_file << std::setprecision(20) << std::scientific;
    if(temp_file.fail()) {
        return DATA_SAVE_FAILED;
    }

    for(j=0; j < row_skip; j++) {
        buffer.clear();
        getline(file,buffer);
        if(buffer.size() > 0) {
            temp_file << buffer;
        }
    }

    for(j = 0; j < (*M).size(); j++) {
        buffer.clear();
        getline(file,buffer);

        if(buffer.size() > 0) {
            buffer.erase(buffer.end()-1,buffer.end());

            temp_file << buffer << " ";
            for(i = 0; i < (*M)[j].size(); i++) {
                temp_file << (((*M)[j])[i]);
                if(i < (*M)[j].size()-1) {
                    temp_file << " ";
                }
                else {
                    temp_file << "\r\n";
                }
            }
        }
    }
    
    while(!file.eof()) {
        buffer.clear();
        getline(file,buffer);
        if(buffer.size() > 0) {
            temp_file << buffer;
        }
    }

    temp_file.close();
    file.close();

    file_write.open(file_name.c_str(), std::ofstream::out | std::ofstream::trunc);
    if(file_write.fail()) {
        return DATA_SAVE_FAILED;
    }

    temp_file_read.open(temp_file_name.c_str(),  std::ios::in);
    if(temp_file_read.fail()) {
        return MATRIX_LOAD_FAILED;
    }

    while(!temp_file_read.eof()) {
        buffer.clear();
        getline(temp_file_read,buffer);
        if(buffer.size() > 0) {
            file_write << buffer;
        }
    }
    temp_file_read.close();
    file_write.close();

    return CALC_OK;
}

int print_mercury6_big(std::string mercury6_big, std::vector<std::vector<double> > *temp_q, std::vector<std::vector<double> > *temp_v, std::vector<double> *temp_m, double DATE) {
	        std::ofstream out_stream;
	        out_stream.open(mercury6_big.c_str(), std::ios::app);
	        std::string line_end = "\n";

	    	if(out_stream.fail()) {
            	return mercury6_big_DATA_SAVE_FAILED;
        	}

            out_stream << std::uppercase << std::setprecision(17) << std::scientific;
            out_stream << ")O+_06 Big-body initial data  (WARNING: Do not delete this line!!)" << line_end;
            out_stream << ") Lines beginning with `)' are ignored." << line_end;
            out_stream << ")---------------------------------------------------------------------" << line_end;
            out_stream << " style (Cartesian, Asteroidal, Cometary) = Cartesian" << line_end;
            out_stream << " epoch (in days) = " << DATE << line_end;
            out_stream << ")---------------------------------------------------------------------" << line_end;
            unsigned int i = 1;
            out_stream << " " << "MERCURY m=" << (*temp_m)[i]/Msol << " r=3.D0 d=5.43" << line_end;
            out_stream << " " << (*temp_q)[i][0] << " " << (*temp_q)[i][1] << " " << (*temp_q)[i][2] << line_end;
            out_stream << " " << (*temp_v)[i][0]/AU*86400 << " " << (*temp_v)[i][1]/AU*86400 << " " << (*temp_v)[i][2]/AU*86400 <<  line_end;          
            out_stream << " " << "0. 0. 0." << line_end; i++;
            out_stream << " " << "VENUS   m=" << (*temp_m)[i]/Msol << " r=3.d0 d=5.24" << line_end;
            out_stream << " " << (*temp_q)[i][0] << " " << (*temp_q)[i][1] << " " << (*temp_q)[i][2] <<  line_end;
            out_stream << " " << (*temp_v)[i][0]/AU*86400 << " " << (*temp_v)[i][1]/AU*86400 << " " << (*temp_v)[i][2]/AU*86400 <<  line_end;          
            out_stream << " " << "0. 0. 0." << line_end; i++;
            out_stream << " " << "EARTH   m=" << (*temp_m)[i]/Msol << " r=3.d0 d=5.52" << line_end;
            out_stream << " " << (*temp_q)[i][0] << " " << (*temp_q)[i][1] << " " << (*temp_q)[i][2] <<  line_end;
            out_stream << " " << (*temp_v)[i][0]/AU*86400 << " " << (*temp_v)[i][1]/AU*86400 << " " << (*temp_v)[i][2]/AU*86400 <<  line_end;          
            out_stream << " " << "0. 0. 0." << line_end; i++;
            out_stream << " " << "MARS    m=" << (*temp_m)[i]/Msol << " r=3.d0 d=3.94" << line_end;
            out_stream << " " << (*temp_q)[i][0] << " " << (*temp_q)[i][1] << " " << (*temp_q)[i][2] <<  line_end;
            out_stream << " " << (*temp_v)[i][0]/AU*86400 << " " << (*temp_v)[i][1]/AU*86400 << " " << (*temp_v)[i][2]/AU*86400 <<  line_end;          
            out_stream << " " << "0. 0. 0." << line_end; i++;
            out_stream << " " << "JUPITER m=" << (*temp_m)[i]/Msol << " r=3.d0 d=1.33" << line_end;
            out_stream << " " << (*temp_q)[i][0] << " " << (*temp_q)[i][1] << " " << (*temp_q)[i][2] <<  line_end;
            out_stream << " " << (*temp_v)[i][0]/AU*86400 << " " << (*temp_v)[i][1]/AU*86400 << " " << (*temp_v)[i][2]/AU*86400 <<  line_end;          
            out_stream << " " << "0. 0. 0." << line_end; i++;
            out_stream << " " << "SATURN  m=" << (*temp_m)[i]/Msol << " r=3.d0 d=0.70" << line_end;
            out_stream << " " << (*temp_q)[i][0] << " " << (*temp_q)[i][1] << " " << (*temp_q)[i][2] <<  line_end;
            out_stream << " " << (*temp_v)[i][0]/AU*86400 << " " << (*temp_v)[i][1]/AU*86400 << " " << (*temp_v)[i][2]/AU*86400 <<  line_end;          
            out_stream << " " << "0. 0. 0." << line_end; i++;
            out_stream << " " << "URANUS  m=" << (*temp_m)[i]/Msol << " r=3.d0 d=1.30" << line_end;
            out_stream << " " << (*temp_q)[i][0] << " " << (*temp_q)[i][1] << " " << (*temp_q)[i][2] <<  line_end;
            out_stream << " " << (*temp_v)[i][0]/AU*86400 << " " << (*temp_v)[i][1]/AU*86400 << " " << (*temp_v)[i][2]/AU*86400 <<  line_end;          
            out_stream << " " << "0. 0. 0." << line_end; i++;
            out_stream << " " << "NEPTUNE m=" << (*temp_m)[i]/Msol << " r=3.d0 d=1.76" << line_end;
            out_stream << " " << (*temp_q)[i][0] << " " << (*temp_q)[i][1] << " " << (*temp_q)[i][2] <<  line_end;
            out_stream << " " << (*temp_v)[i][0]/AU*86400 << " " << (*temp_v)[i][1]/AU*86400 << " " << (*temp_v)[i][2]/AU*86400 <<  line_end;          
            out_stream << " " << "0. 0. 0."; 

            out_stream.close();
            out_stream.clear();

            return CALC_OK;
}

int print_mercury6_small_pb(std::string mercury6_small, std::vector<std::vector<double> > *temp_q, std::vector<std::vector<double> > *temp_v, std::vector<double> *temp_m, std::vector<double> *temp_t, double rho) {
            std::ofstream out_stream;
            out_stream.open(mercury6_small.c_str(), std::ios::app);
            std::string line_end = "\n";
            unsigned int i,used_particle_number_mercury;
            bool TRUTH_CHECKER;

            if(out_stream.fail()) {
                return mercury6_small_DATA_SAVE_FAILED;
            }
            out_stream << ")O+_06 Small-body initial data  (WARNING: Do not delete this line!!)" << line_end;
            out_stream << ") Lines beginning with `)' are ignored." << line_end;
            out_stream << ")---------------------------------------------------------------------" << line_end;
            out_stream << " style (Cartesian, Asteroidal, Cometary) = Car" << line_end;
            out_stream << ")---------------------------------------------------------------------" << line_end;
            
            out_stream << " " << "PB" << "      ep=" << (*temp_t)[0] << " p1=" << (*temp_m)[i] << " p2=" << rho << line_end;
            out_stream << " " << (*temp_q)[0][0]/AU << " " << (*temp_q)[0][1]/AU << " " << (*temp_q)[0][2]/AU;
            out_stream << " " << (*temp_v)[0][0]/AU*86400 << " " << (*temp_v)[0][1]/AU*86400 << " " << (*temp_v)[0][2]/AU*86400;
            out_stream << " " << "0. 0. 0." << line_end;

            out_stream.close();
            out_stream.clear();

            return CALC_OK;
}

int print_mercury6_small(std::string mercury6_small, std::vector<std::vector<double> > *temp_q, std::vector<std::vector<double> > *temp_v, std::vector<double> *temp_m, std::vector<double> *temp_t, double rho, unsigned long long int MAX, unsigned long long int *particle_number_mercury) {
            std::ofstream out_stream;
            out_stream.open(mercury6_small.c_str(), std::ios::app);
            std::string line_end = "\n";
            unsigned int i,used_particle_number_mercury;
            bool TRUTH_CHECKER;

            if(out_stream.fail()) {
                return mercury6_small_DATA_SAVE_FAILED;
            }
            out_stream << ")O+_06 Small-body initial data  (WARNING: Do not delete this line!!)" << line_end;
            out_stream << ") Lines beginning with `)' are ignored." << line_end;
            out_stream << ")---------------------------------------------------------------------" << line_end;
            out_stream << " style (Cartesian, Asteroidal, Cometary) = Car" << line_end;
            out_stream << ")---------------------------------------------------------------------" << line_end;
            
            out_stream << " " << "PB" << "      ep=" << (*temp_t)[0] << " p1=" << (*temp_m)[0] << " p2=" << rho << line_end;
            out_stream << " " << (*temp_q)[0][0]/AU << " " << (*temp_q)[0][1]/AU << " " << (*temp_q)[0][2]/AU;
            out_stream << " " << (*temp_v)[0][0]/AU*86400 << " " << (*temp_v)[0][1]/AU*86400 << " " << (*temp_v)[0][2]/AU*86400;
            out_stream << " " << "0. 0. 0." << line_end;

            (*particle_number_mercury) = 0;
            if((*temp_q).size()-1 > MAX) {
                used_particle_number_mercury = MAX;
            }
            else {
                used_particle_number_mercury = (*temp_q).size()-1;
            }
            for(i=1; i <= used_particle_number_mercury; i++) { //(*temp_q).size()
                TRUTH_CHECKER = TRUE;
                TRUTH_CHECKER = TRUTH_CHECKER && is_finite((*temp_q)[i][0]);
                TRUTH_CHECKER = TRUTH_CHECKER && is_finite((*temp_q)[i][1]);
                TRUTH_CHECKER = TRUTH_CHECKER && is_finite((*temp_q)[i][2]);
                TRUTH_CHECKER = TRUTH_CHECKER && is_finite((*temp_v)[i][0]);
                TRUTH_CHECKER = TRUTH_CHECKER && is_finite((*temp_v)[i][1]);
                TRUTH_CHECKER = TRUTH_CHECKER && is_finite((*temp_v)[i][2]);
                if(TRUTH_CHECKER) {
                    (*particle_number_mercury)++;
                    out_stream << " " << "OBJ" << (i-1) << "     ep=" << (*temp_t)[i] << " p1=" << (*temp_m)[i] << " p2=" << rho << line_end;
                    out_stream << " " << (*temp_q)[i][0]/AU << " " << (*temp_q)[i][1]/AU << " " << (*temp_q)[i][2]/AU;
                    out_stream << " " << (*temp_v)[i][0]/AU*86400 << " " << (*temp_v)[i][1]/AU*86400 << " " << (*temp_v)[i][2]/AU*86400;
                    out_stream << " " << "0. 0. 0." << line_end;
                }
                else {
                    std::cout << "SOMETHING HORRIBLE HAPPEND TO THE TEST PARTICLE DATA! Skipped object in mercury small.in file." << std::endl;
                }
            }
            // p2 = bulk density
            // p1 = mass

            out_stream.close();
            out_stream.clear();

            return CALC_OK;
}

int print_mercury6_param(std::string mercury6_param, OPTIONS *opt, double start_time) {
            std::ofstream out_stream;
            out_stream.open(mercury6_param.c_str(), std::ios::app);
            std::string line_end = "\n";
            double end_time,out_step;
            if(out_stream.fail()) {
                return mercury6_big_DATA_SAVE_FAILED;
            }
            std::string char_bool;
            if((*opt).merc_pr == 1) {
                char_bool = "yes";
            }
            else {
                char_bool = "no";
            }
            std::string integ;
            if((*opt).merc_int == 0) {
                integ = "hyb";
            }
            else if((*opt).merc_int == 1) {
                integ = "BS";
            }
            else {
                integ = "hyb";
            }

            if((*opt).end_date == 0) {
                end_time = start_time + (*opt).integration_time*365.25;
            }
            else {
                end_time = (*opt).end_date;
            }

            if((*opt).snapshot_date != 0) {
                out_step = 4.0;
            }
            else {
                out_step = 365.25;
                //out_step = 1.0;
            }

            out_stream << std::uppercase << std::setprecision(20) << std::showpoint;
            out_stream << ")O+_06 Integration parameters  (WARNING: Do not delete this line!!)" << line_end;
            out_stream << ") Lines beginning with `)' are ignored." << line_end;
            out_stream << ")---------------------------------------------------------------------" << line_end;
            out_stream << ") Important integration parameters:" << line_end;
            out_stream << ")---------------------------------------------------------------------" << line_end;
            out_stream << " algorithm (MVS, BS, BS2, RADAU, HYBRID etc) = " << integ << line_end;
            out_stream << " start time (days)= " << start_time << line_end;
            out_stream << " stop time (days) = " << end_time << line_end;
            out_stream << " output interval (days) = " << out_step << line_end;
            out_stream << " timestep (days) = " << (*opt).merc_step << line_end;
            out_stream << " accuracy parameter=1.d-12" << line_end;
            out_stream << ")---------------------------------------------------------------------" << line_end;
            out_stream << ") Integration options:" << line_end;
            out_stream << ")---------------------------------------------------------------------" << line_end;
            out_stream << " stop integration after a close encounter = no" << line_end;
            out_stream << " allow collisions to occur = no" << line_end;
            out_stream << " include collisional fragmentation = no" << line_end;
            out_stream << " express time in days or years = years" << line_end;
            out_stream << " express time relative to integration start time = no" << line_end;
            out_stream << " output precision = medium" << line_end;
            out_stream << " global parameter1 = 0.d0" << line_end;
            out_stream << " global parameter2 = 0.d0" << line_end;
            out_stream << " include user-defined force = " << char_bool << line_end;
            out_stream << ")---------------------------------------------------------------------" << line_end;
            out_stream << ") These parameters do not need to be adjusted often:" << line_end;
            out_stream << ")---------------------------------------------------------------------" << line_end;
            out_stream << " ejection distance (AU)= " << (*opt).merc_eject << line_end;
            out_stream << " radius of central body (AU) = 0.005" << line_end;
            out_stream << " central mass (solar) = 1.0" << line_end;
            out_stream << " central J2 = 0" << line_end;
            out_stream << " central J4 = 0" << line_end;
            out_stream << " central J6 = 0" << line_end;
            out_stream << " global parameter3 = 0.d0" << line_end;
            out_stream << " global parameter4 = 0.d0" << line_end;
            out_stream << " Hybrid integrator changeover (Hill radii) = " << (*opt).merc_close_int << line_end;
            out_stream << " number of timesteps between data dumps = 500" << line_end;
            out_stream << " number of timesteps between periodic effects = 100" << line_end;

            out_stream.close();
            out_stream.clear();

            return CALC_OK;
}



int print_mercury6_element(std::string mercury6_element, std::vector<std::string> bodies,unsigned int format) {
            std::ofstream out_stream;
            out_stream.open(mercury6_element.c_str(), std::ios::app);
            std::string line_end = "\n";
            unsigned int i;
            
            if(out_stream.fail()) {
                return mercury6_big_DATA_SAVE_FAILED;
            }

            std::string form;
            switch(format) {
                case 0:
                    form = " a8.5 e8.6 i8.4 g8.4 n8.4 l8.4 f8.4 m13e ";
                    break;
                case 1:
                    form = " x13e y13e z13e u13e v13e w13e ";
                    break;
            }

out_stream << ")O+_06 element  (WARNING: Do not delete this line!!)" << line_end;
out_stream << ") Lines beginning with `)' are ignored." << line_end;
out_stream << ")---------------------------------------------------------------------" << line_end;
out_stream << " number of input files = 1" << line_end;
out_stream << ")---------------------------------------------------------------------" << line_end;
out_stream << ") List the input files, one per line" << line_end;
out_stream << " xv.out" << line_end;
out_stream << ")---------------------------------------------------------------------" << line_end;
out_stream << " type of elements (central body, barycentric, Jacobi) = Cen" << line_end;
out_stream << " minimum interval between outputs (days) = 1.0d0" << line_end;
out_stream << " express time in days or years = days" << line_end;
out_stream << " express time relative to integration start time = yes" << line_end;
out_stream << ")---------------------------------------------------------------------" << line_end;
out_stream << ") Output format? (e.g. a8.4 => semi-major axis with 8 digits & 4 dec. places)" << line_end;
out_stream << form << line_end;
out_stream << ")---------------------------------------------------------------------" << line_end;
out_stream << ") Which bodies do you want? (List one per line or leave blank for all bodies)" << line_end;
out_stream << ")" << line_end;
for(i=0; i < bodies.size(); i++) {
    out_stream << bodies[i] << line_end;
}

            out_stream.close();
            out_stream.clear();

            return CALC_OK;
}
