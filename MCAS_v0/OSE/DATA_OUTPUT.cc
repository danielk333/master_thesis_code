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
        out_pos << std::setprecision(10) << std::scientific;
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
    temp_file << std::setprecision(10) << std::scientific;
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

int print_mercury6_big(std::string mercury6_big, std::vector<std::vector<double> > *temp_q, std::vector<std::vector<double> > *temp_v, std::vector<double> *temp_m, double date) {
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
            out_stream << " epoch (in days) = " << date << line_end;
            out_stream << ")---------------------------------------------------------------------" << line_end;
            unsigned int i = 1;
            out_stream << " " << "MERCURY m=" << (*temp_m)[i]/Msol << " r=20.D0 d=5.43" << line_end;
            out_stream << " " << (*temp_q)[i][0] << " " << (*temp_q)[i][1] << " " << (*temp_q)[i][2] << line_end;
            out_stream << " " << (*temp_v)[i][0]/AU*86400 << " " << (*temp_v)[i][1]/AU*86400 << " " << (*temp_v)[i][2]/AU*86400 <<  line_end;          
            out_stream << " " << "0. 0. 0." << line_end; i++;
            out_stream << " " << "VENUS   m=" << (*temp_m)[i]/Msol << " r=20.d0 d=5.24" << line_end;
            out_stream << " " << (*temp_q)[i][0] << " " << (*temp_q)[i][1] << " " << (*temp_q)[i][2] <<  line_end;
            out_stream << " " << (*temp_v)[i][0]/AU*86400 << " " << (*temp_v)[i][1]/AU*86400 << " " << (*temp_v)[i][2]/AU*86400 <<  line_end;          
            out_stream << " " << "0. 0. 0." << line_end; i++;
            out_stream << " " << "EARTH   m=" << (*temp_m)[i]/Msol << " r=20.d0 d=5.52" << line_end;
            out_stream << " " << (*temp_q)[i][0] << " " << (*temp_q)[i][1] << " " << (*temp_q)[i][2] <<  line_end;
            out_stream << " " << (*temp_v)[i][0]/AU*86400 << " " << (*temp_v)[i][1]/AU*86400 << " " << (*temp_v)[i][2]/AU*86400 <<  line_end;          
            out_stream << " " << "0. 0. 0." << line_end; i++;
            out_stream << " " << "MARS    m=" << (*temp_m)[i]/Msol << " r=20.d0 d=3.94" << line_end;
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


int print_mercury6_param(std::string mercury6_param, OPTIONS *opt) {
            std::ofstream out_stream;
            out_stream.open(mercury6_param.c_str(), std::ios::app);
            std::string line_end = "\n";

            if(out_stream.fail()) {
                return mercury6_param_DATA_SAVE_FAILED;
            }

            out_stream << std::uppercase << std::setprecision(17) << std::scientific;

            std::string integ;

            switch((*opt).integrator) {
                case 0:
                    integ = "hyb";
                    break;
                case 1:
                    integ = "BS";
                    break;
                default:
                    integ = "hyb";
            }

            double days_start = (*opt).integration_start;
            double int_time_days = ((*opt).integration_time)*365.25;
            double days_end = days_start + int_time_days;
            double time_step_days = ((*opt).integration_step)*365.25;
            double out_time;

            if(int_time_days/time_step_days >= 5e5) {
                out_time = 10*time_step_days;
            }
            else {
                out_time = time_step_days;   
            }

            out_stream << ")O+_06 Integration parameters  (WARNING: Do not delete this line!!)" << line_end;
            out_stream << ") Lines beginning with `)' are ignored." << line_end;
            out_stream << ")---------------------------------------------------------------------" << line_end;
            out_stream << ") Important integration parameters:" << line_end;
            out_stream << ")---------------------------------------------------------------------" << line_end;
            out_stream << "algorithm (MVS, BS, BS2, RADAU, HYBRID etc) = " << integ << line_end;
            out_stream << "start time (days)= " << days_start << line_end;
            out_stream << "stop time (days) = " << days_end << line_end;
            out_stream << "output interval (days) = " << out_time << line_end;
            out_stream << "timestep (days) = " << time_step_days << line_end;
            out_stream << "accuracy parameter=1.d-12" << line_end;
            out_stream << ")---------------------------------------------------------------------" << line_end;
            out_stream << ") Integration options:" << line_end;
            out_stream << ")---------------------------------------------------------------------" << line_end;
            out_stream << "stop integration after a close encounter = no" << line_end;
            out_stream << "allow collisions to occur = no" << line_end;
            out_stream << "include collisional fragmentation = no" << line_end;
            out_stream << "express time in days or years = years" << line_end;
            out_stream << "express time relative to integration start time = no" << line_end;
            out_stream << "output precision = medium" << line_end;
            out_stream << "< not used at present >" << line_end;
            out_stream << "include relativity in integration= no" << line_end;
            out_stream << "include user-defined force = no" << line_end;
            out_stream << ")---------------------------------------------------------------------" << line_end;
            out_stream << ") These parameters do not need to be adjusted often:" << line_end;
            out_stream << ")---------------------------------------------------------------------" << line_end;
            out_stream << "ejection distance (AU)= 500" << line_end;
            out_stream << "radius of central body (AU) = 0.005" << line_end;
            out_stream << "central mass (solar) = 1.0" << line_end;
            out_stream << "central J2 = 0" << line_end;
            out_stream << "central J4 = 0" << line_end;
            out_stream << " central J6 = 0" << line_end;
            out_stream << "< not used at present >" << line_end;
            out_stream << "< not used at present >" << line_end;
            out_stream << "Hybrid integrator changeover (Hill radii) = 3." << line_end;
            out_stream << "number of timesteps between data dumps = 500" << line_end;
            out_stream << "number of timesteps between periodic effects = 100" << line_end;

            out_stream.close();
            out_stream.clear();

            return CALC_OK;
}


