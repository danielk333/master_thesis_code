/*
 * main.cc
 *
 *  Created on: Jun 17, 2015
 *      Author: dankas
 */

//MCAS Monte Carlo Association Statistics

// Library includes
#include <string>
#include <math.h>
#include <algorithm>
 #include <sys/stat.h>

#include "boost_GS.h"
 
//Define global variables and physical constants
#include "define.hh"

//Functions include
#include "functions.hh"

//Make sim options global
OPTIONS opt;

int main (int argc, char *argv[]) {

/*######################################################################## */
 // DECLARE VARIABLES
/*######################################################################## */

    /* STATUS VARIABLES */
    int mercury6_status,JPL_status;
    int dummy_var;
    /* ITTERATORS */ 
    int i;
    unsigned int ui,uj;

    // Check for linux
    #ifdef __linux__

    #endif /* linux */

    //int k1;

    /* RANDOM SEED */
    srand (time(NULL));

    /* EMPTY_VECTS */
    std::vector<double> EMPTY_DOUBLE_VEC;
    std::vector<std::vector<double> > EMPTY_DOUBLE_MAT;

    /* STRINGS, CHARS AND FILENAMES */
    std::string temp_file_path;
    std::string temp_function_commands;
    std::string execute_command;
    std::string num;
    std::string source;
    std::string destination;
    int TEST_INT;
    std::string LF = "\n";
    std::string CR = "\r";
    std::string CRLF = CR + LF;

    std::string out_file_name;
    std::ofstream out_stream;
    std::string line_end = LF;

    //Streams for copy of file
    std::ifstream  src; 
    std::ofstream  dst;

    //Logfile out
    std::ofstream  log_out;

    //MURMHED 
    std::vector<OBSERVATION_record> objs;
    std::vector<OBSERVATION_record> real_objs;
    std::vector<OBSERVATION_record> real_objs_batch;
    std::vector<OBSERVATION_record> real_objs_batch2;
    OBSERVATION_record OBS_temp;

    OBSERVATION_record PB_temp;
    OBSERVATION_record TP_temp;

    std::string settings_folder = "SETTINGS/";

    /* Generate programs self_path for for dynamics folders */
std::string selfpath = get_selfpath();
selfpath.erase(selfpath.end()-3,selfpath.end());

    if(argc >= 3) {
        INPUT_db_file = argv[2];
    }
    else {
        INPUT_db_file = selfpath + settings_folder + "INPUT_OBJ.txt";
    }

    std::string OSE_config;
    if(argc >= 2) {
        OSE_config = argv[1];
    }
    else {
        OSE_config = selfpath + settings_folder + "OSE_config.cfg";
    }

    std::string mercury6_folder = "mercury6_OSE/";
    std::string mercury6_program = "./mercury6";
    std::string close6_program = "./close6";
    std::string element6_program = "./element6";
    std::string JPL_program = "JPL/./JPL";
    std::string JPL_folder = "JPL/";

    std::string TEMP_file = selfpath + "OUTPUT/TEMP_FILE_OSE.tmp";

    std::string pos_data_file = "pos.data";
    std::string vel_data_file = "vel.data";
    std::string mass_data_file = "mass.data";
    std::string spin_data_file = "spin.data";

    std::string mercury6_big = selfpath + mercury6_folder + "big.in";
    std::string mercury6_small = selfpath + mercury6_folder + "small.in";

    std::string mercury6_xv = selfpath + mercury6_folder + "xv.out";
    std::string mercury6_ce = selfpath + mercury6_folder + "ce.out";
    std::string mercury6_info = selfpath + mercury6_folder + "info.out";
    std::string mercury6_param = selfpath + mercury6_folder + "param.in";

    std::string mercury6_big_dmp = selfpath + mercury6_folder + "big.dmp";
    std::string mercury6_small_dmp = selfpath + mercury6_folder + "small.dmp";
    std::string mercury6_param_dmp = selfpath + mercury6_folder + "param.dmp";
    std::string mercury6_restart_dmp = selfpath + mercury6_folder + "restart.dmp";

    std::string JPL_pos = selfpath + JPL_OUT_folder + pos_data_file;
    std::string JPL_vel = selfpath + JPL_OUT_folder + vel_data_file;

    std::string mercury6_earth_clo = selfpath + mercury6_folder + "EARTH.clo";
    std::string mercury6_earth_aei = selfpath + mercury6_folder + "EARTH.aei";
    std::string mercury6_ORIG_aei = selfpath + mercury6_folder + "ORIG.aei";
    std::string mercury6_element_cfg = selfpath + mercury6_folder + "element.in";
    std::string mercury6_temp_obj_str;

    std::string mass_data_file = "mass.data";
    std::string CMS_m_init = selfpath + settings_folder + mass_data_file;

    std::string log_file = selfpath + "OUTPUT/outlog_OSE.txt";

    //###########
    //CREATE LIST OF ALL CREATED FILES
    // - FOR OPENING, CLOSING, COPY,
    // - DELETE, CREATE, and CLEANUP
    //###########
    std::vector<std::string> list_created_files;
    std::vector<std::vector<std::string> > matrix_created_files;

    list_created_files.push_back(JPL_pos);
    list_created_files.push_back(JPL_vel);
    list_created_files.push_back(mercury6_big);
    matrix_created_files.push_back(list_created_files); //persistant_init_files 0
    list_created_files.clear();

    list_created_files.push_back(mercury6_small);
    matrix_created_files.push_back(list_created_files); //massless_init_files 1
    list_created_files.clear();


    list_created_files.push_back(mercury6_xv);
    list_created_files.push_back(mercury6_ce);
    list_created_files.push_back(mercury6_info);
    list_created_files.push_back(mercury6_earth_clo);
    list_created_files.push_back(mercury6_earth_aei);
    list_created_files.push_back(mercury6_PB_aei);
    matrix_created_files.push_back(list_created_files); //output_files 2
    list_created_files.clear();

    list_created_files.push_back(mercury6_big_dmp);
    list_created_files.push_back(mercury6_small_dmp);
    list_created_files.push_back(mercury6_param_dmp);
    list_created_files.push_back(mercury6_restart_dmp);
    list_created_files.push_back(TEMP_file);
    matrix_created_files.push_back(list_created_files); //temp_files 3
    list_created_files.clear();

    
    list_created_files.push_back(log_file);
    matrix_created_files.push_back(list_created_files); //Output restults 4
    list_created_files.clear();

    //###########
    //CLEAN UP:
    // - IN CASE LAST RUN DIDNT EXIT PROPERLY
    // - TO MAKE ROOM FOR NEW RESULTS
    //###########
    clean_up(matrix_created_files,-1);

    create_file(TEMP_file);

    /* LOAD OPTIONS */
    std::string config_path = selfpath + OSE_config;
    error(load_file(&opt,config_path,DATA_OPTIONS),matrix_created_files);

    // INIT OUTPUT
    log_out.open(log_file.c_str(), std::ios::app);

    output_stream out(std::cout, log_out);
    out.type = opt.logfile;

    out << "Output type " << opt.logfile << " selected."; out.endl(); out.endl();


    //###########
    //TEMPORARY VARIABLES
    //###########

  //double r_vec_kep_state[3], v_vec_kep_state[3], mu_kep_state;
  //double p_kep_state, a_kep_state, ecc_kep_state, incl_kep_state, omega_kep_state, argp_kep_state, nu_kep_state, m_kep_state, arglat_kep_state, truelon_kep_state, lonper_kep_state;

    std::vector<std::vector<double> > time_corrected_vectors_obj,time_corrected_vectors_earth;
    std::vector<double> temp_kepler,temp_state_vec,x_temp,v_temp;
unsigned int propagation_counter;
    unsigned int TEMP_ID_index;
    std::vector<double> temp_obj_pos_a;
    std::vector<double> temp_obj_pos_b;

    std::vector<double> temp_obj_vel_a;
    std::vector<double> temp_obj_vel_b;

    std::vector<std::pair<size_t, double_vector_iterator> > order_D_c_vec;
    //size_t D_c_n;

    std::vector<std::vector<double> > temp_q;
    std::vector<std::vector<double> > temp_v;
    std::vector<std::vector<double> > temp_t;
    std::vector<std::vector<double> > temp_sim;
    std::vector<double> temp_m;

    std::vector<double> temp_id_vec;
    
    std::vector<double> temp_row;
    std::vector<double> temp_data;

    std::vector<std::vector<double> > temp_obj_data;
    std::vector<std::vector<double> > temp_earth_data;

    std::vector<std::vector<double> > temp_mat;
    std::vector<std::vector<double> > temp_mat2;
    std::vector<double> temp_vec;

std::vector<std::vector<unsigned int> > index_matrix;
std::vector<std::vector<unsigned int> > cluster_matrix_next_db;
std::string buffer_temp;

std::vector<unsigned int> EMPTY_UINT_VEC;

    std::vector<double> current_body_xv;

    bool TRUTH_CHECKER;

    std::string temp_batch_name;
    std::string temp_system_command;

    std::vector<std::string> list_of_batches_d_matrix;

    //###########
    //RESULT STORAGE
    //###########
    std::vector<std::vector<double> > encounter_data;
    std::vector<std::vector<double> > earth_encounter_data;
    std::vector<std::vector<double> > encounter_stat_data;

    std::vector<std::vector<double> > filtered_encounter_data;
    std::vector<std::vector<double> > filtered_encounter_data_save;
    std::vector<unsigned int> id_vec;
    std::vector<unsigned int> orig_id_vec;

    unsigned int first_id;
    std::vector<std::vector<double> > body_family;

    //###########
    //Counters
    //###########

    //###########
    //MISC
    //###########
    std::vector<double>::iterator TEMP_ID;

    //###########
    //Physical params
    //###########

    //./rdeph binEphem.200 2457196.70524 (planet) PosVel 0
    execute_command = JPL_program;
    execute_command.insert(0,selfpath);
    execute_command = execute_command + " " + selfpath + JPL_folder + "binEphem.200 2457196.70524 10 PosVel 0 " + JPL_OUT_folder;
    out << "Calling " << execute_command; out.endl();
    JPL_status = system(execute_command.c_str());
    out << "JPL exit status: " << JPL_status; out.endl();

    num = "0";
    for(i=0; i < 8; i++) {
        execute_command = JPL_program;
        execute_command.insert(0,selfpath);
        execute_command = execute_command  + " " + selfpath + JPL_folder + "binEphem.200 2457196.70524 " + num + " PosVel 0 " + JPL_OUT_folder;
        out << "Calling " << execute_command; out.endl();
        JPL_status = system(execute_command.c_str());
        out << "JPL exit status: " << JPL_status; out.endl();

        num[0]++;
    }

    std::vector<double> current_body_kep;
    current_body_kep.push_back(0);current_body_kep.push_back(0);current_body_kep.push_back(0);current_body_kep.push_back(0);current_body_kep.push_back(0);current_body_kep.push_back(0);
    std::vector<double> current_body_char;
    current_body_char.push_back(0);current_body_char.push_back(0);current_body_char.push_back(0);current_body_char.push_back(0);

    out << "##### USING INTEGRATOR " << opt.integrator << " #####"; out.endl();

            temp_q.clear(); temp_v.clear(); temp_m.clear();
            error(load_file(&temp_q,JPL_pos,DATA_MATRIX),matrix_created_files);
            error(load_file(&temp_v,JPL_vel,DATA_MATRIX),matrix_created_files);
            error(load_file(&temp_m,CMS_m_init,DATA_VECTOR),matrix_created_files);

            heliocentric_qv(&temp_q,&temp_v);
            rot_to_invariable_plane(&temp_q,&temp_v,&temp_m);

            out << "# Printing data to " << mercury6_big; out.endl();
            error(print_mercury6_big(mercury6_big,&temp_q,&temp_v,&temp_m,opt.integration_start),matrix_created_files);

            out << "# Printing data to " << mercury6_param; out.endl();
            error(print_mercury6_param(mercury6_param, &opt),matrix_created_files);

            out_stream.open(mercury6_small.c_str(), std::ios::app);
            out_stream << std::uppercase << std::setprecision(17) << std::scientific;


            current_body_xv.clear();
            current_body_kep.clear();

            // Format: 0a 1e 2i 3omega 4Omega 5nu 6m
            error(load_data_matrix(&current_body_kep,INPUT_db_file),matrix_created_files);

            current_body_xv = kepler_to_xv(current_body_kep, G*(current_body_kep[6]+Msol));
            
            temp_q.clear(); temp_v.clear();
            

            out << "# Building Small-body initial data: " << mercury6_small; out.endl();

            out_stream << ")O+_06 Small-body initial data  (WARNING: Do not delete this line!!)" << line_end;
            out_stream << ") Lines beginning with `)' are ignored." << line_end;
            out_stream << ")---------------------------------------------------------------------" << line_end;
            out_stream << " style (Cartesian, Asteroidal, Cometary) = Car" << line_end;
            out_stream << ")---------------------------------------------------------------------" << line_end;
            
            out_stream << " " << "ORIG" << "      ep=" << opt.integration_start << line_end;
            out_stream << " " << current_body_xv[0]/AU << " " << current_body_xv[1]/AU << " " << current_body_xv[2]/AU;
            out_stream << " " << current_body_xv[3]/AU*86400 << " " << current_body_xv[4]/AU*86400 << " " << current_body_xv[5]/AU*86400;
            out_stream << " " << "0. 0. 0." << line_end;

            out_stream.close();
            out_stream.clear();
            temp_q.clear(); temp_v.clear();

            out << "########################################################################"; out.endl();
            out << "# Starting particle propagation"; out.endl();
            out << "########################################################################"; out.endl();

            execute_command = "cd " + selfpath + mercury6_folder + "; " + mercury6_program;
            out << "# Calling " << execute_command; out.endl();
            mercury6_status = system(execute_command.c_str());
            propagation_counter++;

            out << "- mercury6 exit status: " << mercury6_status; out.endl();
            error(mercury6_status,matrix_created_files);


//Remove config file for rewriting
if(file_exists(mercury6_element_cfg)) {
    remove(mercury6_element_cfg.c_str());
    out << "Removed file: " << mercury6_element_cfg; out.endl();
}

out << "# Building Small-body element extract data file: " << mercury6_element_cfg; out.endl();
out_stream.open(mercury6_element_cfg.c_str(), std::ios::app);

//Print more exact data on object that had a close encounter
out_stream << ")O+_06 element  (WARNING: Do not delete this line!!)" << line_end;
out_stream << ") Lines beginning with `)' are ignored." << line_end;
out_stream << ")---------------------------------------------------------------------" << line_end;
out_stream << " number of input files = 1" << line_end;
out_stream << ")---------------------------------------------------------------------" << line_end;
out_stream << ") List the input files, one per line" << line_end;
out_stream << " xv.out" << line_end;
out_stream << ")---------------------------------------------------------------------" << line_end;
out_stream << " type of elements (central body, barycentric, Jacobi) = Cen" << line_end;
out_stream << " minimum interval between outputs (days) = 9.2d0" << line_end;
out_stream << " express time in days or years = years" << line_end;
out_stream << " express time relative to integration start time = yes" << line_end;
out_stream << ")---------------------------------------------------------------------" << line_end;
out_stream << ") Output format? (e.g. a8.4 => semi-major axis with 8 digits & 4 dec. places)" << line_end;
//out_stream << " a8.5 e8.6 i8.4 g8.4 n8.4 l8.4 x13e y13e z13e vx13e vy13e vz13e " << line_end;
out_stream << " a8.5 e8.6 i8.4 g8.4 n8.4 l8.4 f8.4 m13e " << line_end;
out_stream << ")---------------------------------------------------------------------" << line_end;
out_stream << ") Which bodies do you want? (List one per line or leave blank for all bodies)" << line_end;
out_stream << ")" << line_end;
out_stream << "ORIG" << line_end;

out_stream.close();
out_stream.clear();

out << "Created file: " << mercury6_element_cfg; out.endl();

execute_command = "cd " + selfpath + mercury6_folder + "; " + element6_program;
out << "Calling " << execute_command; out.endl();
mercury6_status = system(execute_command.c_str());
out << "element6 exit status: " << mercury6_status; out.endl();
error(mercury6_status,matrix_created_files);

PB_time_series_data.clear();
out << "Attempting to load: " << mercury6_ORIG_aei; out.endl();
error(read_mercury6_matrix(&PB_time_series_data,mercury6_ORIG_aei),matrix_created_files);

w.clear();
w_save.clear();
for(i=0; i < 6; i++) {
    w.push_back(1);
}



for(i=0; i < PB_time_series_data.size(); i++) {


}





clean_up(matrix_created_files,1);
clean_up(matrix_created_files,2);
clean_up(matrix_created_files,3);

    return CALC_OK;
}


void error(int SIM_STATUS, std::vector<std::vector<std::string> > matrix_created_files) {

  if(SIM_STATUS != CALC_OK) {
    std::cout << std::endl << "ERROR ENCOUTERD: " << SIM_STATUS << std::endl << std::endl;
    
    //CLEANUP
    //clean_up(matrix_created_files);

    exit(SIM_STATUS);
  }
}

void clean_up(std::vector<std::vector<std::string> > matrix_created_files, int type) {
    unsigned int i,j;
    bool selector;
    if(type == -1) {
        selector = TRUE;
    }
    else {
        selector = FALSE;
    }
    for(j=0; j < matrix_created_files.size(); j++ ) {
        if((type == j) || selector) {
            for (i = 0; i < matrix_created_files[j].size(); ++i) {
                if(file_exists(matrix_created_files[j][i])) {
                    remove(matrix_created_files[j][i].c_str());
                    std::cout << "Removed file: " << matrix_created_files[j][i] << std::endl;
                }
                else if(folder_exists(matrix_created_files[j][i])) {
                    remove_folder(matrix_created_files[j][i]);
                    std::cout << "Removed file: " << matrix_created_files[j][i] << std::endl;
                }
            }
        }
    }
}

