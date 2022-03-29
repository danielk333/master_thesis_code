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

 
//Define global variables and physical constants
#include "define.hh"

//Functions include
#include "functions.hh"

//Make sim options global
OPTIONS opt;

int main(void) {

/*######################################################################## */
 // DECLARE VARIABLES
/*######################################################################## */

    /* STATUS VARIABLES */
    int PBE_status,CMS_status,mercury6_status,JPL_status,OrbClone_status;
    int dummy_var;
    /* ITTERATORS */ 
    int i,j,k;

    int k1;

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

    //FUNCTIONS TO USE
    std::vector<std::string> FUNCTIONS;
    FUNCTIONS.push_back("SH");
    FUNCTIONS.push_back("D");
    FUNCTIONS.push_back("rho2");
    FUNCTIONS.push_back("varrho1");

    /* Generate programs self_path for for dynamics folders */
std::string selfpath = get_selfpath();
selfpath.erase(selfpath.end()-8,selfpath.end());

    /*  */
#ifdef __linux__
    std::string CMS_program = "CMS/./CMS_0_4";
    std::string mercury6_folder = "mercury6/";
    std::string statistics_folder = "STATISTICS_OUTPUT/";
    std::string mercury6_program = "./mercury6";
    std::string close6_program = "./close6";
    std::string element6_program = "./element6";
    std::string PBE_program = "PBE/./PBE_0_1";
    std::string JPL_program = "JPL/./JPL_0_1";
    std::string settings_folder = "SETTINGS/";
    std::string BOTKKE_folder = "BOTKKE/";
    std::string SS_model_folder = "SS_model/";
    std::string MURMHED_db_folder = "MURMHED/";
    std::string MURMHED_db_file = MURMHED_db_folder + "MURMHED.txt";
    std::string JPL_folder = "JPL/";
    std::string CMS_init_folder = "CMS/INIT/MASSIVE_BODIES/";
    std::string JPL_OUT_folder = "PBE/INIT/MASSIVE_BODIES/";
    std::string PBE_IC_PRODUCER_folder = "PBE/INIT/PRODUCER/";
    std::string PBE_OUT_folder = "PBE/OUT_DATA/";
    std::string CMS_MASSLESS_folder = "CMS/INIT/MASSLESS_BODIES/";
    std::string CMS_OUT_folder = "CMS/OUT_DATA/";
    std::string debug_folder = "DEBUG_DATA/";
    std::string OrbClone_folder = "OrbClone/";
    std::string OrbClone_program = "./OrbClone_0_1";
#endif /* linux */

    std::string TEMP_file = selfpath + "OUTPUT/TEMP_FILE.tmp";

    std::vector<std::string> MURMHED_D_mat_files;
    for(i=0; i < FUNCTIONS.size(); i++) {
        MURMHED_D_mat_files.push_back(selfpath + MURMHED_db_folder + "MURMHED_" + FUNCTIONS[i]);
    }

    std::vector<std::string> MURMHED_profile_files;
    for(i=0; i < FUNCTIONS.size(); i++) {
        MURMHED_profile_files.push_back(selfpath + statistics_folder + "MURMHED_" + FUNCTIONS[i] + "_profile.txt");
    }

    std::string MCAS_config = "MCAS_config.cfg";
    std::string PBE_IC_body = "body_data.data";
    std::string PBE_config = "PBE_config.cfg";
    std::string CMS_config = "CMS_config.cfg";

    std::string pos_data_file = "pos.data";
    std::string vel_data_file = "vel.data";
    std::string mass_data_file = "mass.data";
    std::string spin_data_file = "spin.data";

    std::string mercury6_big = selfpath + mercury6_folder + "big.in";
    std::string mercury6_small = selfpath + mercury6_folder + "small.in";

    std::string mercury6_xv = selfpath + mercury6_folder + "xv.out";
    std::string mercury6_ce = selfpath + mercury6_folder + "ce.out";
    std::string mercury6_info = selfpath + mercury6_folder + "info.out";

    std::string mercury6_big_dmp = selfpath + mercury6_folder + "big.dmp";
    std::string mercury6_small_dmp = selfpath + mercury6_folder + "small.dmp";
    std::string mercury6_param_dmp = selfpath + mercury6_folder + "param.dmp";
    std::string mercury6_restart_dmp = selfpath + mercury6_folder + "restart.dmp";

    std::string JPL_PBE_pos = selfpath + JPL_OUT_folder + pos_data_file;
    std::string JPL_PBE_vel = selfpath + JPL_OUT_folder + vel_data_file;
    std::string JPL_CMS_pos = selfpath + CMS_init_folder + pos_data_file;
    std::string JPL_CMS_vel = selfpath + CMS_init_folder + vel_data_file;

    std::string PBE_q_out = selfpath + PBE_OUT_folder + "particle_q_data.txt";
    std::string PBE_v_out = selfpath + PBE_OUT_folder + "particle_v_data.txt";
    std::string PBE_m_out = selfpath + PBE_OUT_folder + "particle_m_data.txt";
    std::string PBE_t_out = selfpath + PBE_OUT_folder + "particle_t_data.txt";
    std::string PBE_abs_v_out = selfpath + PBE_OUT_folder + "particle_abs_v_data.txt";
    std::string PBE_kep_out = selfpath + PBE_OUT_folder + "particle_kep_data.txt";
    std::string PBE_body_kep_out = selfpath + PBE_OUT_folder + "body_kep_data.txt";
    std::string PBE_M_out = selfpath + PBE_OUT_folder + "body_m_data.txt";
    std::string PBE_sim_out = selfpath + PBE_OUT_folder + "sim_data.txt";

    std::string PBE_weight_statistics = selfpath + PBE_OUT_folder + "particle_weight_stat.txt";
    std::string PBE_abs_v_statistics = selfpath + PBE_OUT_folder + "particle_ejection_v_stat.txt";
    std::string PBE_M_statistics = selfpath + PBE_OUT_folder + "PBE_M_stat.txt";

    std::string PBE_m_dist_in = selfpath + PBE_IC_PRODUCER_folder + "mass_dist.data";

    std::string CMS_q_in = selfpath + CMS_MASSLESS_folder + pos_data_file;
    std::string CMS_v_in = selfpath + CMS_MASSLESS_folder + vel_data_file;
    std::string CMS_m_in = selfpath + CMS_MASSLESS_folder + mass_data_file;

    std::string CMS_m_init = selfpath + CMS_init_folder + mass_data_file;

    std::string CMS_data_out = selfpath + CMS_OUT_folder + "body_data.txt";
    std::string CMS_diag_out = selfpath + CMS_OUT_folder + "energy_diag.txt";
    std::string CMS_sett_out = selfpath + CMS_OUT_folder + "sim_settings.txt";

    std::string mercury6_earth_clo = selfpath + mercury6_folder + "EARTH.clo";
    std::string mercury6_earth_aei = selfpath + mercury6_folder + "EARTH.aei";
    std::string mercury6_element_cfg = selfpath + mercury6_folder + "element.in";
    std::string mercury6_temp_obj_str;

    std::string met_encounter_data_out = selfpath + statistics_folder + "met_encounter_data.txt";
    std::string earth_encounter_data_out = selfpath + statistics_folder + "earth_encounter_data.txt";
    std::string met_before_encounter_data_out = selfpath + statistics_folder + "met_before_encounter_data.txt";
    std::string statistics_PB_kep_out = selfpath + statistics_folder + "ejector_kep_data.txt";
    std::string statistics_PB_char_out = selfpath + statistics_folder + "ejector_char_data.txt";
    std::string statistics_PB_n_out = selfpath + statistics_folder + "ejector_n_data.txt";


    std::vector<std::string> statistics_association_profile;
    std::vector<std::string> statistics_Dc_min_err;
    for(i=0; i < FUNCTIONS.size(); i++) {
        statistics_association_profile.push_back(selfpath + statistics_folder + "association_profile_" + FUNCTIONS[i] + ".txt");
        statistics_Dc_min_err.push_back(selfpath + statistics_folder + "min_error_" + FUNCTIONS[i] + "_Dc.txt");
    }

    std::string dist_debug_a =  selfpath + debug_folder + "a_dist.txt";
    std::string dist_debug_e =  selfpath + debug_folder + "e_dist.txt";
    std::string dist_debug_i =  selfpath + debug_folder + "i_dist.txt";
    std::string dist_debug_omega =  selfpath + debug_folder + "omega_dist.txt";
    std::string dist_debug_Omega =  selfpath + debug_folder + "Omega_dist.txt";

    std::string OrbClone_file =  selfpath + OrbClone_folder + "orb_clones.txt";

    std::string log_file = selfpath + "OUTPUT/outlog.txt";

    char data_out_folder[37] = "_______________________|_simulation/";
    char *time_char;

    time_t rawtime;
    struct tm *timeinfo;

    time(&rawtime);
    timeinfo = localtime(&rawtime);

    time_char = asctime(timeinfo);

    for(i = 0; i < 24; i++) {
        data_out_folder[i] = *(time_char + i);
    }

    data_out_folder[3] = '_'; data_out_folder[7] = '_'; data_out_folder[10] = '_'; data_out_folder[13] = '-'; data_out_folder[16] = '-'; data_out_folder[19] = '_';
    if(data_out_folder[8] = ' ') {
        data_out_folder[8] = '0';
    }
    std::string data_out_folder_string = data_out_folder;
    data_out_folder_string = selfpath + statistics_folder + data_out_folder_string;


    //###########
    //CREATE LIST OF ALL CREATED FILES
    // - FOR OPENING, CLOSING, COPY,
    // - DELETE, CREATE, and CLEANUP
    //###########
    std::vector<std::string> list_created_files;
    std::vector<std::vector<std::string> > matrix_created_files;

    list_created_files.push_back(JPL_PBE_pos);
    list_created_files.push_back(JPL_PBE_vel);
    list_created_files.push_back(JPL_CMS_pos);
    list_created_files.push_back(JPL_CMS_vel);
    list_created_files.push_back(mercury6_big);
    list_created_files.push_back(PBE_m_dist_in);
    list_created_files.push_back(OrbClone_file);
    matrix_created_files.push_back(list_created_files); //persistant_init_files 0
    list_created_files.clear();

    list_created_files.push_back(CMS_q_in);
    list_created_files.push_back(CMS_v_in);
    list_created_files.push_back(CMS_m_in);
    list_created_files.push_back(mercury6_small);
    matrix_created_files.push_back(list_created_files); //massless_init_files 1
    list_created_files.clear();

    list_created_files.push_back(PBE_q_out);
    list_created_files.push_back(PBE_v_out);
    list_created_files.push_back(PBE_m_out);
    list_created_files.push_back(PBE_t_out);
    list_created_files.push_back(PBE_M_out);
    list_created_files.push_back(PBE_sim_out);
    list_created_files.push_back(PBE_abs_v_out);
    list_created_files.push_back(PBE_kep_out);
    list_created_files.push_back(PBE_body_kep_out);
    list_created_files.push_back(mercury6_xv);
    list_created_files.push_back(mercury6_ce);
    list_created_files.push_back(mercury6_info);
    list_created_files.push_back(CMS_data_out);
    list_created_files.push_back(CMS_diag_out);
    list_created_files.push_back(CMS_sett_out);
    list_created_files.push_back(mercury6_earth_clo);
    list_created_files.push_back(mercury6_earth_aei);
    matrix_created_files.push_back(list_created_files); //output_files 2
    list_created_files.clear();

    list_created_files.push_back(mercury6_big_dmp);
    list_created_files.push_back(mercury6_small_dmp);
    list_created_files.push_back(mercury6_param_dmp);
    list_created_files.push_back(mercury6_restart_dmp);
    list_created_files.push_back(TEMP_file);
    matrix_created_files.push_back(list_created_files); //temp_files 3
    list_created_files.clear();

    list_created_files.push_back(met_encounter_data_out);
    list_created_files.push_back(earth_encounter_data_out);
    list_created_files.push_back(met_before_encounter_data_out);
    list_created_files.push_back(statistics_PB_kep_out);
    list_created_files.push_back(statistics_PB_char_out);
    list_created_files.push_back(statistics_PB_n_out);
    list_created_files.push_back(PBE_weight_statistics);
    list_created_files.push_back(PBE_abs_v_statistics);
    list_created_files.push_back(PBE_M_statistics);
    list_created_files.push_back(log_file);
    for(i=0; i < FUNCTIONS.size(); i++) {
        list_created_files.push_back(statistics_association_profile[i]);
        list_created_files.push_back(statistics_Dc_min_err[i]);
    }
    matrix_created_files.push_back(list_created_files); //Output restults 4
    list_created_files.clear();

    list_created_files.push_back(dist_debug_a);
    list_created_files.push_back(dist_debug_e);
    list_created_files.push_back(dist_debug_i);
    list_created_files.push_back(dist_debug_omega);
    list_created_files.push_back(dist_debug_Omega);
    matrix_created_files.push_back(list_created_files); //DEBUG DATA 5
    list_created_files.clear();

    //###########
    //CLEAN UP:
    // - IN CASE LAST RUN DIDNT EXIT PROPERLY
    // - TO MAKE ROOM FOR NEW RESULTS
    //###########
    clean_up(matrix_created_files,-1);

    create_file(TEMP_file);

    /* LOAD OPTIONS */
    std::string config_path = selfpath + settings_folder + MCAS_config;
    error(load_file(&opt,config_path,DATA_OPTIONS),matrix_created_files);

    // INIT OUTPUT
    log_out.open(log_file.c_str(), std::ios::app);

    output_stream out(std::cout, log_out);
    out.type = opt.logfile;

    out << "Output type " << opt.logfile << " selected."; out.endl(); out.endl();

    //Family data
    std::string temp_family;
    std::vector<std::string> family_data;
    std::vector<std::string> family_dist;
    for(i=0; i < 12; i++) {
        temp_family = "family" + IntToString(i) + ".data";
        family_data.push_back(temp_family);
    }
    family_dist.push_back("S0");
    family_dist.push_back("S1.01");
    family_dist.push_back("S1.02");
    family_dist.push_back("S1.03");
    family_dist.push_back("S1.04");
    family_dist.push_back("S1.05");
    family_dist.push_back("S1.06");
    family_dist.push_back("S1.07");
    family_dist.push_back("S1.08");
    family_dist.push_back("S1.09");
    family_dist.push_back("S1.10");
    family_dist.push_back("S2");
    family_dist.push_back("S.hilda");
    family_dist.push_back("Sc");
    family_dist.push_back("SH");
    family_dist.push_back("SL");
    family_dist.push_back("SR");
    family_dist.push_back("SS");
    family_dist.push_back("ST");
    family_dist.push_back("NEO_family.data");

    std::vector<std::vector<double> > orbital_dist_data;

    //###########
    //DATABASE STORAGE
    //###########

    std::vector<std::vector<double> > MURMHED_db;
    std::vector<std::vector<double> > MURMHED_db_distance_matrix_part;
    /*std::vector<std::vector<std::vector<double> > > MURMHED_db_distance_matrix_bundle;
    for(i=0; i < FUNCTIONS.size(); i++) {
        MURMHED_db_distance_matrix_bundle.push_back(EMPTY_DOUBLE_MAT);
    }*/


    //###########
    //SYSTEM VARIABLES
    //###########

    double mem_batch_estimate;
    double faction_of_avalible_mem;
    unsigned int calculation_batches;
    unsigned int items_a_batch;
    unsigned int residue_items;

    unsigned int residue_toggle_loop;
    bool residue_toggle;

    //###########
    //TEMPORARY VARIABLES
    //###########

  //double r_vec_kep_state[3], v_vec_kep_state[3], mu_kep_state;
  //double p_kep_state, a_kep_state, ecc_kep_state, incl_kep_state, omega_kep_state, argp_kep_state, nu_kep_state, m_kep_state, arglat_kep_state, truelon_kep_state, lonper_kep_state;

    std::vector<std::vector<double> > time_corrected_vectors_obj,time_corrected_vectors_earth;
    std::vector<double> temp_kepler,temp_state_vec,x_temp,v_temp;

    unsigned int TEMP_ID_index;
    std::vector<double> temp_obj_pos_a;
    std::vector<double> temp_obj_pos_b;

    std::vector<double> temp_obj_vel_a;
    std::vector<double> temp_obj_vel_b;

    std::vector<std::pair<size_t, double_vector_iterator> > order_D_c_vec;
    size_t D_c_n;

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

    unsigned int association_loop_counter;
    unsigned int association_loop_error_counter;
    unsigned int empty_bin_counter;
    //unsigned int association_loop_counter_max = 10000;
    //bool adaptive_step;
    double association_adaptive_step_threshold;
    double min_d_crit_step;
    double max_d_crit_step;

    double step_multiplyer_inc;

    unsigned int number_of_objects_db_load;

    unsigned int bin_member_counter;

    bool next_step_ok;
    unsigned int visited_zero_progress_counter;
    unsigned int visited_max_progress_counter;

unsigned int char_index_1,char_index_2;
std::vector<std::vector<unsigned int> > index_matrix;
std::vector<std::vector<unsigned int> > cluster_matrix_next_db;
std::string buffer_temp;

std::vector<unsigned int> EMPTY_UINT_VEC;

    double d_crit_bin_width;
    unsigned int d_crit_bin_n;

    double d_crit;
    double d_crit_old;
    double d_crit_step;
    //unsigned int d_crit_bins;
    std::vector<double> association_vec;
    std::vector<double> d_crit_vec;
    double fraction_associated,fraction_associated_old,fraction_associated_next;

    bool TRUTH_CHECKER;

    std::string temp_batch_name;
    std::string temp_system_command;

    std::vector<std::string> list_of_batches_d_matrix;

    //unsigned int vector_counter;
    //###########
    //MODEL STORAGE
    //###########
    std::vector<double> grun_model_mass_flux;

    //###########
    //RESULT STORAGE
    //###########
    std::vector<std::vector<double> > encounter_data;
    std::vector<std::vector<double> > earth_encounter_data;
    std::vector<std::vector<double> > encounter_stat_data;
std::vector<double> error_function;
    std::vector<std::vector<double> > D_matix;
    std::vector<std::vector<unsigned int> > cluster_matrix;
    std::vector<double> fractions;
    std::vector<std::vector<unsigned int> > cluster_matrix_next;
    std::vector<double> fractions_next;
    std::vector<std::vector<double> > association_profile;
    for(i=0; i < FUNCTIONS.size(); i++) {
        association_profile.push_back(EMPTY_DOUBLE_VEC); //dcrit
        association_profile.push_back(EMPTY_DOUBLE_VEC); //fractions
    }

    std::vector<std::vector<double> > association_profile_save;
    std::vector<std::vector<double> > association_profile_save_means;

    std::vector<std::vector<double> > association_profile_db;
    for(i=0; i < FUNCTIONS.size(); i++) {
        association_profile_db.push_back(EMPTY_DOUBLE_VEC); //dcrit
        association_profile_db.push_back(EMPTY_DOUBLE_VEC); //fractions
    }

    std::vector<std::vector<double> > association_profile_db_save;
    std::vector<std::vector<double> > association_profile_db_save_means;

    std::vector<std::vector<double> > filtered_encounter_data;
    std::vector<std::vector<double> > filtered_encounter_data_save;
    std::vector<unsigned int> id_vec;
    std::vector<unsigned int> orig_id_vec;

    unsigned int first_id;
    std::vector<std::vector<double> > body_family;

    //###########
    //Counters
    //###########
    unsigned int total_n_created, total_PB_created, total_met_collided, current_n_created, total_showers_created;

    //###########
    //MISC
    //###########
    bool orbit_ok;
    double body_mass;
    std::vector<double>::iterator TEMP_ID;

    //###########
    //Physical params
    //###########
    double FUDGE = 1.1;
    double earth_hill_r = FUDGE*1.00000261*(1-0.01671123)*cbrt(5.97219e24/(3*1.98855e30));

    /*######################################################################## */
    // LOAD DATA
    /*######################################################################## */

    temp_file_path = selfpath + settings_folder + family_data[opt.family];
    error(load_file(&body_family,temp_file_path,DATA_FAMILY),matrix_created_files);

    out << "PARENT BODY DISTRIBUTION: "; out.endl();
    out << body_family[0][0] << ": " << body_family[0][1] << " " <<  body_family[0][2]; out.endl();
    out << body_family[1][0] << ": " << body_family[1][1] << " " <<  body_family[1][2]; out.endl();
    out << body_family[2][0] << ": " << body_family[2][1] << " " <<  body_family[2][2]; out.endl();
    out << body_family[3][0] << ": " << body_family[3][1] << " " <<  body_family[3][2]; out.endl();

    if(opt.family == 0) {
        temp_file_path = selfpath + SS_model_folder + family_dist[opt.family];
        error(load_file(&orbital_dist_data,temp_file_path,DATA_ORBDIST_SS),matrix_created_files);
    }
    else if(opt.family == 1) {
        for(i=0; i < 9; i++) {
            temp_file_path = selfpath + SS_model_folder + family_dist[opt.family+i];
            error(load_file(&orbital_dist_data,temp_file_path,DATA_ORBDIST_SS),matrix_created_files);
        }
    }
    else if(opt.family > 1 && opt.family < 10) {
        temp_file_path = selfpath + SS_model_folder + family_dist[opt.family+9];
        error(load_file(&orbital_dist_data,temp_file_path,DATA_ORBDIST_SS),matrix_created_files);
    }
    else if(opt.family == 10) {
        temp_file_path = selfpath + BOTKKE_folder + family_dist[opt.family+9];
        error(load_file(&orbital_dist_data,temp_file_path,DATA_ORBDIST),matrix_created_files);
    }
    else if(opt.family == 11) {
        execute_command = selfpath + OrbClone_folder + OrbClone_program + " kep 500 kep";
        out << "Calling " << execute_command; out.endl();
        OrbClone_status = system(execute_command.c_str());
        out << "OrbClone exit status: " << OrbClone_status; out.endl();

        error(load_file(&orbital_dist_data,OrbClone_file,DATA_MATRIX),matrix_created_files);
    }

    /*######################################################################## */
    // CREATE APPROXIAMTE DUST MASS SIZE DISTRIBUTION
    /*######################################################################## */
    //GRUN: 1e-18g to 1g dust
    double min_mass_detectable = 1e-7;
    unsigned int grun_model_resolution = 50;
    grun_model_mass_flux = grun_flux(min_mass_detectable,grun_model_resolution);

    //grun model particles/m^2/year, but year and area does not matter, distribution is important
    // so probabillity is vector devided by number of particles
    
    grun_model_mass_flux = grun_model_mass_flux/(sum_v(grun_model_mass_flux));

    temp_mat2.clear();
    temp_mat2.push_back(grun_model_mass_flux);

    //print_vector(grun_model_mass_flux);

    save_mat(PBE_m_dist_in, &temp_mat2);
    /*######################################################################## */
    // CALCULATE ORBITAL PDF's
    /*######################################################################## */
    statistical_set semi_axis;
    if(opt.family >= 0 && opt.family <= 10) {
        semi_axis.data = (extract_col(orbital_dist_data,0)/((double)1 - extract_col(orbital_dist_data,1)));
    }
    else {
        semi_axis.data = (extract_col(orbital_dist_data,0));
    }

    if(opt.family >= 0 && opt.family <= 9) {
        orbital_dist_data = swap_cols(orbital_dist_data,3,4);
    }
    statistical_set ecc(extract_col(orbital_dist_data,1));
    statistical_set inc(extract_col(orbital_dist_data,2));
    statistical_set omega(extract_col(orbital_dist_data,3));
    statistical_set Omega(extract_col(orbital_dist_data,4));
    
    semi_axis.convert_to_PDF_and_clear(0);
    ecc.convert_to_PDF_and_clear(0);
    inc.convert_to_PDF_and_clear(0);
    omega.convert_to_PDF_and_clear(0);
    Omega.convert_to_PDF_and_clear(0);

    semi_axis.save_histogram_to_file(dist_debug_a);
    ecc.save_histogram_to_file(dist_debug_e);
    inc.save_histogram_to_file(dist_debug_i);
    omega.save_histogram_to_file(dist_debug_omega);
    Omega.save_histogram_to_file(dist_debug_Omega);

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
    switch(opt.integrator) {
        case 0:
            source = selfpath + JPL_OUT_folder + pos_data_file;
            destination = selfpath + CMS_init_folder + pos_data_file;
            src.open(source.c_str(), std::ios::binary);
            dst.open(destination.c_str(), std::ios::binary);
            dst << src.rdbuf();
            src.close(); src.clear();
            dst.close(); dst.clear();
 
            source = selfpath + JPL_OUT_folder + vel_data_file;
            destination = selfpath + CMS_init_folder + vel_data_file;
            src.open(source.c_str(), std::ios::binary);
            dst.open(destination.c_str(), std::ios::binary);
            dst << src.rdbuf();
            src.close(); src.clear();
            dst.close(); dst.clear();
            
            break;
        case 1:
            temp_q.clear(); temp_v.clear(); temp_m.clear();
            error(load_file(&temp_q,JPL_PBE_pos,DATA_MATRIX),matrix_created_files);
            error(load_file(&temp_v,JPL_PBE_vel,DATA_MATRIX),matrix_created_files);
            error(load_file(&temp_m,CMS_m_init,DATA_VECTOR),matrix_created_files);

            heliocentric_qv(&temp_q,&temp_v);
            rot_to_invariable_plane(&temp_q,&temp_v,&temp_m);

            error(print_mercury6_big(mercury6_big,&temp_q,&temp_v,&temp_m),matrix_created_files);
            break;
    }

    /*######################################################################## */
    // LOAD MURMHED DATA 
    /*######################################################################## */

out << "###################################"; out.endl();
out << "###### CHECKING MURMHED DATA ######"; out.endl();
out << "###################################"; out.endl();


bool DO_THEY_EXIST;
DO_THEY_EXIST = TRUE;

for(i=0; i < MURMHED_D_mat_files.size(); i++) {
    DO_THEY_EXIST = DO_THEY_EXIST && folder_exists(MURMHED_D_mat_files[i]);
}

if(DO_THEY_EXIST) {
    out << "# All distance matrices pre-calculated"; out.endl();
}
else {
    out << "# Distance matrices missing: Filling the gaps"; out.endl();
    temp_file_path = selfpath + MURMHED_db_file;
    error(load_file(&MURMHED_db,temp_file_path,DATA_MATRIX),matrix_created_files);

    out << "# Building filtered database distance matrix calculation"; out.endl();
    for(i=0; i < 200; i++) { //MURMHED_db.size()
        if(FALSE) {
            //PUT SHOWER SKIP HERE
        }
        else {
            OBS_temp.a = MURMHED_db[i][26];
            OBS_temp.e = MURMHED_db[i][27];
            OBS_temp.i = MURMHED_db[i][30];
            OBS_temp.omega = MURMHED_db[i][31];
            OBS_temp.Omega = MURMHED_db[i][29];

            OBS_temp.ra = MURMHED_db[i][9];
            OBS_temp.dec = MURMHED_db[i][10];
            OBS_temp.v_g = MURMHED_db[i][17];
            OBS_temp.lambda = MURMHED_db[i][36];

            real_objs.push_back(OBS_temp);
        }
    }
    out << "# " << real_objs.size() << " Objects built, " << (MURMHED_db.size() - real_objs.size()) << " shower members skipped"; out.endl();
    MURMHED_db.clear();

    out << "# Avalible RAM: " << opt.mem_allowed << " Mb"; out.endl();
    faction_of_avalible_mem = 0.001;
    mem_batch_estimate = ((double)real_objs.size())*((double)real_objs.size() - (double)1)*0.5*((double)sizeof(double));
    mem_batch_estimate = mem_batch_estimate/((double)(8*1024*1024)); //Mb

    calculation_batches = (unsigned int)ceil(mem_batch_estimate/(opt.mem_allowed*faction_of_avalible_mem));
    items_a_batch = (unsigned int)floor(((double)real_objs.size())/((double)calculation_batches));
    residue_items = real_objs.size() - calculation_batches*items_a_batch;

    if(residue_items > 0) {
        residue_toggle_loop = 1;
    }
    else {
        residue_toggle_loop = 0;
    }

    out << "# Estimating memory usage to " << mem_batch_estimate << " Mb, allowed usage is " << opt.mem_allowed*faction_of_avalible_mem << "Mb"; out.endl();
    out << "# Splitting work into " << calculation_batches << " batches, eatch of " << items_a_batch << " items, leaving " << residue_items << " residue items"; out.endl();
    for(j=0; j < FUNCTIONS.size(); j++) {
        if(!(folder_exists(MURMHED_D_mat_files[j]))) {

            boost::filesystem::create_directory(MURMHED_D_mat_files[j]);

            out << "# Calculating " << FUNCTIONS[j] << " distance matrix"; out.endl();
            for(k=0; k < calculation_batches+residue_toggle_loop; k++) {

                temp_batch_name = MURMHED_D_mat_files[j] + "/" + IntToString((unsigned int)(k*items_a_batch)) + "_" + IntToString((unsigned int)((k+1)*items_a_batch)) + "_" + FUNCTIONS[j] + ".txt";

                for(k1=0; k1 < (calculation_batches+residue_toggle_loop-k); k1++) {
                    real_objs_batch.clear();
                    real_objs_batch2.clear();
                    MURMHED_db_distance_matrix_part.clear();
                    std::cout << "- Doing batch coordinate (" << k << "," << k1 << "):" << std::endl;
                    if(k != calculation_batches+1) {
                        real_objs_batch.insert( real_objs_batch.begin(), real_objs.begin() + k*items_a_batch, real_objs.begin() + (k+1)*items_a_batch);
                    }
                    else {
                        real_objs_batch.insert( real_objs_batch.begin(), real_objs.begin() + k*items_a_batch, real_objs.begin() + k*items_a_batch + residue_items);
                    }
                    
                    if(k1 == 0) {
                        MURMHED_db_distance_matrix_part = calculate_metric_matrix(&real_objs_batch,FUNCTIONS[j]);
                        error(save_mat(temp_batch_name, &MURMHED_db_distance_matrix_part),matrix_created_files);
                    }
                    else if(k1 == (calculation_batches+1-k) && k != calculation_batches+1) {
                        real_objs_batch2.insert( real_objs_batch2.begin(), real_objs.begin() + k1*items_a_batch, real_objs.begin() + k1*items_a_batch + residue_items);
                        MURMHED_db_distance_matrix_part = augment_metric_matrix(&real_objs_batch,&real_objs_batch2,FUNCTIONS[j]);
                        error(save_batch(temp_batch_name, TEMP_file, &MURMHED_db_distance_matrix_part,k*items_a_batch),matrix_created_files);
                    }
                    else {
                        real_objs_batch2.insert( real_objs_batch2.begin(), real_objs.begin() + k1*items_a_batch, real_objs.begin() + (k1+1)*items_a_batch);
                        MURMHED_db_distance_matrix_part = augment_metric_matrix(&real_objs_batch,&real_objs_batch2,FUNCTIONS[j]);
                        error(save_batch(temp_batch_name, TEMP_file, &MURMHED_db_distance_matrix_part,k*items_a_batch),matrix_created_files);
                    }
                }
            }
        out << "# " << FUNCTIONS[j] << " distance matrix calculated"; out.endl();
        out << "# Saved " << FUNCTIONS[j] << " to file: " << MURMHED_D_mat_files[j]; out.endl();
        }
    }
    
    real_objs.clear();

    out << "# Database distance matrix calculation done, clearing data"; out.endl(); out.endl();   
}

DO_THEY_EXIST = TRUE;
for(i=0; i < MURMHED_profile_files.size(); i++) {
    DO_THEY_EXIST = DO_THEY_EXIST && file_exists(MURMHED_profile_files[i]);
}

if(DO_THEY_EXIST) {
    out << "# All database association profiles pre-calculated"; out.endl();
}
else {
    out << "# Missing database association profiles: Filling the gaps"; out.endl();

    for(j=0; j < FUNCTIONS.size(); j++) {
        if(!(file_exists(MURMHED_profile_files[j]))) {
            out << "# -- Creating database association profile " << MURMHED_profile_files[j]; out.endl();
            list_of_batches_d_matrix.clear();
            list_of_batches_d_matrix = list_files(MURMHED_D_mat_files[j]);

            index_matrix.clear();
            for(k=0; k < list_of_batches_d_matrix.size(); k++) {
                list_of_batches_d_matrix[k].erase(0,MURMHED_D_mat_files[j].size()+1);
                out << "Index range recorded for file " << list_of_batches_d_matrix[k]; out.endl();

                char_index_1 = list_of_batches_d_matrix[k].find("_", 0);
                char_index_2 = list_of_batches_d_matrix[k].find("_", char_index_1+1);

                index_matrix.push_back(EMPTY_UINT_VEC);
                buffer_temp = list_of_batches_d_matrix[k].substr(0,char_index_1);
                index_matrix[k].push_back((unsigned int)strtod(buffer_temp.c_str(), NULL));

                buffer_temp = list_of_batches_d_matrix[k].substr(char_index_1+1,char_index_2);
                index_matrix[k].push_back((unsigned int)strtod(buffer_temp.c_str(), NULL));
            }
            out << "Completed index matrix:"; out.endl();
            print_matrix(index_matrix,out);
            number_of_objects_db_load = index_matrix.back().back()-1;

            out << "Running cluster analysis"; out.endl();


                min_d_crit_step = 0.00001;
    max_d_crit_step = 1;
    association_adaptive_step_threshold = 0.2;//0.01;
    step_multiplyer_inc = 1.7;
    d_crit = 0;
    d_crit_step = 0.001;
    fraction_associated = 0;
    fraction_associated_next = 0;
    
    association_vec.clear();
    out << "# Starting association loop."; out.endl();
    association_loop_counter = 0;
    fraction_associated_old = 0;
    d_crit_vec.clear();
    while(fraction_associated < 1) {
        association_loop_counter++;
        
        cluster_matrix.clear();
        fractions.clear();

        cluster_matrix = single_linkage_clustering_files(MURMHED_D_mat_files[j], list_of_batches_d_matrix, index_matrix, d_crit);
        for(i=0; i < cluster_matrix.size(); i++) {
            if((double)cluster_matrix[i].size() > 1) {
                fractions.push_back((double)cluster_matrix[i].size());
            }
        }
        if(fractions.size() > 0) {
            fractions = fractions/((double)number_of_objects_db_load);
            fraction_associated = sum_v(fractions);
        }
        else {
            fraction_associated = 0;
        }
        
        if(fraction_associated - fraction_associated_old != 0) {
            d_crit_vec.push_back(d_crit);
            association_vec.push_back(fraction_associated);
        }

        fraction_associated_old = fraction_associated;
        
        d_crit_old = d_crit;
        d_crit += d_crit_step;

        visited_zero_progress_counter = 0;
        visited_max_progress_counter = 0;
        do {
            cluster_matrix_next.clear();
            fractions_next.clear();

            cluster_matrix_next = single_linkage_clustering_files(MURMHED_D_mat_files[j], list_of_batches_d_matrix, index_matrix, d_crit);
            for(i=0; i < cluster_matrix_next.size(); i++) {
                if((double)cluster_matrix_next[i].size() > 1) {
                    fractions_next.push_back((double)cluster_matrix_next[i].size());
                }
            }
            fractions_next = fractions_next/((double)number_of_objects_db_load);

            fraction_associated_next = sum_v(fractions_next);
            out << "# Checking next step: " << "frac = " << fraction_associated << ", frac_next = " << fraction_associated_next << ", with D_c = " << d_crit; out.endl();
            if((fraction_associated_next - fraction_associated) > association_adaptive_step_threshold && opt.adaptive_step) {
                visited_max_progress_counter++;
                next_step_ok = FALSE;
                d_crit = d_crit_old;
                if((d_crit_step*0.5) > min_d_crit_step) {
                    d_crit_step = d_crit_step*0.5;
                }
                else {
                    out << "# Minimum step reached, continuing: " << "frac = " << fraction_associated << ", frac_next = " << fraction_associated_next << ", with D_c = " << d_crit; out.endl();
                    next_step_ok = TRUE;
                }
                d_crit += d_crit_step;
            }
            else if((fraction_associated_next - fraction_associated) == 0 && opt.adaptive_step) {
                visited_zero_progress_counter++;
                next_step_ok = FALSE;
                d_crit = d_crit_old;
                if((d_crit_step*step_multiplyer_inc) < max_d_crit_step) {
                    d_crit_step = d_crit_step*step_multiplyer_inc;
                }
                else {
                    out << "# Maximum step reached, continuing: " << "frac = " << fraction_associated << ", frac_next = " << fraction_associated_next << ", with D_c = " << d_crit; out.endl();
                    next_step_ok = TRUE;
                }
                d_crit += d_crit_step;
            }
            else {
                next_step_ok = TRUE;
            }

            if(visited_zero_progress_counter > visited_max_progress_counter && visited_max_progress_counter > 10) {
                out << "# Cluster too small, continuing: " << "frac = " << fraction_associated << ", frac_next = " << fraction_associated_next << ", with D_c = " << d_crit; out.endl();
                next_step_ok = TRUE;
            }
        } while(!next_step_ok && fraction_associated != 1);

        out << "Sample associated " << fraction_associated*100 << " percent with " << d_crit_old << " threshold."; out.endl();
    }
    
    out << "# Association loop exited."; out.endl();
    
    out << "# Calculating new means."; out.endl();
    
    //j*2       d_crit vectors
    //j*2+1     frac vectors

    association_profile_db[j*2].insert( association_profile_db[j*2].end(), d_crit_vec.begin(), d_crit_vec.end() );
    association_profile_db[j*2+1].insert( association_profile_db[j*2+1].end(), association_vec.begin(), association_vec.end() );

    order_D_c_vec.clear();
    order_D_c_vec.resize(association_profile_db[j*2].size());

    D_c_n = 0;
    for (double_vector_iterator it = association_profile_db[j*2].begin(); it != association_profile_db[j*2].end(); ++it, ++D_c_n) {
        order_D_c_vec[D_c_n] = std::make_pair(D_c_n, it);
    }

    std::sort(order_D_c_vec.begin(), order_D_c_vec.end(), double_vec_ordering_struct());

    association_profile_db[j*2] = sort_from_ref((association_profile_db[j*2]),order_D_c_vec);
    association_profile_db[j*2+1] = sort_from_ref((association_profile_db[j*2+1]),order_D_c_vec);
    bin_member_counter = 1;
    for(i=association_profile_db[j*2].size()-1; i > 0; i--) {
        if(association_profile_db[j*2][i] == association_profile_db[j*2][i-1]) {
            bin_member_counter++;
            association_profile_db[j*2+1][i-1] = association_profile_db[j*2+1][i-1] + association_profile_db[j*2+1][i];

            association_profile_db[j*2].erase(association_profile_db[j*2].begin()+i);
            association_profile_db[j*2+1].erase(association_profile_db[j*2+1].begin()+i);
        }
        else {
            std::cout << "Merging duplicate data points to mean: bin total=" << association_profile_db[j*2+1][i-1] << ", n=" << bin_member_counter <<  ", result=" << association_profile_db[j*2+1][i-1]/((double)bin_member_counter) <<  ", d_c=" << association_profile_db[j*2][i-1] << std::endl;
            if(bin_member_counter != 1) {
                association_profile_db[j*2+1][i-1] = association_profile_db[j*2+1][i-1]/((double)bin_member_counter);
                bin_member_counter = 1;
            }
        }
    }

    out << "########################################################################"; out.endl();

       }     
    }

    out << "# Saving data."; out.endl();
for(i=0; i < FUNCTIONS.size(); i++) {
    if(!(file_exists(MURMHED_profile_files[i]))) {
    association_profile_db_save.clear();
    association_profile_db_save.push_back(association_profile_db[i*2]);
    association_profile_db_save.push_back(association_profile_db[i*2+1]);

    association_profile_db_save_means.clear();
    association_profile_db_save_means.push_back(EMPTY_DOUBLE_VEC);
    association_profile_db_save_means.push_back(EMPTY_DOUBLE_VEC);

    if(opt.d_crit_bins == 0) {
        d_crit_bin_n = round(sqrt((double)(association_profile_db_save[0].size())));
    }
    else {
        d_crit_bin_n = opt.d_crit_bins;
    }
    d_crit_bin_width = (association_profile_db_save[0].back() - association_profile_db_save[0][0])/((double)d_crit_bin_n);


    for(j=0; j < d_crit_bin_n; j++) {
        association_profile_db_save_means[0].push_back(association_profile_db_save[0][0]+(j+0.5)*d_crit_bin_width);
        association_profile_db_save_means[1].push_back(-1);
    }

    /*std::cout << "PRINTING MATRIX DEBUG MODE: " << std::endl;
    print_matrix(association_profile_db_save_means);
    std::cout << std::endl << std::endl;*/

    out << "## Starting bin avrage calculation"; out.endl();
    out << "Bin number: " << d_crit_bin_n << ", max val: " << association_profile_db_save[0].back(); out.endl();
    k = 0;
    for(j=0; j < association_profile_db_save[1].size(); j++) {
        if(association_profile_db_save[0][j] > (association_profile_db_save[0][0] + (k+1)*d_crit_bin_width)) {
            out << "-- Ending bin at " << association_profile_db_save_means[1][k] << ", to " << association_profile_db_save_means[1][k]/((double)bin_member_counter) << ", using n=" << bin_member_counter; out.endl();
            association_profile_db_save_means[1][k] = association_profile_db_save_means[1][k]/((double)bin_member_counter);
            k++;

            while(association_profile_db_save[0][j] > (association_profile_db_save[0][0] + (k+1)*d_crit_bin_width)) {
                out << "-- Skipping bin at " << (association_profile_db_save[0][0] + (k+1)*d_crit_bin_width); out.endl();
                association_profile_db_save_means[1][k] = 0;
                k++;
            }
        }

        //std::cout << "current bin val: " << association_profile_db_save_means[1][k] << ", j: " << j <<  ", k: " << k << std::endl;
        if(association_profile_db_save_means[1][k] == -1) {
            bin_member_counter = 1;
            out << "-- Starting new bin at: (" << association_profile_db_save[0][j] << "," << association_profile_db_save[1][j] << "), bin center and width: (" << association_profile_db_save_means[0][k] << "," << d_crit_bin_width << ")"; out.endl();
            association_profile_db_save_means[1][k] = association_profile_db_save[1][j];
        }
        else {
            out << "-- Adding value to bin: (" << association_profile_db_save[0][j] << "," << association_profile_db_save[1][j] << "), bin value: " << association_profile_db_save_means[1][k]; out.endl();
            bin_member_counter++;
            association_profile_db_save_means[1][k] = association_profile_db_save_means[1][k] + association_profile_db_save[1][j];
        }

    }
    out << "-- Ending bin at " << association_profile_db_save_means[1][k] << ", to " << association_profile_db_save_means[1][k]/((double)bin_member_counter) << ", using n=" << bin_member_counter; out.endl();
    association_profile_db_save_means[1][k] = association_profile_db_save_means[1][k]/((double)bin_member_counter);

    for(j=0; j < (d_crit_bin_n-1); j++) {
        if(association_profile_db_save_means[1][j] == 0 && j > 0) {
            empty_bin_counter=j;
            while(association_profile_db_save_means[1][empty_bin_counter] == 0) {
                empty_bin_counter++;
            }
            out << "-- Discovered empty bins, extrapolating from " << association_profile_db_save_means[1][j-1] << " to " << association_profile_db_save_means[1][empty_bin_counter] << ", missing bins=" << empty_bin_counter-j;out.endl();
            for(k=j; k < empty_bin_counter; k++) {
                association_profile_db_save_means[1][k] = association_profile_db_save_means[1][j-1] + ( ((double)(k-j+1))/((double)(empty_bin_counter-j+1)) )*(association_profile_db_save_means[1][empty_bin_counter] - association_profile_db_save_means[1][j-1]);
                out << "-- Assigning (" << association_profile_db_save_means[0][k] << "," << association_profile_db_save_means[1][k] << "), ";
            }
            out.endl();
        }
    }

    out << "- Transposing data matrix."; out.endl();
    association_profile_db_save_means = matrix_transpose(association_profile_db_save_means);

    save_mat(MURMHED_profile_files[i], &association_profile_db_save_means);
    out << "- Saved file: " << MURMHED_profile_files[i]; out.endl();
    }
}
}

unsigned int last_total_PB_created;

    // FROM HERE NO MORE DECLARATIONS -----------------------------------
out << "##############################"; out.endl();
out << "###### Entering MC mode ######"; out.endl();
out << "##############################"; out.endl();

total_showers_created = 0;
total_n_created = 0;
total_PB_created = 0;
total_met_collided = 0;

bool MC_not_done = TRUE;
bool not_done = TRUE;
do {
out << ". Progress_check_mc_loop_start"; out.endl();
not_done = TRUE;
do {
out << ". Progress_check_PB_gen_start"; out.endl();
    //GENERATE A IC FOR PBE

    for(i=0; i < 4; i++) {
        switch((unsigned int)body_family[i][0]) {
            case 0: 
                current_body_char[i] = draw_from_uniform(body_family[i][1],body_family[i][2]);
                break;
        }
    }
    body_mass = 0.75*PI*pow(current_body_char[0],3)*current_body_char[1];

    out << ". Progress_check_PB_ok_orbit_check_start"; out.endl();

    last_total_PB_created = total_PB_created;
    do {
    total_PB_created++;


    current_body_kep[0] = semi_axis.hist_bins[draw_from_dist(semi_axis.PDF)];
    current_body_kep[1] = ecc.hist_bins[draw_from_dist(ecc.PDF)];
    current_body_kep[2] = inc.hist_bins[draw_from_dist(inc.PDF)];
    current_body_kep[3] = omega.hist_bins[draw_from_dist(omega.PDF)];
    current_body_kep[4] = Omega.hist_bins[draw_from_dist(Omega.PDF)];

      if(current_body_kep[0]*(1-current_body_kep[1])*1.2 < current_body_char[3]) {
        orbit_ok = TRUE;
      }
      else {
        orbit_ok = FALSE;

        switch((unsigned int)body_family[3][0]) {
            case 0: 
                current_body_char[3] = draw_from_uniform(body_family[3][1],body_family[3][2]);
                break;
        }
      }

    std::cout << "\r-- Checking PB " << total_PB_created - last_total_PB_created << ", a=" << current_body_kep[0] << ", e=" << current_body_kep[1] << ", i=" << current_body_kep[2] << ", om=" << current_body_kep[3] << ", Om=" << current_body_kep[4] << ", r_c=" << current_body_char[3] << ", 1.2*q=" << current_body_kep[0]*(1-current_body_kep[1])*1.2 << std::flush;
    } while(!orbit_ok);
    std::cout << "\r-- Checked PB " << total_PB_created - last_total_PB_created << "                      " << std::endl;
    out << ". Progress_check_PB_selected after " << total_PB_created - last_total_PB_created << " loops"; out.endl();
    current_body_char[3] = current_body_char[3]*AU;
    current_body_kep[0] = current_body_kep[0]*AU; //convert from AU to m (for calc consistency)
    current_body_kep[2] = current_body_kep[2]*(PI/180); //radians
    current_body_kep[3] = current_body_kep[3]*(PI/180);
    current_body_kep[4] = current_body_kep[4]*(PI/180);
    current_body_kep[5] = PI; //true anoamly, start at 0=peri, pi=appihelion
    
    temp_vec = kepler_to_xv(current_body_kep, G*(body_mass+Msol)); //sun mass

    out << "-------- IC for PBE configuring to: "; out.endl();
    out << "a = " << current_body_kep[0]/AU << ", e = " << current_body_kep[1] << ", i = " << current_body_kep[2]; out.endl();
    out << "omega = " << current_body_kep[3] << ", Omega = " << current_body_kep[4] << ", nu = " << current_body_kep[5]; out.endl();
    out << "R_c = " << current_body_char[0] << ", m_c = " << body_mass << ", sigma = " << current_body_char[1]; out.endl();
    out << "alpha = " << current_body_char[2] << ", r_c = " << current_body_char[3]/AU; out.endl(); out.endl();
    out << "-------- IC for PBE (q,v): " << sqrt(pow(temp_vec[0],2) + pow(temp_vec[1],2) + pow(temp_vec[2],2))/AU << " " << sqrt(pow(temp_vec[3],2) + pow(temp_vec[4],2) + pow(temp_vec[5],2))*0.001 ; out.endl();
    out << "q: " << temp_vec[0] << " " << temp_vec[1] << " " << temp_vec[2] ; out.endl();
    out << "v: " << temp_vec[3]/1000 << " " << temp_vec[4]/1000 << " " << temp_vec[5]/1000; out.endl() ; out.endl();
    

    out_file_name = PBE_IC_body;
    out_file_name.insert(0,PBE_IC_PRODUCER_folder);

    if(file_exists(out_file_name)) {
      remove(out_file_name.c_str());
      out << "Removed file: " << out_file_name; out.endl();
    }

    out_stream.open(out_file_name.c_str(), std::ios::app);
    //out_stream << std::setprecision((*opt).precision);
    if(out_stream.fail()) {
        error(DATA_SAVE_FAILED,matrix_created_files);
    }

    out_stream << temp_vec[0] << " " << temp_vec[1] << " " << temp_vec[2];
    out_stream << "\r\n";
    out_stream << temp_vec[3] << " " << temp_vec[4] << " " << temp_vec[5];
    out_stream << "\r\n";
    out_stream << current_body_char[0] << " " << current_body_char[1] << " " << current_body_char[2] << " " << current_body_char[3];
    out_stream << "\r\n";

    out_stream.close();
    out_stream.clear();
    out << "File " << out_file_name << " created"; out.endl();


    execute_command = PBE_program;
    execute_command.insert(0,selfpath);
    execute_command = execute_command + " " + selfpath + settings_folder + PBE_config;
    out << "Calling " << execute_command; out.endl();
    PBE_status = system(execute_command.c_str());
    
    out << "PBE exit status: " << PBE_status; out.endl();
    error(PBE_status,matrix_created_files);

    out << "Using integrator: " << opt.integrator; out.endl();
    switch(opt.integrator) {
        case 0:
            rename(PBE_q_out.c_str(), CMS_q_in.c_str());
            rename(PBE_v_out.c_str(), CMS_v_in.c_str());
            rename(PBE_m_out.c_str(), CMS_m_in.c_str());

            execute_command = CMS_program;
            execute_command.insert(0,selfpath);
            execute_command = execute_command + " " + selfpath + settings_folder + CMS_config + " 1";
            out << "Calling " << execute_command; out.endl();
            CMS_status = system(execute_command.c_str());

            out << "CMS exit status: " << CMS_status; out.endl();
            error(CMS_status,matrix_created_files);

            break;
        case 1:
            out_stream.open(mercury6_small.c_str(), std::ios::app);
            out_stream << std::uppercase << std::setprecision(17) << std::scientific;

            temp_q.clear(); temp_v.clear(); temp_m.clear(); temp_t.clear(); temp_sim.clear();
            error(load_data_matrix(&temp_q,PBE_q_out),matrix_created_files);
            error(load_data_matrix(&temp_v,PBE_v_out),matrix_created_files);
            error(load_data_matrix(&temp_t,PBE_t_out),matrix_created_files);
            error(load_data_matrix(&temp_sim,PBE_sim_out),matrix_created_files);
            
            //error(load_data_matrix(&temp_m,PBE_m_out));

            out_stream << ")O+_06 Small-body initial data  (WARNING: Do not delete this line!!)" << line_end;
            out_stream << ") Lines beginning with `)' are ignored." << line_end;
            out_stream << ")---------------------------------------------------------------------" << line_end;
            out_stream << " style (Cartesian, Asteroidal, Cometary) = Car" << line_end;
            out_stream << ")---------------------------------------------------------------------" << line_end;
            total_n_created += temp_q.size();
            current_n_created = temp_q.size();
            if(temp_q.size() > 4000) {
                dummy_var = 4000;
            }
            else {
                dummy_var = temp_q.size();
            }
            for(i=0; i < dummy_var; i++) { //temp_q.size()
                TRUTH_CHECKER = TRUE;
                TRUTH_CHECKER = TRUTH_CHECKER && is_finite(temp_q[i][0]);
                TRUTH_CHECKER = TRUTH_CHECKER && is_finite(temp_q[i][1]);
                TRUTH_CHECKER = TRUTH_CHECKER && is_finite(temp_q[i][2]);
                TRUTH_CHECKER = TRUTH_CHECKER && is_finite(temp_v[i][0]);
                TRUTH_CHECKER = TRUTH_CHECKER && is_finite(temp_v[i][1]);
                TRUTH_CHECKER = TRUTH_CHECKER && is_finite(temp_v[i][2]);
                if(TRUTH_CHECKER) {
                    out_stream << " " << "OBJ" << i << "     ep=" << 2451179.5 - (temp_sim[0][1] - temp_t[0][i])/(3600*24*365) << line_end;
                    out_stream << " " << temp_q[i][0]/AU << " " << temp_q[i][1]/AU << " " << temp_q[i][2]/AU;
                    out_stream << " " << temp_v[i][0]/AU*86400 << " " << temp_v[i][1]/AU*86400 << " " << temp_v[i][2]/AU*86400;
                    out_stream << " " << "0. 0. 0." << line_end;
                }
                else {
                    out << "SOMETHING HORRIBLE HAPPEND TO THE TEST PARTICLE DATA! Skipped object in mercury small.in file."; out.endl();
                }
            }

            out_stream.close();
            out_stream.clear();
            temp_q.clear(); temp_v.clear(); temp_m.clear();

            out << "########################################################################"; out.endl();
            out << "# Total_met_collided = " << total_met_collided; out.endl();
            out << "# Total_showers_created = " << total_showers_created << ", " << opt.showers - total_showers_created << " left"; out.endl();
            out << "# Starting new particle propagation"; out.endl();
            out << "########################################################################"; out.endl();

            execute_command = "cd " + selfpath + mercury6_folder + "; " + mercury6_program;
            out << "Calling " << execute_command; out.endl();
            mercury6_status = system(execute_command.c_str());

            out << "mercury6 exit status: " << mercury6_status; out.endl();
            error(mercury6_status,matrix_created_files);
            //exit(-1);
            execute_command = "cd " + selfpath + mercury6_folder + "; " + close6_program;
            out << "Calling " << execute_command; out.endl();
            mercury6_status = system(execute_command.c_str());

            out << "close6 exit status: " << mercury6_status; out.endl();
            error(mercury6_status,matrix_created_files);
            //exit(-1);
            error(read_mercury6_matrix(&encounter_data,mercury6_earth_clo),matrix_created_files);

            if(encounter_data.size() > 0) {
                id_vec.clear();
                orig_id_vec.clear();
                filtered_encounter_data.clear();
                out << "Encounters detected, checking hill radii of: " << encounter_data.size() << " objects"; out.endl();
                //Check for encounters inside hill radii
                for(i=0; i < encounter_data.size(); i=i+2) {
                    if(encounter_data[i][2] <= earth_hill_r) {
                        id_vec.push_back((unsigned int)encounter_data[i][1]);
                        orig_id_vec.push_back((unsigned int)encounter_data[i][1]);
                        filtered_encounter_data.push_back(encounter_data[i]);
                        filtered_encounter_data.back().insert(filtered_encounter_data.back().end(),encounter_data[i+1].begin(),encounter_data[i+1].end());
                    }
                }

                encounter_data.clear();
                out <<"Encounters filtered down to: " << filtered_encounter_data.size() << " objects"; out.endl();
                if(filtered_encounter_data.size() > 1) {

                    temp_mat.clear(); temp_vec.clear();
                    error(load_data_matrix(&temp_mat,PBE_abs_v_out),matrix_created_files);
                    save_mat(PBE_abs_v_statistics, &temp_mat);

                    temp_mat.clear(); temp_vec.clear();
                    error(load_data_matrix(&temp_mat,PBE_sim_out),matrix_created_files);
                    temp_vec.push_back(temp_mat[0][0]);
                    temp_mat.clear();
                    temp_mat.push_back(temp_vec);
                    save_mat(PBE_weight_statistics, &temp_mat);

                    temp_mat.clear(); temp_vec.clear();
                    error(load_data_matrix(&temp_mat,PBE_M_out),matrix_created_files);
                    temp_mat = matrix_transpose(temp_mat);
                    save_mat(PBE_M_statistics, &temp_mat);

                    temp_mat.clear(); temp_vec.clear();
                    total_showers_created++;
                    not_done = FALSE;
                }
                else {
                    clean_up(matrix_created_files,1);
                    clean_up(matrix_created_files,2);
                    clean_up(matrix_created_files,3);
                }
            }
            else {
                clean_up(matrix_created_files,1);
                clean_up(matrix_created_files,2);
                clean_up(matrix_created_files,3);
            }

total_met_collided+= filtered_encounter_data.size();

/*std::cout << "filtered_encounter_data: " << std::endl;
print_matrix(filtered_encounter_data);

exit(-1);*/

encounter_stat_data.clear();
temp_vec.clear();
temp_vec.push_back(total_PB_created); // PB ID
temp_vec.push_back(total_n_created);
temp_vec.push_back(current_n_created);
temp_vec.push_back(total_met_collided);
temp_vec.push_back(filtered_encounter_data.size());

encounter_stat_data.push_back(temp_vec);
save_mat(statistics_PB_n_out, &encounter_stat_data);

            break;
    }

} while(not_done);



//Filter our unique objects
std::sort(id_vec.begin(),id_vec.end());
id_vec.erase(std::unique(id_vec.begin(),id_vec.end()),id_vec.end());

out << "DEBUG PRINTING id_vec: "; out.endl();
print_vector(id_vec,out);

//Remove config file for rewriting
if(file_exists(mercury6_element_cfg)) {
    remove(mercury6_element_cfg.c_str());
    out << "Removed file: " << mercury6_element_cfg; out.endl();
}
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
out_stream << "EARTH" << line_end;

for(i=0; i < id_vec.size(); i++) {
    out_stream << "OBJ" << id_vec[i] << line_end;
}
out_stream.close();
out_stream.clear();
out << "Created file: " << mercury6_element_cfg; out.endl();

execute_command = "cd " + selfpath + mercury6_folder + "; " + element6_program;
out << "Calling " << execute_command; out.endl();
mercury6_status = system(execute_command.c_str());
out << "element6 exit status: " << mercury6_status; out.endl();
error(mercury6_status,matrix_created_files);

encounter_data.clear();
earth_encounter_data.clear();
filtered_encounter_data_save.clear();

for(i=0; i < id_vec.size(); i++) {
    temp_obj_data.clear();
    mercury6_temp_obj_str = selfpath + mercury6_folder + "OBJ" + IntToString(id_vec[i]) + ".aei";
    out << "Attempting to load: " << mercury6_temp_obj_str; out.endl();
    error(read_mercury6_matrix(&temp_obj_data,mercury6_temp_obj_str),matrix_created_files);

    for(j=0; j < filtered_encounter_data.size(); j++ ) {
        if(orig_id_vec[j] == id_vec[i]) {
            first_id = j;
            break;
        }
    }

    temp_data = abs_element_v(extract_col(temp_obj_data,0)-(filtered_encounter_data[first_id][0]));
    TEMP_ID = std::min_element(temp_data.begin(),temp_data.end());
    TEMP_ID_index = (unsigned int)(TEMP_ID - temp_data.begin());

    if((filtered_encounter_data[first_id][0] - temp_obj_data[TEMP_ID_index][0]) < 0 && TEMP_ID_index >= 1) {
        TEMP_ID_index = TEMP_ID_index - 1;
    }

    temp_vec.clear();
    temp_vec.push_back(filtered_encounter_data[first_id][0]);
    temp_vec.push_back(filtered_encounter_data[first_id][2]);
    temp_vec.push_back(filtered_encounter_data[first_id][9]);
    temp_vec.push_back(filtered_encounter_data[first_id][10]);
    temp_vec.push_back(filtered_encounter_data[first_id][11]);
    temp_vec.push_back(filtered_encounter_data[first_id][12]);
    temp_vec.push_back(filtered_encounter_data[first_id][13]);
    temp_vec.push_back(filtered_encounter_data[first_id][14]);

    filtered_encounter_data_save.push_back(temp_vec);

    encounter_data.push_back(temp_obj_data[TEMP_ID_index]);

    temp_vec.clear();
    temp_vec.push_back(filtered_encounter_data[first_id][0]);
    temp_vec.push_back(filtered_encounter_data[first_id][2]);
    temp_vec.push_back(filtered_encounter_data[first_id][3]);
    temp_vec.push_back(filtered_encounter_data[first_id][4]);
    temp_vec.push_back(filtered_encounter_data[first_id][5]);
    temp_vec.push_back(filtered_encounter_data[first_id][6]);
    temp_vec.push_back(filtered_encounter_data[first_id][7]);
    temp_vec.push_back(filtered_encounter_data[first_id][8]);

    earth_encounter_data.push_back(temp_vec);

}

for(i=0; i < id_vec.size(); i++) {
    mercury6_temp_obj_str = selfpath + mercury6_folder + "OBJ" + IntToString(id_vec[i]) + ".aei";
    if(file_exists(mercury6_temp_obj_str)) {
        remove(mercury6_temp_obj_str.c_str());
        out << "Removed file: " << mercury6_temp_obj_str; out.endl();
    }
}

encounter_stat_data.clear();
encounter_stat_data = filtered_encounter_data_save;
for(i=0; i < encounter_stat_data.size(); i++) {
    encounter_stat_data[i].insert(encounter_stat_data[i].begin(),total_PB_created);
}
save_mat(met_encounter_data_out, &encounter_stat_data);

encounter_stat_data.clear();
encounter_stat_data = encounter_data;
for(i=0; i < encounter_stat_data.size(); i++) {
    encounter_stat_data[i].insert(encounter_stat_data[i].begin(),total_PB_created);
}
save_mat(met_before_encounter_data_out, &encounter_stat_data);

encounter_stat_data.clear();
encounter_stat_data = earth_encounter_data;
for(i=0; i < encounter_stat_data.size(); i++) {
    encounter_stat_data[i].insert(encounter_stat_data[i].begin(),total_PB_created);
}
save_mat(earth_encounter_data_out, &encounter_stat_data);

encounter_stat_data.clear();
encounter_stat_data.push_back(current_body_char);
save_mat(statistics_PB_char_out, &encounter_stat_data);

encounter_stat_data.clear();
encounter_stat_data.push_back(current_body_kep);
    encounter_stat_data[0][0] = encounter_stat_data[0][0]/AU; //convert from m to AU (for output consistency)
    encounter_stat_data[0][2] = encounter_stat_data[0][2]/(PI/180); //deg
    encounter_stat_data[0][3] = encounter_stat_data[0][3]/(PI/180);
    encounter_stat_data[0][4] = encounter_stat_data[0][4]/(PI/180);
save_mat(statistics_PB_kep_out, &encounter_stat_data);

 encounter_stat_data.clear();
//print_matrix(encounter_data);

//print_matrix(encounter_stat_data);


out << "########################################################################"; out.endl();
out << "# STARTING ASSOCIATION                                                 #"; out.endl();
out << "########################################################################"; out.endl();



//out; out.endl() << "##################### total_met_collided = " << total_met_collided << " ########################## "; out.endl(); out.endl();
//std::cin.get();

out << "# Objects created."; out.endl();
objs.clear();
for(i=0; i < encounter_data.size(); i++) {
    OBS_temp.a = encounter_data[i][1];
    OBS_temp.e = encounter_data[i][2];
    OBS_temp.i = encounter_data[i][3];
    OBS_temp.omega = encounter_data[i][4];
    OBS_temp.Omega = encounter_data[i][5];

    objs.push_back(OBS_temp);
}

min_d_crit_step = 0.00001;
max_d_crit_step = 1;
if(encounter_data.size() < 10) {
    association_adaptive_step_threshold = 2/encounter_data.size();
}
else if(encounter_data.size() < 100 && encounter_data.size() >= 10) {
    association_adaptive_step_threshold = 0.2;
}
else {
    association_adaptive_step_threshold = 0.05;
}
    out << "########################################################################"; out.endl();
    out << "# Total_met_collided = " << total_met_collided; out.endl();
    out << "########################################################################"; out.endl();
for(j=0; j < FUNCTIONS.size(); j++) {

    //D_matix = calculate_metric_matrix(&objs, d_SH);
    D_matix.clear(); error_function.clear();
    D_matix = calculate_metric_matrix(&objs,FUNCTIONS[j]);

    out << "DEBUG PRINTING MATRIX: "; out.endl();
    print_matrix(D_matix,out);

    out << "# D_matrix with D_" << FUNCTIONS[j] << " calculated."; out.endl();
    d_crit = 0;
    d_crit_step = 0.001;
    fraction_associated = 0;
    fraction_associated_next = 0;
    
    association_vec.clear();
    out << "# Starting association loop."; out.endl();
    association_loop_counter = 0;
    fraction_associated_old = 0;
    association_loop_error_counter = 0;
    d_crit_vec.clear();
    do {
        association_loop_counter++;
        /*if(association_loop_counter >= association_loop_counter_max) {
            out << "# LOOP MAY BE STUCK, 1000 ITTERATIONS PASSED."; out.endl();
            out << "# D_c = " << d_crit << ", fraction_associated = " << fraction_associated; out.endl();
        }*/
        /*if(association_loop_counter % 1000 == 0) {
            out << "# 1000 loops passed: " << "fraction_associated = " << fraction_associated << ", with D_c = " << d_crit; out.endl();
        }*/
        /*if(association_loop_counter > 10000) {
            out << "# 10000 loops passed: SOMETHING WENT WRONG " << "fraction_associated = " << fraction_associated << ", with D_c = " << d_crit; out.endl();
            out << "Distance matrix debug printout: "; out.endl();
            print_matrix(D_matix);
            break;
        }*/

        
        cluster_matrix.clear();
        fractions.clear();

        cluster_matrix = single_linkage_clustering_matrix(&D_matix, d_crit);

        error_function.push_back( (double)cluster_matrix.size()/((double)encounter_data.size()) );

        for(i=0; i < cluster_matrix.size(); i++) {
            if((double)cluster_matrix[i].size() > 1) {
                fractions.push_back((double)cluster_matrix[i].size());
            }
        }
        if(fractions.size() > 0) {
            fractions = fractions/((double)encounter_data.size());
            fraction_associated = sum_v(fractions);
        }
        else {
            fraction_associated = 0;
        }
        
        if(fraction_associated - fraction_associated_old != 0) {
            d_crit_vec.push_back(d_crit);
            association_vec.push_back(fraction_associated);
        }

        fraction_associated_old = fraction_associated;
        
        d_crit_old = d_crit;
        d_crit += d_crit_step;

        visited_zero_progress_counter = 0;
        visited_max_progress_counter = 0;
        do {
            cluster_matrix_next.clear();
            fractions_next.clear();

            cluster_matrix_next = single_linkage_clustering_matrix(&D_matix, d_crit);
            for(i=0; i < cluster_matrix_next.size(); i++) {
                if((double)cluster_matrix_next[i].size() > 1) {
                    fractions_next.push_back((double)cluster_matrix_next[i].size());
                }
            }
            fractions_next = fractions_next/((double)encounter_data.size());

            fraction_associated_next = sum_v(fractions_next);
            //out << "# Checking next step: " << "frac = " << fraction_associated << ", frac_next = " << fraction_associated_next << ", with D_c = " << d_crit; out.endl();
            if((fraction_associated_next - fraction_associated) > association_adaptive_step_threshold && opt.adaptive_step) {
                visited_max_progress_counter++;
                next_step_ok = FALSE;
                d_crit = d_crit_old;
                if((d_crit_step*0.5) > min_d_crit_step) {
                    d_crit_step = d_crit_step*0.5;
                }
                else {
                    out << "# Minimum step reached, continuing: " << "frac = " << fraction_associated << ", frac_next = " << fraction_associated_next << ", with D_c = " << d_crit; out.endl();
                    next_step_ok = TRUE;
                }
                d_crit += d_crit_step;
            }
            else if((fraction_associated_next - fraction_associated) == 0 && opt.adaptive_step) {
                visited_zero_progress_counter++;
                next_step_ok = FALSE;
                d_crit = d_crit_old;
                if((d_crit_step*1.2) < max_d_crit_step) {
                    d_crit_step = d_crit_step*1.2;
                }
                else {
                    out << "# Maximum step reached, continuing: " << "frac = " << fraction_associated << ", frac_next = " << fraction_associated_next << ", with D_c = " << d_crit; out.endl();
                    next_step_ok = TRUE;
                }
                d_crit += d_crit_step;
            }
            else {
                next_step_ok = TRUE;
            }

            if(visited_zero_progress_counter > visited_max_progress_counter && visited_max_progress_counter > 10) {
                out << "# Cluster too small, continuing: " << "frac = " << fraction_associated << ", frac_next = " << fraction_associated_next << ", with D_c = " << d_crit; out.endl();
                next_step_ok = TRUE;
            }

            if(fraction_associated_next-fraction_associated == 0) {
                association_loop_error_counter++;
            }
            else {
                association_loop_error_counter = 0;
            }
            if(association_loop_error_counter > 100) {
                out << "# AN ERROR MIGHT HAVE OCCURED, EXITING LOOP: " << "association_loop_error_counter = " << association_loop_error_counter; out.endl();
                out << "#                                          : " << "cluster_matrix.size() = " << cluster_matrix.size(); out.endl();
                out << "Printing D_matrix:"; out.endl();
                print_matrix(D_matix,out);
                out << "Printing encounter_data:"; out.endl();
                print_matrix(encounter_data,out);
                break;
            }
        } while(!next_step_ok && cluster_matrix.size() > 1);

        out << "Sample associated " << fraction_associated*100 << " percent with " << d_crit_old << " threshold."; out.endl();
        out << "Cluster fraction " << (1 - ((double)cluster_matrix.size())/((double)encounter_data.size()))*100 << " percent with " << d_crit_old << " threshold."; out.endl();
        if(association_loop_error_counter > 100) {
                break;
        }
    } while(cluster_matrix.size() > 1);
    
    out << "# Association loop exited."; out.endl();

    out << "# Saving min err function Dc to file " << statistics_Dc_min_err[j]; out.endl();
    temp_mat.clear(); temp_vec.clear();
    temp_vec.push_back(d_crit_vec.back());
    temp_mat.push_back(temp_vec);
    save_mat(statistics_Dc_min_err[j], &temp_mat);
    
    
    out << "# Calculating new means."; out.endl();
    
    //j*2       d_crit vectors
    //j*2+1     frac vectors

    association_profile[j*2].insert( association_profile[j*2].end(), d_crit_vec.begin(), d_crit_vec.end() );
    association_profile[j*2+1].insert( association_profile[j*2+1].end(), association_vec.begin(), association_vec.end() );

    order_D_c_vec.clear();
    order_D_c_vec.resize(association_profile[j*2].size());

    D_c_n = 0;
    for (double_vector_iterator it = association_profile[j*2].begin(); it != association_profile[j*2].end(); ++it, ++D_c_n) {
        order_D_c_vec[D_c_n] = std::make_pair(D_c_n, it);
    }

    std::sort(order_D_c_vec.begin(), order_D_c_vec.end(), double_vec_ordering_struct());

    association_profile[j*2] = sort_from_ref((association_profile[j*2]),order_D_c_vec);
    association_profile[j*2+1] = sort_from_ref((association_profile[j*2+1]),order_D_c_vec);
    bin_member_counter = 1;
    for(i=association_profile[j*2].size()-1; i > 0; i--) {
        if(association_profile[j*2][i] == association_profile[j*2][i-1]) {
            bin_member_counter++;
            association_profile[j*2+1][i-1] = association_profile[j*2+1][i-1] + association_profile[j*2+1][i];

            association_profile[j*2].erase(association_profile[j*2].begin()+i);
            association_profile[j*2+1].erase(association_profile[j*2+1].begin()+i);
        }
        else {
            std::cout << "Merging duplicate data points to mean: bin total=" << association_profile[j*2+1][i-1] << ", n=" << bin_member_counter <<  ", result=" << association_profile[j*2+1][i-1]/((double)bin_member_counter) <<  ", d_c=" << association_profile[j*2][i-1] << std::endl;
            if(bin_member_counter != 1) {
                association_profile[j*2+1][i-1] = association_profile[j*2+1][i-1]/((double)bin_member_counter);
                bin_member_counter = 1;
            }
        }
    }

    out << "########################################################################"; out.endl();
}

out << "# Saving data."; out.endl();
for(i=0; i < FUNCTIONS.size(); i++) {
    association_profile_save.clear();
    association_profile_save.push_back(association_profile[i*2]);
    association_profile_save.push_back(association_profile[i*2+1]);

    association_profile_save_means.clear();
    association_profile_save_means.push_back(EMPTY_DOUBLE_VEC);
    association_profile_save_means.push_back(EMPTY_DOUBLE_VEC);

    if(opt.d_crit_bins == 0) {
        d_crit_bin_n = round(sqrt((double)(association_profile_save[0].size())));
    }
    else {
        d_crit_bin_n = opt.d_crit_bins;
    }
    d_crit_bin_width = (association_profile_save[0].back() - association_profile_save[0][0])/((double)d_crit_bin_n);


    for(j=0; j < d_crit_bin_n; j++) {
        association_profile_save_means[0].push_back(association_profile_save[0][0]+(j+0.5)*d_crit_bin_width);
        association_profile_save_means[1].push_back(-1);
    }

    /*std::cout << "PRINTING MATRIX DEBUG MODE: " << std::endl;
    print_matrix(association_profile_save_means);
    std::cout << std::endl << std::endl;*/

    out << "## Starting bin avrage calculation"; out.endl();
    out << "Bin number: " << d_crit_bin_n << ", max val: " << association_profile_save[0].back(); out.endl();
    k = 0;
    for(j=0; j < association_profile_save[1].size(); j++) {
        if(association_profile_save[0][j] > (association_profile_save[0][0] + (k+1)*d_crit_bin_width)) {
            out << "-- Ending bin at " << association_profile_save_means[1][k] << ", to " << association_profile_save_means[1][k]/((double)bin_member_counter) << ", using n=" << bin_member_counter; out.endl();
            association_profile_save_means[1][k] = association_profile_save_means[1][k]/((double)bin_member_counter);
            k++;

            while(association_profile_save[0][j] > (association_profile_save[0][0] + (k+1)*d_crit_bin_width)) {
                out << "-- Skipping bin at " << (association_profile_save[0][0] + (k+1)*d_crit_bin_width); out.endl();
                association_profile_save_means[1][k] = 0;
                k++;
            }
        }

        //std::cout << "current bin val: " << association_profile_save_means[1][k] << ", j: " << j <<  ", k: " << k << std::endl;
        if(association_profile_save_means[1][k] == -1) {
            bin_member_counter = 1;
            out << "-- Starting new bin at: (" << association_profile_save[0][j] << "," << association_profile_save[1][j] << "), bin center and width: (" << association_profile_save_means[0][k] << "," << d_crit_bin_width << ")"; out.endl();
            association_profile_save_means[1][k] = association_profile_save[1][j];
        }
        else {
            out << "-- Adding value to bin: (" << association_profile_save[0][j] << "," << association_profile_save[1][j] << "), bin value: " << association_profile_save_means[1][k]; out.endl();
            bin_member_counter++;
            association_profile_save_means[1][k] = association_profile_save_means[1][k] + association_profile_save[1][j];
        }

    }
    out << "-- Ending bin at " << association_profile_save_means[1][k] << ", to " << association_profile_save_means[1][k]/((double)bin_member_counter) << ", using n=" << bin_member_counter; out.endl();
    association_profile_save_means[1][k] = association_profile_save_means[1][k]/((double)bin_member_counter);

    for(j=0; j < (d_crit_bin_n-1); j++) {
        if(association_profile_save_means[1][j] == 0 && j > 0) {
            empty_bin_counter=j;
            while(association_profile_save_means[1][empty_bin_counter] == 0) {
                empty_bin_counter++;
            }
            out << "-- Discovered empty bins, extrapolating from " << association_profile_save_means[1][j-1] << " to " << association_profile_save_means[1][empty_bin_counter] << ", missing bins=" << empty_bin_counter-j;out.endl();
            for(k=j; k < empty_bin_counter; k++) {
                association_profile_save_means[1][k] = association_profile_save_means[1][j-1] + ( ((double)(k-j+1))/((double)(empty_bin_counter-j+1)) )*(association_profile_save_means[1][empty_bin_counter] - association_profile_save_means[1][j-1]);
                out << "-- Assigning (" << association_profile_save_means[0][k] << "," << association_profile_save_means[1][k] << "), ";
            }
            out.endl();
        }
    }

    out << "- Transposing data matrix."; out.endl();
    association_profile_save_means = matrix_transpose(association_profile_save_means);

    if(file_exists(statistics_association_profile[i])) {
        remove(statistics_association_profile[i].c_str());
        out << "- Removed file: " << statistics_association_profile[i]; out.endl();
    }

    save_mat(statistics_association_profile[i], &association_profile_save_means);
    out << "- Saved file: " << statistics_association_profile[i]; out.endl();
}

clean_up(matrix_created_files,1);
clean_up(matrix_created_files,2);
clean_up(matrix_created_files,3);

if(total_showers_created >= opt.showers) {
    MC_not_done = FALSE;

    out << "########################################################################"; out.endl();
    out << "# Goal met of total_showers_created = " << total_showers_created << ", exiting loop"; out.endl();
    out << "########################################################################"; out.endl();
}
else {
    out << "# Goal NOT met, continuing loop"; out.endl();
}

} while(MC_not_done);


//out << "CLUSTER FOR ENCOUNTERS: "; out.endl();
//print_matrix(cluster_matrix);

    out.endl();
    out << "# Moving results to collected simulation folder " << data_out_folder_string; out.endl();
    
    boost::filesystem::path dir(data_out_folder_string);
    if(boost::filesystem::create_directory(dir)) {
        out << "# Folder create success"; out.endl();
    }

    out << "# Moving file " << met_encounter_data_out; out.endl();
    execute_command = data_out_folder_string + "met_encounter_data.txt";
    rename(met_encounter_data_out.c_str(), execute_command.c_str());

    out << "# Moving file " << earth_encounter_data_out; out.endl();
    execute_command = data_out_folder_string + "earth_encounter_data.txt";
    rename(earth_encounter_data_out.c_str(), execute_command.c_str());

    out << "# Moving file " << met_before_encounter_data_out; out.endl();
    execute_command = data_out_folder_string + "met_before_encounter_data.txt";
    rename(met_before_encounter_data_out.c_str(), execute_command.c_str());

    out << "# Moving file " << statistics_PB_kep_out; out.endl();
    execute_command = data_out_folder_string + "ejector_kep_data.txt";
    rename(statistics_PB_kep_out.c_str(), execute_command.c_str());

    out << "# Moving file " << statistics_PB_char_out; out.endl();
    execute_command = data_out_folder_string + "ejector_char_data.txt";
    rename(statistics_PB_char_out.c_str(), execute_command.c_str());

    out << "# Moving file " << statistics_PB_n_out; out.endl();
    execute_command = data_out_folder_string + "ejector_n_data.txt";
    rename(statistics_PB_n_out.c_str(), execute_command.c_str());

    out << "# Moving file " << PBE_weight_statistics; out.endl();
    execute_command = data_out_folder_string + "particle_weight_stat.txt";
    rename(PBE_weight_statistics.c_str(), execute_command.c_str());

    out << "# Moving file " << PBE_abs_v_statistics; out.endl();
    execute_command = data_out_folder_string + "particle_ejection_v_stat.txt";
    rename(PBE_abs_v_statistics.c_str(), execute_command.c_str());

    out << "# Moving file " << PBE_M_statistics; out.endl();
    execute_command = data_out_folder_string + "PBE_M_stat.txt";
    rename(PBE_M_statistics.c_str(), execute_command.c_str());

    for(i=0; i < FUNCTIONS.size(); i++) {
        out << "# Moving file " << statistics_association_profile[i]; out.endl();
        execute_command = data_out_folder_string + "association_profile_" + FUNCTIONS[i] + ".txt";
        rename(statistics_association_profile[i].c_str(), execute_command.c_str()); 

        out << "# Moving file " << statistics_Dc_min_err[i]; out.endl();
        execute_command = data_out_folder_string + "min_error_" + FUNCTIONS[i] + "_Dc.txt";
        rename(statistics_Dc_min_err[i].c_str(), execute_command.c_str()); 
    }

    out << "# Moving file " << dist_debug_a; out.endl();
    execute_command = data_out_folder_string + "a_dist.txt";
    rename(dist_debug_a.c_str(), execute_command.c_str());

    out << "# Moving file " << dist_debug_e; out.endl();
    execute_command = data_out_folder_string + "e_dist.txt";
    rename(dist_debug_e.c_str(), execute_command.c_str());

    out << "# Moving file " << dist_debug_i; out.endl();
    execute_command = data_out_folder_string + "i_dist.txt";
    rename(dist_debug_i.c_str(), execute_command.c_str());

    out << "# Moving file " << dist_debug_omega; out.endl();
    execute_command = data_out_folder_string + "omega_dist.txt";
    rename(dist_debug_omega.c_str(), execute_command.c_str());

    out << "# Moving file " << dist_debug_Omega; out.endl();
    execute_command = data_out_folder_string + "Omega_dist.txt";
    rename(dist_debug_Omega.c_str(), execute_command.c_str());

    out << "# Moving file " << log_file; out.endl();
    execute_command = data_out_folder_string + "outlog.txt";
    rename(log_file.c_str(), execute_command.c_str());

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
    int i,j;
    bool selector;
    if(type == -1) {
        selector = TRUE;
    }
    else {
        selector = FALSE;
    }
    for(j=0; j < matrix_created_files.size(); j++ ) {
        if(type == j || selector) {
            for (i = 0; i < matrix_created_files[j].size(); ++i) {
                if(file_exists(matrix_created_files[j][i])) {
                    remove(matrix_created_files[j][i].c_str());
                    std::cout << "Removed file: " << matrix_created_files[j][i] << std::endl;
                }
            }
        }
    }
}

