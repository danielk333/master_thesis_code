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
    int PBE_status,CMS_status,mercury6_status,JPL_status,OrbClone_status,OAA_status;
    //long long int dummy_var;
    //std::cin >> dummy_var;
    /* ITTERATORS */ 
    long long int i;
    unsigned long long int ui,uj,uk;

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

    std::ostringstream strs;
    std::string strs_str;
    //strs << NUMBER; strs_str = strs.str();

    std::string out_file_name;
    std::ofstream out_stream;
    std::string line_end = LF;

    //String for module return
    std::string MODULE_RETURN;

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
    OBSERVATION_record TP_temp2;

    //FUNCTIONS TO USE
    std::vector<std::string> FUNCTIONS;
    FUNCTIONS.push_back("SH");
    FUNCTIONS.push_back("D");
    FUNCTIONS.push_back("rho2");
    FUNCTIONS.push_back("varrho1");

    /* Generate programs self_path for for dynamics folders */
std::string selfpath = get_selfpath();
selfpath.erase(selfpath.end()-11,selfpath.end());
    
    std::string MCAS_config = "MCAS_config.cfg";
    std::string PBE_IC_body = "body_data.data";
    std::string PBE_config = "PBE_config.cfg";
    std::string CMS_config = "CMS_config.cfg";
    std::string SUOC_config = "SUOC_config.cfg";
    std::string OAA_config = "OAA_config.cfg";
    std::string OAA_config_simple = "OAA_config_simple.cfg";

    std::vector<double> OAA_config_vector;

    std::string CMS_program = "CMS/./CMS";
    std::string mercury6_folder = "mercury6/";
    std::string statistics_folder = "STATISTICS_OUTPUT/";
    std::string mercury6_program = "./mercury6";
    std::string close6_program = "./close6";
    std::string element6_program = "./element6";
    std::string PBE_program = "PBE/./PBE";
    std::string JPL_program = "JPL/./JPL";
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
    std::string SUOC_folder = "SUOC/";
    std::string SUOC_program = "./SUOC";
    std::string OAA_folder = "OAA/";
    std::string OAA_program = "./OAA";
    std::string SPICE_folder = "SPICE/";
    std::string SPICE_program = "./JPL_SPICE";

    std::string STAT_DEBUG = selfpath + "OUTPUT/DEBUG_STAT.txt";

    std::string temp_file_path_2;
    std::string TEMP_file = selfpath + "OUTPUT/TEMP_FILE.tmp";

    std::string pos_data_file = "pos.data";
    std::string vel_data_file = "vel.data";
    std::string mass_data_file = "mass.data";
    std::string prop_data_file = "prop.data";

    std::string mercury6_big = selfpath + mercury6_folder + "big.in";
    std::string mercury6_small = selfpath + mercury6_folder + "small.in";
    std::string mercury6_param = selfpath + mercury6_folder + "param.in";

    std::string mercury6_xv = selfpath + mercury6_folder + "xv.out";
    std::string mercury6_ce = selfpath + mercury6_folder + "ce.out";
    std::string mercury6_info = selfpath + mercury6_folder + "info.out";

    std::string mercury6_big_dmp = selfpath + mercury6_folder + "big.dmp";
    std::string mercury6_small_dmp = selfpath + mercury6_folder + "small.dmp";
    std::string mercury6_param_dmp = selfpath + mercury6_folder + "param.dmp";
    std::string mercury6_restart_dmp = selfpath + mercury6_folder + "restart.dmp";

    std::string mercury6_big_tmp = selfpath + mercury6_folder + "big.tmp";
    std::string mercury6_small_tmp = selfpath + mercury6_folder + "small.tmp";
    std::string mercury6_param_tmp = selfpath + mercury6_folder + "param.tmp";
    std::string mercury6_restart_tmp = selfpath + mercury6_folder + "restart.tmp";

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
    std::string PBE_M_comet_out = selfpath + PBE_OUT_folder + "body_m_data.txt";
    std::string PBE_sim_out = selfpath + PBE_OUT_folder + "sim_data.txt";

    std::string PBE_weight_statistics = selfpath + PBE_OUT_folder + "particle_weight_stat.txt";
    std::string PBE_abs_v_statistics = selfpath + PBE_OUT_folder + "particle_ejection_v_stat.txt";
    std::string PBE_all_m_statistics = selfpath + PBE_OUT_folder + "particle_ejection_m_stat.txt";
    std::string PBE_M_statistics = selfpath + PBE_OUT_folder + "PBE_M_stat.txt";

    std::string PBE_m_dist_in = selfpath + PBE_IC_PRODUCER_folder + "mass_dist.data";
    std::string PBE_mass_config = selfpath + settings_folder + "mass_dist.txt";

    std::string CMS_q_in = selfpath + CMS_MASSLESS_folder + pos_data_file;
    std::string CMS_v_in = selfpath + CMS_MASSLESS_folder + vel_data_file;
    std::string CMS_prop_in = selfpath + CMS_MASSLESS_folder + prop_data_file;

    std::string planet_m_init = selfpath + settings_folder + mass_data_file;
    std::string PBE_m_init = selfpath + JPL_OUT_folder + mass_data_file;
    std::string CMS_m_init = selfpath + CMS_init_folder + mass_data_file;

    std::string CMS_data_out = selfpath + CMS_OUT_folder + "body_data.txt";
    std::string CMS_diag_out = selfpath + CMS_OUT_folder + "energy_diag.txt";
    std::string CMS_sett_out = selfpath + CMS_OUT_folder + "sim_settings.txt";

    std::string mercury6_earth_clo = selfpath + mercury6_folder + "EARTH.clo";
    std::string mercury6_earth_aei = selfpath + mercury6_folder + "EARTH.aei";
    std::string mercury6_PB_aei = selfpath + mercury6_folder + "PB.aei";
    std::string mercury6_element_cfg = selfpath + mercury6_folder + "element.in";
    std::string mercury6_temp_obj_str;

    std::string particle_snapshot_vel = selfpath + statistics_folder + "snapshot_vel.txt";
    std::string particle_snapshot_pos = selfpath + statistics_folder + "snapshot_pos.txt";
    std::string particle_snapshot_earth = selfpath + statistics_folder + "snapshot_earth.txt";

    std::string met_encounter_data_out = selfpath + statistics_folder + "met_encounter_data.txt";
    std::string earth_encounter_data_out = selfpath + statistics_folder + "earth_encounter_data.txt";
    std::string met_before_encounter_data_out = selfpath + statistics_folder + "met_before_encounter_data.txt";
    std::string statistics_PB_kep_out = selfpath + statistics_folder + "ejector_kep_data.txt";
    std::string statistics_all_PB_kep_out = selfpath + statistics_folder + "PB_kep_data.txt";
    std::string statistics_time_series_PB_kep_out = selfpath + statistics_folder + "PB_kep_timeS_data.txt";
    std::string statistics_PB_char_out = selfpath + statistics_folder + "ejector_char_data.txt";
    std::string statistics_PB_n_out = selfpath + statistics_folder + "ejector_n_data.txt";
    std::string simulation_time_data = selfpath + statistics_folder + "time_data.txt";

    std::string OAA_input_file = selfpath + OAA_folder + "OAA_input.txt";
    std::string OAA_arguments = selfpath + settings_folder + OAA_config_simple + " " + OAA_input_file + " 1";
    
    std::vector<std::string> OAA_OUTPUT_D_mat_folders; 
    std::vector<std::string> OAA_OUTPUT_D_mat_files;
    std::vector<std::string> OAA_OUTPUT_profile_files;
    std::vector<std::string> OAA_OUTPUT_error_files;
    std::vector<std::string> OAA_OUTPUT_cluster_files;
    for(ui=0; ui < FUNCTIONS.size(); ui++) {
        OAA_OUTPUT_D_mat_folders.push_back(selfpath + OAA_folder + "OUTPUT_" + FUNCTIONS[ui]);
        OAA_OUTPUT_D_mat_files.push_back(selfpath + OAA_folder + "OUTPUT_" + FUNCTIONS[ui] + "_mat.txt");
        OAA_OUTPUT_profile_files.push_back(selfpath + OAA_folder + "OUTPUT_" + FUNCTIONS[ui] + "_association_profile.txt");
        OAA_OUTPUT_error_files.push_back(selfpath + OAA_folder + "OUTPUT_" + FUNCTIONS[ui] + "_error_profile.txt");
        OAA_OUTPUT_cluster_files.push_back(selfpath + OAA_folder + "OUTPUT_" + FUNCTIONS[ui] + "_clusters.txt");
    }

    std::vector<std::string> statistics_association_profile;
    std::vector<std::string> statistics_error_profile;
    std::vector<std::string> statistics_function_divergance;
    std::vector<std::string> stream_dissipation;
    for(ui=0; ui < FUNCTIONS.size(); ui++) {
        statistics_association_profile.push_back(selfpath + statistics_folder + "association_profile_" + FUNCTIONS[ui] + ".txt");
        statistics_error_profile.push_back(selfpath + statistics_folder + "error_profile_" + FUNCTIONS[ui] + ".txt");
        statistics_function_divergance.push_back(selfpath + statistics_folder + "function_divergance_" + FUNCTIONS[ui] + ".txt");
        stream_dissipation.push_back(selfpath + statistics_folder + "stream_dissipation_data_" + FUNCTIONS[ui] + ".txt");
    }

    std::string PB_kepler_dist =  selfpath + debug_folder + "PB_kep_dist.txt";
    std::string DEBUG_PBE_data_q_out = selfpath + PBE_OUT_folder + "debug_pos_data.txt";
    std::string DEBUG_PBE_energy_out = selfpath + PBE_OUT_folder + "debug_energy_data.txt";

    std::string OrbClone_file =  selfpath + SUOC_folder + "orb_clones.txt";

    std::string log_file = selfpath + "OUTPUT/outlog.txt";

    std::vector<std::vector<double> > EXEC_TIME_ALL;
    for(i = 0; i < 10; i++) {
        EXEC_TIME_ALL.push_back(EMPTY_DOUBLE_VEC);
    }
    clock_t tStart;
    clock_t time0 = clock_t();
    
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
    if(data_out_folder[8] == ' ') {
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
    list_created_files.push_back(PBE_m_init);
    list_created_files.push_back(CMS_m_init);
    list_created_files.push_back(JPL_CMS_pos);
    list_created_files.push_back(JPL_CMS_vel);
    list_created_files.push_back(PBE_m_dist_in);
    list_created_files.push_back(OrbClone_file);
    matrix_created_files.push_back(list_created_files); //persistant_init_files 0
    list_created_files.clear();

    list_created_files.push_back(CMS_q_in);
    list_created_files.push_back(CMS_v_in);
    list_created_files.push_back(CMS_prop_in);
    list_created_files.push_back(mercury6_small);
    list_created_files.push_back(mercury6_big);
    list_created_files.push_back(mercury6_param);
    matrix_created_files.push_back(list_created_files); //massless_init_files 1
    list_created_files.clear();

    //list_created_files.push_back("");

    list_created_files.push_back(PBE_q_out);
    list_created_files.push_back(PBE_v_out);
    list_created_files.push_back(PBE_m_out);
    list_created_files.push_back(PBE_t_out);
    list_created_files.push_back(PBE_M_comet_out);
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
    list_created_files.push_back(mercury6_PB_aei);
    for(ui=0; ui < FUNCTIONS.size(); ui++) {
        list_created_files.push_back(OAA_OUTPUT_D_mat_files[ui]);
        list_created_files.push_back(OAA_OUTPUT_D_mat_folders[ui]);
        list_created_files.push_back(OAA_OUTPUT_profile_files[ui]);
        list_created_files.push_back(OAA_OUTPUT_error_files[ui]);
        list_created_files.push_back(OAA_OUTPUT_cluster_files[ui]);
    }
    list_created_files.push_back(OAA_input_file);
    matrix_created_files.push_back(list_created_files); //output_files 2
    list_created_files.clear();

    list_created_files.push_back(mercury6_big_dmp);
    list_created_files.push_back(mercury6_small_dmp);
    list_created_files.push_back(mercury6_param_dmp);
    list_created_files.push_back(mercury6_restart_dmp);
    list_created_files.push_back(mercury6_big_tmp);
    list_created_files.push_back(mercury6_small_tmp);
    list_created_files.push_back(mercury6_param_tmp);
    list_created_files.push_back(mercury6_restart_tmp);
    list_created_files.push_back(DEBUG_PBE_data_q_out);
    list_created_files.push_back(DEBUG_PBE_energy_out);
    list_created_files.push_back(TEMP_file);
    matrix_created_files.push_back(list_created_files); //temp_files 3
    list_created_files.clear();

    list_created_files.push_back(met_encounter_data_out);
    list_created_files.push_back(earth_encounter_data_out);
    list_created_files.push_back(met_before_encounter_data_out);
    list_created_files.push_back(statistics_time_series_PB_kep_out);
    list_created_files.push_back(statistics_PB_kep_out);
    list_created_files.push_back(statistics_all_PB_kep_out);
    list_created_files.push_back(statistics_PB_char_out);
    list_created_files.push_back(statistics_PB_n_out);
    list_created_files.push_back(PBE_weight_statistics);
    list_created_files.push_back(PBE_abs_v_statistics);
    list_created_files.push_back(PBE_M_statistics);
    list_created_files.push_back(simulation_time_data);
    list_created_files.push_back(particle_snapshot_vel);
    list_created_files.push_back(particle_snapshot_pos);
    list_created_files.push_back(particle_snapshot_earth);
    list_created_files.push_back(log_file);
    for(ui=0; ui < FUNCTIONS.size(); ui++) {
        list_created_files.push_back(statistics_association_profile[ui]);
        list_created_files.push_back(statistics_error_profile[ui]);
        list_created_files.push_back(statistics_function_divergance[ui]);
        list_created_files.push_back(stream_dissipation[ui]);
    }
    matrix_created_files.push_back(list_created_files); //Output restults 4
    list_created_files.clear();

    list_created_files.push_back(PB_kepler_dist);
    list_created_files.push_back(STAT_DEBUG);
    matrix_created_files.push_back(list_created_files); //DEBUG DATA 5
    list_created_files.clear();

    //###########
    //CLEAN UP:
    // - IN CASE LAST RUN DIDNT EXIT PROPERLY
    // - TO MAKE ROOM FOR NEW RESULTS
    //###########
    clean_up(matrix_created_files,-1);

    create_file(TEMP_file);

     // INIT OUTPUT
    log_out.open(log_file.c_str(), std::ios::app);

    output_stream out(std::cout, log_out);
    out.type = opt.logfile;

    out << "Output type " << opt.logfile << " selected."; out.endl(); out.endl();

    std::vector<std::string> mercury6_element_bodies;
    out << "CHECKING FOR OLD AIE FILES"; out.endl();
mercury6_element_bodies.clear();
mercury6_element_bodies.push_back("MERCURY");
mercury6_element_bodies.push_back("VENUS");
mercury6_element_bodies.push_back("EARTH");
mercury6_element_bodies.push_back("MARS");
    mercury6_element_bodies.push_back("JUPITER");
    mercury6_element_bodies.push_back("SATURN");
    mercury6_element_bodies.push_back("URANUS");
    mercury6_element_bodies.push_back("NEPTUNE");

    for(ui=0; ui < 9999; ui++) {
    mercury6_temp_obj_str = selfpath + mercury6_folder + "OBJ" + IntToString(ui) + ".aei";
    if(file_exists(mercury6_temp_obj_str)) {
        remove(mercury6_temp_obj_str.c_str());
        out << "# Removed file: " << mercury6_temp_obj_str; out.endl();
    }
    }
    for(ui=0; ui < mercury6_element_bodies.size(); ui++) {
    mercury6_temp_obj_str = selfpath + mercury6_folder + mercury6_element_bodies[ui] + ".aei";
    if(file_exists(mercury6_temp_obj_str)) {
        remove(mercury6_temp_obj_str.c_str());
        out << "# Removed file: " << mercury6_temp_obj_str; out.endl();
    }
    }


    std::vector<std::vector<double> > orbital_dist_data;

    //###########
    //DATABASE STORAGE
    //###########

    std::vector<std::vector<double> > MURMHED_db;
    std::vector<std::vector<double> > MURMHED_db_distance_matrix_part;

    std::vector<std::vector<double> > time_corrected_vectors_obj,time_corrected_vectors_earth;
    std::vector<double> temp_kepler,temp_state_vec,x_temp,v_temp;
    double snapshot_offset;
    unsigned long long int propagation_counter;
    unsigned long long int TEMP_ID_index;
    std::vector<double> temp_obj_pos_a;
    std::vector<double> temp_obj_pos_b;

    std::vector<double> temp_obj_vel_a;
    std::vector<double> temp_obj_vel_b;
    double mean_stream_dissipation;
    unsigned long long int internal_dissipation_counter;
    std::vector<unsigned long long int> internal_dissipation_vector;
    double XNp1;
    double std_stream_dissipation;

    std::vector<std::pair<size_t, double_vector_iterator> > order_D_c_vec;

    std::vector<std::vector<double> > body_norm_dist;
    std::vector<std::vector<double> > body_norm_dist_bins;
    unsigned long long int max_particle_number_mercury;
    unsigned long long int particle_number_mercury;
    unsigned long long int used_particle_number_mercury;

    std::vector<std::vector<double> > temp_q;
    std::vector<std::vector<double> > temp_v;

    std::vector<std::vector<double> > temp_q_major; 
    std::vector<std::vector<double> > temp_v_major; 
    std::vector<double> temp_m_major; 

    std::vector<double> time_vec;

    std::vector<std::vector<double> > snapshot_data_pos;
    std::vector<std::vector<double> > snapshot_data_vel;
    std::vector<std::vector<double> > snapshot_data_earth;

    std::vector<std::vector<double> > temp_t;
    std::vector<std::vector<double> > temp_sim;
    std::vector<std::vector<double> > temp_m_mat;
    std::vector<double> temp_m;
    std::vector<double> PBE_particle_m;
    long long int earth_snap_range;
    std::vector<double> temp_id_vec;
    
    std::vector<double> used_t;
    std::vector<double> temp_row;
    std::vector<double> temp_data;

    std::vector<std::vector<double> > temp_obj_data;
    std::vector<std::vector<double> > temp_obj_data2;
    std::vector<std::vector<double> > temp_earth_data;

    std::vector<std::vector<double> > temp_mat;
    std::vector<std::vector<double> > temp_mat2;
    std::vector<double> temp_vec;

    double rand_orb;
    unsigned long long int rand_orb_ind;
    std::vector<std::vector<unsigned int> > index_matrix;
    std::vector<std::vector<unsigned int> > cluster_matrix_next_db;
    std::string buffer_temp;

    std::vector<unsigned int> EMPTY_UINT_VEC;
    /*double J_date;
    double J_date_diff;
    double theta_earth;*/

    double total_PBE_sim_time; 
    double TP_diff;
    double PB_period;
    double PBE_time_sync;
    double mercury6_start_time;
    unsigned long long int total_particles_integrated;
    unsigned long long int total_close_encouters_investigated;

    std::vector<double> current_body_xv;

    std::vector<double> association_vec;
    std::vector<double> d_crit_vec;
    double mercury6_stop_time;

    bool TRUTH_CHECKER;

    std::string temp_batch_name;
    std::string temp_system_command;

    std::vector<std::string> list_of_batches_d_matrix;

    //###########
    //MODEL STORAGE
    //###########
    std::vector<double> grun_model_mass_flux;

    //###########
    //RESULT STORAGE
    //###########
    unsigned int period_n;

    std::vector<std::vector<double> > encounter_data;
    std::vector<std::vector<double> > PB_time_series_data;
    std::vector<std::vector<double> > earth_encounter_data;
    std::vector<std::vector<double> > encounter_stat_data;

    std::vector<std::vector<double> > PB_kep_dist_data;

    std::vector<std::vector<double> > D_SH_stream_dissipation_time_series;
    std::vector<std::vector<double> > D_D_stream_dissipation_time_series;
    std::vector<std::vector<double> > D_rho2_stream_dissipation_time_series;
    std::vector<std::vector<double> > D_varrho1_stream_dissipation_time_series;

    std::vector<std::vector<double> > D_SH_diff_mat_time_series;
    std::vector<std::vector<double> > D_D_diff_mat_time_series;
    std::vector<std::vector<double> > D_rho2_diff_mat_time_series;
    std::vector<std::vector<double> > D_varrho1_diff_mat_time_series;

    std::vector<double> error_function;
    std::vector<std::vector<double> > D_matix;
    std::vector<std::vector<unsigned int> > cluster_matrix;
    std::vector<double> fractions;
    std::vector<std::vector<unsigned int> > cluster_matrix_next;
    std::vector<double> fractions_next;
    std::vector<std::vector<double> > association_profile;
    for(ui=0; ui < FUNCTIONS.size(); ui++) {
        association_profile.push_back(EMPTY_DOUBLE_VEC); //dcrit
        association_profile.push_back(EMPTY_DOUBLE_VEC); //fractions
    }

    std::vector<std::vector<double> > association_profile_save;
    std::vector<std::vector<double> > association_profile_save_means;

    std::vector<std::vector<double> > association_profile_db;
    for(ui=0; ui < FUNCTIONS.size(); ui++) {
        association_profile_db.push_back(EMPTY_DOUBLE_VEC); //dcrit
        association_profile_db.push_back(EMPTY_DOUBLE_VEC); //fractions
    }

    std::vector<std::vector<double> > association_profile_db_save;
    std::vector<std::vector<double> > association_profile_db_save_means;

    std::vector<std::vector<double> > filtered_encounter_data;
    std::vector<std::vector<double> > filtered_encounter_data_save;
    std::vector<unsigned int> id_vec;
    std::vector<unsigned int> orig_id_vec;

    unsigned long long int first_id;
    std::vector<std::vector<double> > body_family;

    //###########
    //Counters
    //###########
    unsigned long long int total_n_created, total_PB_created, total_met_collided, current_n_created, total_showers_created,total_X_created;

    //###########
    //MISC
    //###########
    bool orbit_ok;
    double body_mass;
    std::vector<double>::iterator TEMP_ID;

    //###########
    //Physical params
    //###########
    //double FUDGE = 1.05;
    double earth_hill_r = 1.00000261*(1-0.01671123)*cbrt(5.97219e24/(3*1.98855e30));
    std::string unit;

    //double ecliptic_inc = 1.5787*PI/180.0;
    //double ecliptic_lon = 107.5822*PI/180.0;

    out << "# earth_hill_r    : " << earth_hill_r; out.endl();
    //out << "# ecliptic_inc    : " << ecliptic_inc; out.endl();
    //out << "# ecliptic_lon    : " << ecliptic_lon; out.endl();

    /*######################################################################## */
    // LOAD DATA
    /*######################################################################## */

   

    /*######################################################################## */
    // CALCULATE ORBITAL PDF's
    /*######################################################################## */

    // ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤


    /* LOAD OPTIONS */
    std::string config_path = selfpath + settings_folder + MCAS_config;
    out << "Loading " << config_path; out.endl();
    error(load_file(&opt,config_path,DATA_OPTIONS),matrix_created_files);

    double time_j2000_diff;

    /*######################################################################## */
    // CALCULATE SOLARSYSTEM IC
    /*######################################################################## */
    std::vector<std::string> BODIES_NAMES;
    BODIES_NAMES.push_back("SUN CENTER");
    BODIES_NAMES.push_back("MERCURY CENTER");
    BODIES_NAMES.push_back("VENUS CENTER");
    BODIES_NAMES.push_back("EARTH CENTER");
    BODIES_NAMES.push_back("MARS BARYCENTER");
    BODIES_NAMES.push_back("JUPITER BARYCENTER");
    BODIES_NAMES.push_back("SATURN BARYCENTER");
    BODIES_NAMES.push_back("URANUS BARYCENTER");
    BODIES_NAMES.push_back("NEPTUNE BARYCENTER");

    std::string FILES_LIST;
    std::string LEAPSECOND_FILE;
    std::vector<std::string> EPHEMERIS_FILES;
    LEAPSECOND_FILE = selfpath + SPICE_folder + "naif0011.tls";
    EPHEMERIS_FILES.push_back(selfpath + SPICE_folder + "de430.bsp");
    EPHEMERIS_FILES.push_back(selfpath + SPICE_folder + "de431_part-1.bsp");
    EPHEMERIS_FILES.push_back(selfpath + SPICE_folder + "de431_part-2.bsp");
    EPHEMERIS_FILES.push_back(selfpath + SPICE_folder + "de430.bsp");

    for(i=0; i < EPHEMERIS_FILES.size(); i++) {
        FILES_LIST = FILES_LIST + EPHEMERIS_FILES[i] + " ";
    }
    FILES_LIST = FILES_LIST + LEAPSECOND_FILE;


    temp_file_path = selfpath + MURMHED_db_file;
    error(load_file(&MURMHED_db,temp_file_path,DATA_MATRIX),matrix_created_files);

    std::vector<std::vector<double> > DATA_ORB;
    std::vector<double> DATA_MJD;
    std::vector<double> DATA_JD;
    std::vector<double> temp_orb;
    std::vector<std::vector<double> > DATA_ID;
    double MJD_mod = 2400000.5;//
    double j2000_mod = 2451545.0; //
    out << "Extracting orbital elements"; out.endl();
    unsigned int DATA_ADD;
    unsigned int DATA_START;
    if(opt.MODE == 0) {
        DATA_ADD = MURMHED_db.size();
        DATA_START = opt.dist_type;
    }
    else {
        DATA_ADD = opt.MODE;
        DATA_START = opt.dist_type;
    }



    for(i=DATA_START; i < (DATA_START+DATA_ADD); i++) { //MURMHED_db.size()
        if(MURMHED_db[i][27] < 1.0) {
            //out << "MJD: " << MURMHED_db[i][1] << ", Y: " << MURMHED_db[i][2]; out.endl();
            DATA_MJD.push_back(MURMHED_db[i][1]);
            DATA_JD.push_back(MURMHED_db[i][1] + MJD_mod);
            DATA_ID.push_back(EMPTY_DOUBLE_VEC);
            DATA_ID.back().push_back(MURMHED_db[i][0]);
            
            DATA_ORB.push_back(EMPTY_DOUBLE_VEC);
            DATA_ORB.back().push_back(MURMHED_db[i][26]);
            DATA_ORB.back().push_back(MURMHED_db[i][27]);
            DATA_ORB.back().push_back(MURMHED_db[i][30]);
            DATA_ORB.back().push_back(MURMHED_db[i][31]);
            DATA_ORB.back().push_back(MURMHED_db[i][29]);
        }
    }
    out << "# " << DATA_ORB.size() << " orbits added, " << (MURMHED_db.size() - DATA_ORB.size()) << " orbits excluded"; out.endl();
    MURMHED_db.clear();
    double MAX_MJD,MIN_MJD;
    MAX_MJD = *std::max_element(DATA_MJD.begin(),DATA_MJD.end());
    MIN_MJD = *std::min_element(DATA_MJD.begin(),DATA_MJD.end());

    mercury6_start_time = MAX_MJD + MJD_mod;
    mercury6_stop_time = mercury6_start_time - opt.integration_time*365.25;
    time_j2000_diff = (mercury6_start_time - j2000_mod)*(3600.0*24.0); //sec

            strs << time_j2000_diff; strs_str = strs.str();

            for(i=0; i < BODIES_NAMES.size(); i++) {
                execute_command = selfpath + SPICE_folder + SPICE_program + " " + BODIES_NAMES[i] + " " + strs_str + " " + JPL_OUT_folder + " " + "NONE" + " " + FILES_LIST;
                out << "Calling " << execute_command; out.endl();
                out.endl(); out << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" ; out.endl();
                MODULE_RETURN = exec(execute_command.c_str());
                out << MODULE_RETURN;
                out << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" ; out.endl();out.endl(); 
            }

    

    temp_q_major.clear(); temp_v_major.clear(); temp_m_major.clear();
    out << "# Loaded Planet position data"; out.endl();
    error(load_file(&temp_q_major,JPL_PBE_pos,DATA_MATRIX),matrix_created_files);
    out << "# Loaded Planet velocity data"; out.endl();
    error(load_file(&temp_v_major,JPL_PBE_vel,DATA_MATRIX),matrix_created_files);
    out << "# Loaded Planet mass data"; out.endl();
    error(load_file(&temp_m_major,planet_m_init,DATA_VECTOR),matrix_created_files);


            out << "# Converting mass to kg"; out.endl();
            out << "--------- AU**2/D**3 -> kg ----------"; out.endl();
            for(i=0; i<temp_m_major.size(); i++) {
                out << "MASS " << i << ": " << temp_m_major[i] << " -> "; 
                temp_m_major[i] = temp_m_major[i]/pow(G_gauss,2)*Msol;
                out << temp_m_major[i]; out.endl();
            }


            out << "# Removing old Planet state data"; out.endl();
            if(file_exists(mercury6_big)) {
                remove(mercury6_big.c_str());
                out << "Removed file: " << mercury6_big; out.endl();
            }
            out << "# Print big planet data to mercury6"; out.endl();
            error(print_mercury6_big(mercury6_big,&temp_q_major,&temp_v_major,&temp_m_major,mercury6_start_time),matrix_created_files);
            out << "Saved file: " << mercury6_big; out.endl();
           
temp_q.clear();temp_v.clear();temp_m.clear();
            for(i=0; i<DATA_ORB.size(); i++) {
                temp_orb = DATA_ORB[i];
                temp_orb[0] = temp_orb[0]*AU;
                temp_orb[2] = temp_orb[2]*(PI/180.0);
                temp_orb[3] = temp_orb[3]*(PI/180.0);
                temp_orb[4] = temp_orb[4]*(PI/180.0);
                temp_orb.push_back(-temp_orb[3]);
                temp_vec = kepler_to_xv(temp_orb,Msol*G);

                temp_q.push_back(EMPTY_DOUBLE_VEC);
                temp_q.back().push_back(temp_vec[0]);
                temp_q.back().push_back(temp_vec[1]);
                temp_q.back().push_back(temp_vec[2]);

                temp_v.push_back(EMPTY_DOUBLE_VEC);
                temp_v.back().push_back(temp_vec[3]);
                temp_v.back().push_back(temp_vec[4]);
                temp_v.back().push_back(temp_vec[5]);

                temp_m.push_back(0.0);
            }


if(file_exists(mercury6_param)) {
                remove(mercury6_param.c_str());
                out << "Removed file: " << mercury6_param; out.endl();
            }
            error(print_mercury6_param(mercury6_param, &opt, mercury6_start_time),matrix_created_files);
            out << "Saved file: " << mercury6_param; out.endl();
            
            out << "# mercury6 start : " << JDtoJ(mercury6_start_time) << " y"; out.endl();
            out << "# mercury6 stop  : " << JDtoJ(mercury6_stop_time) << " y"; out.endl();
            out << "# mercury6 start : " << (mercury6_start_time) << " JD"; out.endl();
            out << "# mercury6 stop  : " << (mercury6_stop_time) << " JD"; out.endl();
            out << "# mercury6 diff  : " << (mercury6_stop_time-mercury6_start_time) << " days"; out.endl();
max_particle_number_mercury = 10000;
            out << "# Building Small-body initial data: " << mercury6_small; out.endl();
            if(file_exists(mercury6_small)) {
                remove(mercury6_small.c_str());
                out << "Removed file: " << mercury6_small; out.endl();
            }
            error(print_mercury6_small(mercury6_small, &temp_q, &temp_v, &temp_m, &DATA_JD, 0.0, max_particle_number_mercury, &particle_number_mercury),matrix_created_files);
            out << "Saved file: " << mercury6_small; out.endl();
            temp_q.clear(); temp_v.clear(); temp_m.clear();

            execute_command = "cd " + selfpath + mercury6_folder + "; " + mercury6_program;
            out << "# Calling " << execute_command; out.endl();
            
        /*out.endl(); out << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" ; out.endl();
        MODULE_RETURN = exec(execute_command.c_str());
        out << MODULE_RETURN;
        out << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" ; out.endl();out.endl(); */
        system(execute_command.c_str());



//Remove config file for rewriting
if(file_exists(mercury6_element_cfg)) {
    remove(mercury6_element_cfg.c_str());
    out << "Removed file: " << mercury6_element_cfg; out.endl();
}

out << "# Building Small-body element extract data file: " << mercury6_element_cfg; out.endl();

mercury6_element_bodies.clear();
mercury6_element_bodies.push_back("MERCURY");
mercury6_element_bodies.push_back("VENUS");
mercury6_element_bodies.push_back("EARTH");
mercury6_element_bodies.push_back("MARS");
    mercury6_element_bodies.push_back("JUPITER");
    mercury6_element_bodies.push_back("SATURN");
    mercury6_element_bodies.push_back("URANUS");
    mercury6_element_bodies.push_back("NEPTUNE");

    for(ui=0; ui < DATA_ORB.size(); ui++) {
        mercury6_element_bodies.push_back("OBJ" + IntToString(ui));
    }


error(print_mercury6_element(mercury6_element_cfg,mercury6_element_bodies,0),matrix_created_files);

out << "Created file: " << mercury6_element_cfg; out.endl();

execute_command = "cd " + selfpath + mercury6_folder + "; " + element6_program;
out << "Calling " << execute_command; out.endl();

out.endl(); out << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" ; out.endl();
        MODULE_RETURN = exec(execute_command.c_str());
        out << MODULE_RETURN;
        out << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" ; out.endl();out.endl();

std::vector<std::vector<double> > JUPITER_data;
std::vector<std::vector<double> > T_param;
/*
std::vector<std::vector<double> > TEMP_PLANET_DATA1;
std::vector<std::vector<double> > TEMP_PLANET_DATA2;
std::vector<std::vector<double> > TEMP_PLANET_DATA3;
std::vector<std::vector<double> > TEMP_PLANET_DATA4;
std::vector<std::vector<double> > TEMP_PLANET_DATA5;
std::vector<std::vector<double> > TEMP_PLANET_DATA6;

mercury6_temp_obj_str = selfpath + mercury6_folder + mercury6_element_bodies[0] + ".aei";
error(read_mercury6_matrix(&TEMP_PLANET_DATA1,mercury6_temp_obj_str),matrix_created_files);
mercury6_temp_obj_str = selfpath + mercury6_folder + mercury6_element_bodies[1] + ".aei";
error(read_mercury6_matrix(&TEMP_PLANET_DATA2,mercury6_temp_obj_str),matrix_created_files);
mercury6_temp_obj_str = selfpath + mercury6_folder + mercury6_element_bodies[3] + ".aei";
error(read_mercury6_matrix(&TEMP_PLANET_DATA3,mercury6_temp_obj_str),matrix_created_files);
mercury6_temp_obj_str = selfpath + mercury6_folder + mercury6_element_bodies[5] + ".aei";
error(read_mercury6_matrix(&TEMP_PLANET_DATA4,mercury6_temp_obj_str),matrix_created_files);
mercury6_temp_obj_str = selfpath + mercury6_folder + mercury6_element_bodies[6] + ".aei";
error(read_mercury6_matrix(&TEMP_PLANET_DATA5,mercury6_temp_obj_str),matrix_created_files);
mercury6_temp_obj_str = selfpath + mercury6_folder + mercury6_element_bodies[7] + ".aei";
error(read_mercury6_matrix(&TEMP_PLANET_DATA6,mercury6_temp_obj_str),matrix_created_files);

*/

double TEMP_TVAL;

mercury6_temp_obj_str = selfpath + mercury6_folder + "JUPITER" + ".aei";
error(read_mercury6_matrix(&JUPITER_data,mercury6_temp_obj_str),matrix_created_files);
            T_param.push_back(EMPTY_DOUBLE_VEC);
            for(uj=0; uj < JUPITER_data.size(); uj++) {
                T_param.back().push_back(JUPITER_data[uj][0]);
            }

    for(ui=0; ui < DATA_ORB.size(); ui++) {
        temp_obj_data.clear();
        mercury6_temp_obj_str = selfpath + mercury6_folder + "OBJ" + IntToString(ui) + ".aei";
        //out << "# Attempting to load: " << mercury6_temp_obj_str; out.endl();
        error(read_mercury6_matrix(&temp_obj_data,mercury6_temp_obj_str),matrix_created_files);



        T_param.push_back(EMPTY_DOUBLE_VEC);
        for(uj=0; uj < JUPITER_data.size(); uj++) {
            if(uj >=  temp_obj_data.size()) {
                T_param.back().push_back( 0.0 );
            }
            else {
                TEMP_TVAL = JUPITER_data[uj][1]/temp_obj_data[uj][1] + 2.0*sqrt((temp_obj_data[uj][1]/JUPITER_data[uj][1])*(1 - pow(temp_obj_data[uj][2],2))) *cos(temp_obj_data[uj][3]*(PI/180.0));
                if(TEMP_TVAL < 0) {
                    TEMP_TVAL = 0.0;
                }
                T_param.back().push_back( TEMP_TVAL );
            }
            
        }

    }
/*
        mercury6_temp_obj_str = selfpath + mercury6_folder + "OBJ0.aei";
        //out << "# Attempting to load: " << mercury6_temp_obj_str; out.endl();
        error(read_mercury6_matrix(&temp_obj_data,mercury6_temp_obj_str),matrix_created_files);
        
save_mat("debug1.txt", &temp_obj_data);

save_mat("JUP_data.txt", &JUPITER_data);
*/
std::string temp_str2;
temp_str2 = "TISSERAND_OUT_" + IntToString(DATA_START) + "_to_" + IntToString(DATA_START+DATA_ADD) + ".txt";
save_mat(temp_str2, &T_param);
temp_str2 = "DATA_ID_OUT_" + IntToString(DATA_START) + "_to_" + IntToString(DATA_START+DATA_ADD) + ".txt";
save_mat(temp_str2, &DATA_ID);


    return CALC_OK;
}


void error(int SIM_STATUS, std::vector<std::vector<std::string> > matrix_created_files) {

  if(SIM_STATUS != CALC_OK) {
    std::cout << std::endl << "ERROR ENCOUTERD: " << SIM_STATUS << std::endl << std::endl;

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
                    std::cout << "Removed folder: " << matrix_created_files[j][i] << std::endl;
                }
            }
        }
    }
}

