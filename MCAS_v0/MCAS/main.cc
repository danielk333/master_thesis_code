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
selfpath.erase(selfpath.end()-4,selfpath.end());
    
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

    /* LOAD OPTIONS */
    std::string config_path = selfpath + settings_folder + MCAS_config;
    error(load_file(&opt,config_path,DATA_OPTIONS),matrix_created_files);

    // INIT OUTPUT
    log_out.open(log_file.c_str(), std::ios::app);

    output_stream out(std::cout, log_out);
    out.type = opt.logfile;

    out << "Output type " << opt.logfile << " selected."; out.endl(); out.endl();

    out.endl();
    out << "------- LOADED OPTIONS ----------"; out.endl();
    out << "MODE             :" << opt.MODE; out.endl();
    out << "integrator       :" << opt.integrator; out.endl();
    out << "IC_gen           :" << opt.IC_gen; out.endl();
    out << "start_date       :" << opt.start_date; out.endl();
    out << "integration_time :" << opt.integration_time; out.endl();
    out << "end_date         :" << opt.end_date; out.endl();
    out << "family           :" << opt.family; out.endl();
    out << "dist_type        :" << opt.dist_type; out.endl();
    out << "calc_clusters    :" << opt.calc_clusters; out.endl();
    out << "mass_dist        :" << opt.mass_dist; out.endl();
    out << "mass_dist_res    :" << opt.mass_dist_res; out.endl();
    out << "mass_min         :" << opt.mass_min; out.endl();
    out << "mass_max         :" << opt.mass_max; out.endl();
    out << "mem_allowed      :" << opt.mem_allowed; out.endl();
    out << "X_type           :" << opt.X_type; out.endl();
    out << "X                :" << opt.X; out.endl();
    out << "close_match      :" << opt.close_match; out.endl();
    out << "mass format      :" << opt.mass_format; out.endl();
    out << "logfile          :" << opt.logfile; out.endl();
    out << "merc_int         :" << opt.merc_int; out.endl();
    out << "merc_step        :" << opt.merc_step; out.endl();
    out << "merc_pr          :" << opt.merc_pr; out.endl();
    out << "merc_eject       :" << opt.merc_eject; out.endl();
    out << "merc_close_int   :" << opt.merc_close_int; out.endl();
    out << "mercury6_max     :" << opt.mercury6_max; out.endl();
    out << "histogram_n      :" << opt.histogram_n; out.endl();
    out << "SUOC_sigma       :" << opt.SUOC_sigma; out.endl();
    out << "SUOC_hist_res    :" << opt.SUOC_hist_res; out.endl();
    out.endl();

    out << "CHECKING FOR OLD AIE FILES"; out.endl();
    for(ui=0; ui < 9999; ui++) {
    mercury6_temp_obj_str = selfpath + mercury6_folder + "OBJ" + IntToString(ui) + ".aei";
    if(file_exists(mercury6_temp_obj_str)) {
        remove(mercury6_temp_obj_str.c_str());
        out << "# Removed file: " << mercury6_temp_obj_str; out.endl();
    }
    }

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
    std::vector<std::string> mercury6_element_bodies;

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

    temp_file_path = selfpath + settings_folder + family_data[opt.family];
    error(load_file(&body_family,temp_file_path,DATA_FAMILY),matrix_created_files);
    out.endl();
    out << "# PARENT BODY DISTRIBUTION: "; out.endl();
    out << "Distribution type     | Parameter 1, Parameter 2 |"; out.endl();
    out << "----------------------|--------------------------|"; out.endl();
    for(ui=0; ui < body_family.size(); ui++) {
        switch((unsigned int)(body_family[ui][0])) {
            case 0:
                out << "Uniform distribution  |";
                break;
            case 1:
                out << "Normal distribution   |";
                break;
        }
        switch(ui) {
            case 0:
                unit = "m";
                break;
            case 1:
                unit = "kg/m^3";
                break;
            case 2:
                unit = "";
                break;
            case 3:
                unit = "AU";
                break;
            case 4:
                unit = "kg";
                break;
        }
        
        out << body_family[ui][1] << unit << ", " << body_family[ui][2] << unit; out.endl();
    }
     out.endl();
    /*out << body_family[1][0] << ": " << body_family[1][1] << " " <<  body_family[1][2]; out.endl();
    out << body_family[2][0] << ": " << body_family[2][1] << " " <<  body_family[2][2]; out.endl();
    out << body_family[3][0] << ": " << body_family[3][1] << " " <<  body_family[3][2]; out.endl();*/

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
        temp_mat.clear(); temp_vec.clear();
        temp_vec.push_back(opt.SUOC_sigma);
        temp_vec.push_back((double)(opt.SUOC_hist_res));

        temp_mat.push_back(temp_vec);
        temp_file_path_2 = SUOC_folder + SUOC_config;
        if(file_exists(temp_file_path_2)) {
            remove(temp_file_path_2.c_str());
            out << "Removed file: " << temp_file_path_2; out.endl();
        }
        save_mat(temp_file_path_2, &temp_mat);
        out << "Saved input file: " << temp_file_path_2; out.endl(); 
        out << "Genrating " << opt.histogram_n << " orbital clones"; out.endl();
        execute_command = selfpath + SUOC_folder + SUOC_program + " kep " + IntToString(opt.histogram_n) + " kep";
        out << "Calling " << execute_command; out.endl();
        out.endl(); out << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" ; out.endl();
        MODULE_RETURN = exec(execute_command.c_str());
        out << MODULE_RETURN;
        out << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" ;
        out.endl();
        out.endl();
        //OrbClone_status = system(execute_command.c_str());
        //out << "SUOC exit status: " << OrbClone_status; out.endl();

        error(load_file(&orbital_dist_data,OrbClone_file,DATA_MATRIX),matrix_created_files);
    }

    /*######################################################################## */
    // CREATE APPROXIAMTE DUST MASS SIZE DISTRIBUTION
    /*######################################################################## */
    //GRUN: 1e-18g to 1g dust
    std::vector<double> mass_vec;
    std::vector<double> mass_prob;
    double max_mass_grun = 1; //grams
    double D_mass;
    double max_mass_used;

    out << "# Using mass distribution type : " << opt.mass_dist; out.endl();

    switch(opt.mass_dist) {
        case 0:
            grun_model_mass_flux = grun_flux(opt.mass_min,opt.mass_dist_res);
            //grun model particles/m^2/year, but year and area does not matter, distribution is important
            // so probabillity is vector devided by number of particles
            
            grun_model_mass_flux = grun_model_mass_flux/(sum_v(grun_model_mass_flux));

            D_mass = (max_mass_grun - opt.mass_min)/opt.mass_dist_res;
            mass_vec.clear();
            for(ui=0; ui < opt.mass_dist_res; ui++) {
                mass_vec.push_back(opt.mass_min + ((double)ui + 0.5)*D_mass);
            } 

            temp_mat2.clear();
            mass_vec = mass_vec*1e-3; //grams to kg
            temp_mat2.push_back(mass_vec); 
            temp_mat2.push_back(grun_model_mass_flux);
            max_mass_used = max_mass_grun;
            save_mat(PBE_m_dist_in, &temp_mat2);
            break;
        case 1:
            out << "# Mass min : " << opt.mass_min << " g"; out.endl();
            out << "# Mass max : " << opt.mass_max << " g"; out.endl();
            if(opt.mass_max == opt.mass_min) {
                mass_vec.clear();
                mass_prob.clear();
                temp_mat2.clear();
                mass_vec.push_back(opt.mass_min*1e-3);  //grams to kg
                mass_vec.push_back(opt.mass_max*1e-3);
                mass_prob.push_back(0.5);
                mass_prob.push_back(0.5);
                temp_mat2.push_back(mass_vec);
                temp_mat2.push_back(mass_prob);
            }
            else {
                D_mass = (opt.mass_max - opt.mass_min)/opt.mass_dist_res;
                mass_vec.clear();
                mass_prob.clear();
                for(ui=0; ui < opt.mass_dist_res; ui++) {
                    mass_vec.push_back(opt.mass_min + ((double)ui + 0.5)*D_mass);
                    mass_prob.push_back(1.0/opt.mass_dist_res);
                } 

                temp_mat2.clear();
                mass_vec = mass_vec*1e-3;  //grams to kg
                temp_mat2.push_back(mass_vec);
                temp_mat2.push_back(mass_prob);
            }
            max_mass_used = opt.mass_max;
            save_mat(PBE_m_dist_in, &temp_mat2);
            break;
        case 2:
            out << "# Mass min : " << opt.mass_min << " g"; out.endl();
            out << "# Mass max : " << opt.mass_max << " g"; out.endl();
            D_mass = (log10(opt.mass_max) - log10(opt.mass_min))/opt.mass_dist_res;
            mass_vec.clear();
            mass_prob.clear();
            for(ui=0; ui < opt.mass_dist_res; ui++) {
                mass_vec.push_back(pow(10.0,log10(opt.mass_min) + ((double)ui + 0.5)*D_mass));
                mass_prob.push_back(1.0/opt.mass_dist_res);
            } 

            temp_mat2.clear();
            mass_vec = mass_vec*1e-3;  //grams to kg
            temp_mat2.push_back(mass_vec);
            temp_mat2.push_back(mass_prob);
            max_mass_used = opt.mass_max;
            save_mat(PBE_m_dist_in, &temp_mat2);
            break;
        case 3:
            out << "# Mass data loaded from: " << PBE_mass_config; out.endl();
            temp_mat2.clear();
            error(load_file(&temp_mat2,PBE_mass_config,DATA_MATRIX),matrix_created_files);
            temp_mat2[0] = temp_mat2[0]*1e-3; //grams to kg
            temp_mat2[1] = temp_mat2[1]/(sum_v(temp_mat2[1]));

            max_mass_used = 0;
            for(ui=0; ui < temp_mat2[0].size(); ui++) {
                if(temp_mat2[0][ui] > max_mass_used) {
                    max_mass_used = temp_mat2[0][ui];
                }
            }

            save_mat(PBE_m_dist_in, &temp_mat2);
            break;
        default:
            error(NOT_SUPPORTED_MASS_DIST,matrix_created_files);
    }

    /*######################################################################## */
    // CALCULATE ORBITAL PDF's
    /*######################################################################## */

    statistical_set semi_axis;
    statistical_set ecc;
    statistical_set inc;
    statistical_set omega;
    statistical_set Omega;
    statistical_set tp;

    std::vector<std::vector<double> > DIST_DATA;
    std::vector<std::vector<double> > DIST_DATA_BINS;
    std::vector<double> DIST_DATA_HIST;
    std::vector<double> low_edges;
    std::vector<double> high_edges;
    std::vector<double> DIST_DATA_DX;
    std::vector<double> temp_DIST_data;
    std::vector<double> temp_orb_data; 
    double DIST_DATA_temp_val;
    unsigned long long int DIST_DATA_DIM, bin_N;
    double side_N;
    std::vector<std::vector<double> > RND_DATA;


    if(opt.family >= 0 && opt.family <= 9) {
        orbital_dist_data = swap_cols(orbital_dist_data,3,4);
    }

    out << "Orbital sample of " << orbital_dist_data.size() << " orbits loaded."; out.endl();
    switch(opt.dist_type) {
        case 0://----------------------------------------------------------------------------------
            out << "Creating sample data structure"; out.endl();
            for(i=0; i < orbital_dist_data.size(); i++) {
                temp_vec.clear();
                if(opt.family >= 0 && opt.family <= 10) {
                    temp_vec.push_back(orbital_dist_data[i][0]/(1.0 - orbital_dist_data[i][1]));
                }
                else {
                    temp_vec.push_back(orbital_dist_data[i][0]);
                }
                temp_vec.push_back(orbital_dist_data[i][1]);
                temp_vec.push_back(orbital_dist_data[i][2]);
                temp_vec.push_back(orbital_dist_data[i][3]);
                temp_vec.push_back(orbital_dist_data[i][4]);
                if(opt.family == 11) {
                    temp_vec.push_back(orbital_dist_data[i][5]);
                }
                DIST_DATA.push_back(temp_vec);
            }

            PB_kep_dist_data = DIST_DATA;
            break;
        case 1://----------------------------------------------------------------------------------
    
    if(opt.family >= 0 && opt.family <= 10) {
        semi_axis.data = (extract_col(orbital_dist_data,0)/((double)1 - extract_col(orbital_dist_data,1)));
    }
    else {
        semi_axis.data = (extract_col(orbital_dist_data,0));
    }

    ecc.data = (extract_col(orbital_dist_data,1));
    inc.data = (extract_col(orbital_dist_data,2));
    omega.data = (extract_col(orbital_dist_data,3));
    Omega.data = (extract_col(orbital_dist_data,4));
    
    if(opt.family == 11) {
        tp.data = extract_col(orbital_dist_data,5);
        out.endl();out << "Creating histogram for perihelion passage"; out.endl();
        tp.convert_to_PDF_and_clear(0);
    }
     out.endl();out << "Creating histogram for semi major axis"; out.endl();
    semi_axis.convert_to_PDF_and_clear(0);
     out.endl();out << "Creating histogram for eccentricity"; out.endl();
    ecc.convert_to_PDF_and_clear(0);
     out.endl();out << "Creating histogram for inlcination"; out.endl();
    inc.convert_to_PDF_and_clear(0);
     out.endl();out << "Creating histogram for argument of perihelion"; out.endl();
    omega.convert_to_PDF_and_clear(0);
     out.endl();out << "Creating histogram for longitude of ascending node"; out.endl();
    Omega.convert_to_PDF_and_clear(0);



    out << "Generating distribution file"; out.endl();
    PB_kep_dist_data.clear();
    temp_mat.clear();
    temp_mat = semi_axis.save_histogram_to_file();
    PB_kep_dist_data.push_back(temp_mat[0]);
    PB_kep_dist_data.push_back(temp_mat[1]);

        temp_mat.clear();
    temp_mat = ecc.save_histogram_to_file();
    PB_kep_dist_data.push_back(temp_mat[0]);
    PB_kep_dist_data.push_back(temp_mat[1]);

        temp_mat.clear();
    temp_mat = inc.save_histogram_to_file();
    PB_kep_dist_data.push_back(temp_mat[0]);
    PB_kep_dist_data.push_back(temp_mat[1]);

        temp_mat.clear();
    temp_mat = omega.save_histogram_to_file();
    PB_kep_dist_data.push_back(temp_mat[0]);
    PB_kep_dist_data.push_back(temp_mat[1]);

        temp_mat.clear();
    temp_mat = Omega.save_histogram_to_file();
    PB_kep_dist_data.push_back(temp_mat[0]);
    PB_kep_dist_data.push_back(temp_mat[1]);

    if(opt.family == 11) {
        temp_mat.clear();
        temp_mat = tp.save_histogram_to_file();
        PB_kep_dist_data.push_back(temp_mat[0]);
        PB_kep_dist_data.push_back(temp_mat[1]);
    }
            break;
        case 2: //----------------------------------------------------------------------------------
    out << "Creating sample data structure"; out.endl();
    for(i=0; i < orbital_dist_data.size(); i++) {
        temp_vec.clear();
        if(opt.family >= 0 && opt.family <= 10) {
            temp_vec.push_back(orbital_dist_data[i][0]/(1.0 - orbital_dist_data[i][1]));
        }
        else {
            temp_vec.push_back(orbital_dist_data[i][0]);
        }
        temp_vec.push_back(orbital_dist_data[i][1]);
        temp_vec.push_back(orbital_dist_data[i][2]);
        temp_vec.push_back(orbital_dist_data[i][3]);
        temp_vec.push_back(orbital_dist_data[i][4]);
        DIST_DATA.push_back(temp_vec);
    }
    if(opt.family == 11) {
        tp.data = extract_col(orbital_dist_data,5);
        out.endl();out << "Creating histogram for perihelion passage"; out.endl();
        tp.convert_to_PDF_and_clear(0);
    }
    DIST_DATA_DIM = 5;

    out << "Generating 5 dimensional hypercube bins"; out.endl();
    DIST_DATA_BINS = generate_N_bin_list(DIST_DATA, 0);
    bin_N = DIST_DATA_BINS.size();
    side_N = round(pow((double)bin_N,1.0/(double)DIST_DATA_DIM));
    out << "Hypercube split into " << bin_N << " bins where each side is " << side_N << " segments"; out.endl();

    out << "Finding corners and bin volumes"; out.endl();
    for(i=0; i < DIST_DATA_DIM; i++) {
        temp_DIST_data = extract_col(DIST_DATA, i);
        DIST_DATA_temp_val = *std::min_element(temp_DIST_data.begin(),temp_DIST_data.end());
        low_edges.push_back(DIST_DATA_temp_val);
        
        DIST_DATA_temp_val = *std::max_element(temp_DIST_data.begin(),temp_DIST_data.end());
        high_edges.push_back(DIST_DATA_temp_val);
    }
    temp_DIST_data.clear();

    for(i=0; i < DIST_DATA_DIM; i++) {
        DIST_DATA_DX.push_back((high_edges[i] - low_edges[i])/side_N);
    }
    out << "Generating multi dimensional histogram"; out.endl();
    DIST_DATA_HIST = generate_N_hist_list(DIST_DATA, DIST_DATA_BINS, DIST_DATA_DX);

    out << "Generating distribtuion file"; out.endl();
    PB_kep_dist_data.clear();
    PB_kep_dist_data = DIST_DATA_BINS;
    for(i=0; i < PB_kep_dist_data.size(); i++) {
        PB_kep_dist_data[i].insert(PB_kep_dist_data[i].begin(),DIST_DATA_HIST[i]);
    }

            break;
    }
    // ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤

    save_mat(PB_kepler_dist, &PB_kep_dist_data);
     out.endl();out << "- Saved file: " << PB_kepler_dist; out.endl();
    temp_mat.clear();
    PB_kep_dist_data.clear();
    
    // ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤

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

    switch(opt.IC_gen) {
        case 0:
            strs << opt.start_date; strs_str = strs.str();

            execute_command = JPL_program;
            execute_command.insert(0,selfpath);
            execute_command = execute_command + " " + selfpath + JPL_folder + "bin1550to2550.430 " + strs_str + " 10 PosVel 0 " + JPL_OUT_folder;
            out << "Calling " << execute_command; out.endl();
        
            out.endl(); out << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" ; out.endl();
            MODULE_RETURN = exec(execute_command.c_str());
            out << MODULE_RETURN;
            out << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" ; out.endl();out.endl(); 

            num = "0";
            for(i=0; i < 8; i++) {
                execute_command = JPL_program;
                execute_command.insert(0,selfpath);
                execute_command = execute_command  + " " + selfpath + JPL_folder + "bin1550to2550.430 " + strs_str + " " + num + " PosVel 0 " + JPL_OUT_folder;
                out << "Calling " << execute_command; out.endl();
                out.endl(); out << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" ; out.endl();
                MODULE_RETURN = exec(execute_command.c_str());
                out << MODULE_RETURN;
                out << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" ; out.endl();out.endl(); 
                num[0]++;
            }
            break;
        case 1:

            time_j2000_diff = (opt.start_date-2451545.0)*(3600.0*24.0); //sec

            strs << time_j2000_diff; strs_str = strs.str();

            for(i=0; i < BODIES_NAMES.size(); i++) {
                execute_command = selfpath + SPICE_folder + SPICE_program + " " + BODIES_NAMES[i] + " " + strs_str + " " + JPL_OUT_folder + " " + "NONE" + " " + FILES_LIST;
                out << "Calling " << execute_command; out.endl();
                out.endl(); out << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" ; out.endl();
                MODULE_RETURN = exec(execute_command.c_str());
                out << MODULE_RETURN;
                out << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" ; out.endl();out.endl(); 
            }
            break;
    }
    


    strs.clear();
    std::vector<double> current_body_kep;
    current_body_kep.push_back(0);
    current_body_kep.push_back(0);
    current_body_kep.push_back(0);
    current_body_kep.push_back(0);
    current_body_kep.push_back(0);
    current_body_kep.push_back(0);

    std::vector<double> current_body_char;
    current_body_char.push_back(0);
    current_body_char.push_back(0);
    current_body_char.push_back(0);
    current_body_char.push_back(0);


    temp_q.clear(); temp_v.clear(); temp_m.clear();
    out << "# Loaded Planet position data"; out.endl();
    error(load_file(&temp_q,JPL_PBE_pos,DATA_MATRIX),matrix_created_files);
    out << "# Loaded Planet velocity data"; out.endl();
    error(load_file(&temp_v,JPL_PBE_vel,DATA_MATRIX),matrix_created_files);
    out << "# Loaded Planet mass data"; out.endl();
    error(load_file(&temp_m,planet_m_init,DATA_VECTOR),matrix_created_files);

    switch(opt.mass_format) {
        case 0:
            out << "# Mass in kg"; out.endl();
            break;
        case 1:
            out << "# Converting mass to kg"; out.endl();
            out << "--------- AU**2/D**3 -> kg ----------"; out.endl();
            for(i=0; i<temp_m.size(); i++) {
                out << "MASS " << i << ": " << temp_m[i] << " -> "; 
                temp_m[i] = temp_m[i]/pow(G_gauss,2)*Msol;
                out << temp_m[i]; out.endl();
            }

            break;
    }
    temp_mat.clear();
    temp_mat.push_back(temp_m);
    out << "# Saving mass data to PBE"; out.endl();
    save_mat(PBE_m_init, &temp_mat);


    switch(opt.IC_gen) {
        case 0:
            out << "# Converting to heliocentric coordinate system"; out.endl();
            heliocentric_qv(&temp_q,&temp_v);
            out << "# Finding and rotating to invariable plane (rot z then rot y)"; out.endl();
            rot_to_invariable_plane(&temp_q,&temp_v,&temp_m);

            out << "# Removing old PBE Planet position data"; out.endl();
            if(file_exists(JPL_PBE_pos)) {
                remove(JPL_PBE_pos.c_str());
                out << "Removed file: " << JPL_PBE_pos; out.endl();
            }
            out << "# Removing old PBE Planet velocity data"; out.endl();
            if(file_exists(JPL_PBE_vel)) {
                remove(JPL_PBE_vel.c_str());
                out << "Removed file: " << JPL_PBE_vel; out.endl();
            }
            out << "# Saving invariable plane relative Planet position data to PBE"; out.endl();
            save_mat(JPL_PBE_pos, &temp_q);
            out << "# Saving invariable plane relative Planet velocity data to PBE"; out.endl();
            save_mat(JPL_PBE_vel, &temp_v);
            break;
        case 1:
            out << "# No coordinate transformation to be done, using J2000.0"; out.endl();
            break;

    }

    /*######################################################################## */
    // LOAD MURMHED DATA 
    /*######################################################################## */

unsigned long long int last_total_PB_created;

    // FROM HERE NO MORE DECLARATIONS -----------------------------------
out << "##############################"; out.endl();
out << "###### Entering MC mode ######"; out.endl();
out << "##############################"; out.endl();

total_showers_created = 0;
total_X_created = 0;
total_n_created = 0;
total_PB_created = 0;
total_met_collided = 0;
propagation_counter = 0;
total_particles_integrated = 0;
total_close_encouters_investigated = 0;

max_particle_number_mercury = opt.mercury6_max;

body_norm_dist.clear();
body_norm_dist_bins.clear();

for(i=0; i < body_family.size(); i++) {
    switch((unsigned int)body_family[i][0]) {
        case 1:
            out << "# Creating normal distributions for variable " << i; out.endl();
            body_norm_dist.push_back(EMPTY_DOUBLE_VEC);
            body_norm_dist_bins.push_back(EMPTY_DOUBLE_VEC);

            body_norm_dist[i] = normal_distribution_density_vector_above_zero(body_family[i][1],body_family[i][2],opt.SUOC_sigma, opt.SUOC_hist_res);
            body_norm_dist_bins[i] = normal_distribution_bin_mid_vector_above_zero(body_family[i][1],body_family[i][2],opt.SUOC_sigma, opt.SUOC_hist_res);
            
            break;
    }
}

bool MC_not_done = TRUE;
bool not_done = TRUE;
do {
out << ". Progress_check_mc_loop_start"; out.endl();
not_done = TRUE;
do {
out << ". Progress_check_PB_gen_start"; out.endl();

    tStart = clock();

    //GENERATE A IC FOR PBE

    for(i=0; i < 4; i++) {
        switch((unsigned int)body_family[i][0]) {
            case 0: 
                current_body_char[i] = draw_from_uniform(body_family[i][1],body_family[i][2]);
                break;
            case 1: 
                current_body_char[i] = body_norm_dist_bins[i][draw_from_dist(body_norm_dist[i])];
                break;
        }
    }
    if(body_family.size() < 5) {
        body_mass = 4.0/3.0*PI*pow(current_body_char[0],3)*current_body_char[1];
    }
    else {
        switch((unsigned int)body_family[4][0]) {
            case 0: 
                body_mass = draw_from_uniform(body_family[4][1],body_family[4][2]);
                break;
            case 1: 
                body_mass = body_norm_dist_bins[4][draw_from_dist(body_norm_dist[4])];
                break;
        }
    }
    

    out << ". Progress_check_PB_ok_orbit_check_start"; out.endl();

    last_total_PB_created = total_PB_created;
    do {
    total_PB_created++;

    switch(opt.dist_type) {
        case 0:
            rand_orb = ((double)rand() / RAND_MAX)*((double)(DIST_DATA.size() - 1));
            rand_orb_ind = (unsigned int)round(rand_orb);
            current_body_kep[0] = DIST_DATA[rand_orb_ind][0];
            current_body_kep[1] = DIST_DATA[rand_orb_ind][1];
            current_body_kep[2] = DIST_DATA[rand_orb_ind][2];
            current_body_kep[3] = DIST_DATA[rand_orb_ind][3];
            current_body_kep[4] = DIST_DATA[rand_orb_ind][4];
            break;
        case 1:
            current_body_kep[0] = draw_from_smoth_dist(semi_axis.PDF,semi_axis.hist_bins);
            current_body_kep[1] = draw_from_smoth_dist(ecc.PDF,ecc.hist_bins);
            current_body_kep[2] = draw_from_smoth_dist(inc.PDF,inc.hist_bins);
            current_body_kep[3] = draw_from_smoth_dist(omega.PDF,omega.hist_bins);
            current_body_kep[4] = draw_from_smoth_dist(Omega.PDF,Omega.hist_bins);
            break;
        case 2:
            temp_orb_data.clear();
            temp_orb_data = draw_from_N_dist(DIST_DATA_HIST,DIST_DATA_BINS,DIST_DATA_DX);
            current_body_kep[0] = temp_orb_data[0];
            current_body_kep[1] = temp_orb_data[1];
            current_body_kep[2] = temp_orb_data[2];
            current_body_kep[3] = temp_orb_data[3];
            current_body_kep[4] = temp_orb_data[4];
            break;
    }


    PB_period = 2*PI*sqrt(pow(current_body_kep[0]*AU,3)/(G*(body_mass+Msol)))/(3600.0*24.0); //days

    if(opt.family == 11) {
        switch(opt.dist_type) {
            case 0:
                TP_diff = DIST_DATA[rand_orb_ind][5] - opt.start_date; //diff from peri
                break;
            case 1:
                TP_diff = tp.hist_bins[draw_from_dist(tp.PDF)] - opt.start_date; //diff from not peri
                break;
            case 2:
                TP_diff = tp.hist_bins[draw_from_dist(tp.PDF)] - opt.start_date; //diff from not peri
                break;
        }

    }
    else {
        TP_diff = 0;
    }

      if(current_body_kep[0]*(1-current_body_kep[1])*1.2 < current_body_char[3]) {
        orbit_ok = TRUE;
      }
      else {
        orbit_ok = FALSE;

        switch((unsigned int)body_family[3][0]) {
            case 0: 
                current_body_char[3] = draw_from_uniform(body_family[3][1],body_family[3][2]);
                break;
            case 1:
                current_body_char[3] = body_norm_dist_bins[3][draw_from_dist(body_norm_dist[3])];
                break;
        }
      }

      temp_mat.clear(); temp_vec.clear();
      temp_vec = current_body_kep;
      temp_vec.insert(temp_vec.begin(),total_PB_created);
      temp_mat.push_back(temp_vec);
      save_mat(statistics_all_PB_kep_out, &temp_mat);
      

    //std::cout << "\r-- Checking PB " << total_PB_created - last_total_PB_created << ", a=" << current_body_kep[0] << ", e=" << current_body_kep[1] << ", i=" << current_body_kep[2] << ", om=" << current_body_kep[3] << ", Om=" << current_body_kep[4] << ", r_c=" << current_body_char[3] << ", 1.2*q=" << current_body_kep[0]*(1-current_body_kep[1])*1.2 << std::flush;
    } while(!orbit_ok);
    //std::cout << "\r-- Checked PB " << total_PB_created - last_total_PB_created << "                      " << std::endl;
    out << ". Progress_check_PB_selected after " << total_PB_created - last_total_PB_created << " loops"; out.endl();
    current_body_char[3] = current_body_char[3]*AU;
    current_body_kep[0] = current_body_kep[0]*AU; //convert from AU to m (for calc consistency)
    current_body_kep[2] = current_body_kep[2]*(PI/180); //radians
    current_body_kep[3] = current_body_kep[3]*(PI/180);
    current_body_kep[4] = current_body_kep[4]*(PI/180);
    current_body_kep[5] = 0; //true anoamly, start at 0=peri, pi=appihelion
    temp_vec.clear();
    temp_vec = kepler_to_xv(current_body_kep, G*(body_mass+Msol)); //sun mass

    current_body_xv.clear();
    current_body_xv = temp_vec;

    out.endl(); out << "-------- IC for PBE configuring to: "; out.endl();
    out << std::setprecision(3) << std::fixed;
    out << "Period                      : " << PB_period/365.25 << " y "; out.endl();
    out << "Appihelion passage          : " << round(opt.start_date + TP_diff - PB_period*0.5) << " JD "; out.endl();
    out << "Perihelion passage          : " << round(opt.start_date + TP_diff) << " JD "; out.endl();
    out << "Sync time                   : " << TP_diff << " JD "; out.endl();
    out << "Perihelion                  : " << (1 - current_body_kep[1])*current_body_kep[0]/AU << " AU "; out.endl();
    out << "Appihelion                  : " << (1 + current_body_kep[1])*current_body_kep[0]/AU << " AU "; out.endl();
    out << "Semi major axis             : " << current_body_kep[0]/AU << " AU "; out.endl();
    out << "Eccentricity                : " << current_body_kep[1]  << " "; out.endl();
    out << "Inclination                 : " << current_body_kep[2]*180.0/PI << " deg"; out.endl();
    out << "Argument of periapsis       : " << current_body_kep[3]*180.0/PI  << " deg "; out.endl();
    out << "Longitude of ascending node : " << current_body_kep[4]*180.0/PI  << " deg "; out.endl();
    out << "True anoamly                : " << current_body_kep[5]*180.0/PI << " deg "; out.endl();
    out << "Comet radius                : " << current_body_char[0]  << " m "; out.endl();
    out << std::scientific;
    out << "Comet mass                  : " << body_mass  << " kg "; out.endl();
    out << std::fixed;
    out << "Bulk density                : " << current_body_char[1] << " kg/m^3"; out.endl();
    out << "Activity factor             : " << current_body_char[2]  << "  "; out.endl();
    out << "Sublimation radius          : " << current_body_char[3]/AU << " AU "; out.endl(); out.endl();
    out << std::setprecision(10) << std::scientific;
    //out << "-------- IC for PBE (q,v): " << sqrt(pow(temp_vec[0],2) + pow(temp_vec[1],2) + pow(temp_vec[2],2))/AU << " " << sqrt(pow(temp_vec[3],2) + pow(temp_vec[4],2) + pow(temp_vec[5],2))*0.001 ; out.endl();
    //out << "q: " << temp_vec[0] << " " << temp_vec[1] << " " << temp_vec[2] ; out.endl();
    //out << "v: " << temp_vec[3]/1000 << " " << temp_vec[4]/1000 << " " << temp_vec[5]/1000; out.endl() ; out.endl();
    

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
    out_stream << current_body_char[0] << " " << current_body_char[1] << " " << current_body_char[2] << " " << current_body_char[3] << " " << body_mass << " " << TP_diff << " " << opt.start_date;
    out_stream << "\r\n";

    out_stream.close();
    out_stream.clear();
    out << "File " << out_file_name << " created"; out.endl();

    EXEC_TIME_ALL[0].push_back((double)(clock() - tStart)/CLOCKS_PER_SEC);
    tStart = clock();

//exit(-1);
    execute_command = PBE_program;
    execute_command.insert(0,selfpath);
    execute_command = execute_command + " " + selfpath + settings_folder + PBE_config;
    out << "Calling " << execute_command; out.endl();

        /*out.endl(); out << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" ; out.endl();
        MODULE_RETURN = exec(execute_command.c_str());
        out << MODULE_RETURN;
        out << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" ; out.endl();out.endl();
*/
        system(execute_command.c_str());
        exit(-1);
    EXEC_TIME_ALL[1].push_back((double)(clock() - tStart)/CLOCKS_PER_SEC);
    tStart = clock();
    


        
        //system(execute_command.c_str());
    //PBE_status = system(execute_command.c_str());
    //out << "PBE exit status: " << PBE_status; out.endl();
    //error(PBE_status,matrix_created_files);

    out << "Loading data files from PBE:"; out.endl();
    out << PBE_q_out; out.endl();
    out << PBE_v_out; out.endl();
    out << PBE_sim_out; out.endl();
    out << PBE_m_out; out.endl();
    temp_q.clear(); temp_v.clear(); temp_m_mat.clear(); temp_t.clear(); temp_sim.clear();
    error(load_data_matrix(&temp_q,PBE_q_out),matrix_created_files);
    error(load_data_matrix(&temp_v,PBE_v_out),matrix_created_files);
    //error(load_data_matrix(&temp_t,PBE_t_out),matrix_created_files);
    error(load_data_matrix(&temp_sim,PBE_sim_out),matrix_created_files);
    error(load_data_matrix(&temp_m_mat,PBE_m_out),matrix_created_files);
    PBE_particle_m.clear();
    PBE_particle_m = temp_m_mat[0];

    out << "Using integrator: " << opt.integrator; out.endl();
    switch(opt.integrator) {
        case 0:
            temp_vec.clear(); out.endl();
            
            temp_q_major.clear();
            temp_v_major.clear();
            temp_m_major.clear();
            out << "#------ Extracting major body data ----- "; out.endl();
            for(i=0; i < 9; i++) {
                temp_q_major.push_back(temp_q[i]/AU);
                temp_v_major.push_back(temp_v[i]);
                temp_m_major.push_back(PBE_particle_m[i]);
            }
            temp_q.erase(temp_q.begin(),temp_q.begin()+9);
            temp_v.erase(temp_v.begin(),temp_v.begin()+9);
            PBE_particle_m.erase(PBE_particle_m.begin(),PBE_particle_m.begin()+9);

            if(file_exists(PBE_m_out)) {
                remove(PBE_m_out.c_str());
                out << "Removed file: " << PBE_m_out; out.endl();
            }
            out << "Saving reduced file: " << PBE_m_out; out.endl();
            temp_mat.clear();
            temp_mat.push_back(PBE_particle_m);
            save_mat(PBE_m_out, &temp_mat);

            if(file_exists(JPL_CMS_pos)) {
                remove(JPL_CMS_pos.c_str());
                out << "Removed file: " << JPL_CMS_pos; out.endl();
            }
            if(file_exists(JPL_CMS_vel)) {
                remove(JPL_CMS_vel.c_str());
                out << "Removed file: " << JPL_CMS_vel; out.endl();
            }
            if(file_exists(CMS_m_init)) {
                remove(CMS_m_init.c_str());
                out << "Removed file: " << CMS_m_init; out.endl();
            }

            out << "Saving major bodies to CMS"; out.endl();
            out << JPL_CMS_pos; out.endl();
            save_mat(JPL_CMS_pos, &temp_q_major);
            out << JPL_CMS_vel; out.endl();
            save_mat(JPL_CMS_vel, &temp_v_major);
            out << CMS_m_init; out.endl();
            temp_mat.clear();
            temp_mat.push_back(temp_m_major);
            save_mat(CMS_m_init, &temp_mat);

            out << "Saving small bodies to CMS"; out.endl();
            out << CMS_prop_in; out.endl();
            temp_mat.clear();
            temp_mat = fill_mat(2,PBE_particle_m.size(),temp_sim[0][2]);
            temp_mat[0].clear(); temp_mat[0] = PBE_particle_m;
            save_mat(CMS_prop_in, &temp_mat);
            out << CMS_q_in; out.endl();
            save_mat(CMS_q_in, &temp_q);
            out << CMS_v_in; out.endl();
            save_mat(CMS_v_in, &temp_v);

            //exit(-1);

            execute_command = CMS_program;
            execute_command.insert(0,selfpath);
            execute_command = execute_command + " " + selfpath + settings_folder + CMS_config;
            out << "Calling " << execute_command; out.endl();
            CMS_status = system(execute_command.c_str());

            out << "CMS exit status: " << CMS_status; out.endl();
            
            break;
        case 1:

            temp_vec.clear(); out.endl();

            out << "# Converting to heliocentric coordinate system"; out.endl();
            heliocentric_qv(&temp_q,&temp_v);

            temp_q_major.clear();
            temp_v_major.clear();
            temp_m_major.clear();
            out << "#------ Extracting major body data ----- "; out.endl();
            for(i=0; i < 9; i++) {
                temp_q_major.push_back(temp_q[i]/AU);
                temp_v_major.push_back(temp_v[i]);
                temp_m_major.push_back(PBE_particle_m[i]);
            }
            temp_q.erase(temp_q.begin(),temp_q.begin()+9);
            temp_v.erase(temp_v.begin(),temp_v.begin()+9);
            PBE_particle_m.erase(PBE_particle_m.begin(),PBE_particle_m.begin()+9);

            if(file_exists(PBE_m_out)) {
                remove(PBE_m_out.c_str());
                out << "Removed file: " << PBE_m_out; out.endl();
            }
            out << "Saving reduced file: " << PBE_m_out; out.endl();
            temp_mat.clear();
            temp_mat.push_back(PBE_particle_m);
            save_mat(PBE_m_out, &temp_mat);

            total_PBE_sim_time = temp_sim[0][1]/(3600.0*24.0);
            PBE_time_sync = temp_sim[0][4]/(3600.0*24.0);
            mercury6_start_time = opt.start_date + total_PBE_sim_time;
            if(opt.end_date == 0) {
                mercury6_stop_time = mercury6_start_time + opt.integration_time*365.25;
            }
            else {
                mercury6_stop_time = opt.end_date;
            }
            

            out << "# Removing old Planet state data"; out.endl();
            if(file_exists(mercury6_big)) {
                remove(mercury6_big.c_str());
                out << "Removed file: " << mercury6_big; out.endl();
            }
            out << "# Print big planet data to mercury6"; out.endl();
            error(print_mercury6_big(mercury6_big,&temp_q_major,&temp_v_major,&temp_m_major,mercury6_start_time),matrix_created_files);
            out << "Saved file: " << mercury6_big; out.endl();
            period_n = (unsigned int)temp_sim[0][5];
            
            time_vec.clear();
            time_vec.push_back(PBE_time_sync);
            time_vec.push_back(opt.start_date + PBE_time_sync);
            time_vec.push_back(mercury6_start_time);
            time_vec.push_back(mercury6_start_time);
            time_vec.push_back(mercury6_stop_time);
            for(i=0; i < period_n; i++) {
                time_vec.push_back(opt.start_date + temp_sim[0][6+i]);
            }

            out.endl();
            out << std::setprecision(2) << std::fixed;
            out << "#------ Time variables ----- "; out.endl();
            out << "# PBE sync time  : " << PBE_time_sync/365.25 << " y"; out.endl();
            out << "# PBE start      : " << JDtoJ(opt.start_date + PBE_time_sync) << " y"; out.endl();
            for(i=0; i < period_n; i++) {
            out << "# Perihelion " << i+1 << "   : " << JDtoJ(opt.start_date + temp_sim[0][6+i]) << " y"; out.endl();
            }
            out << "# PBE stop       : " << JDtoJ(mercury6_start_time) << " y"; out.endl();
            out << "# mercury6 start : " << JDtoJ(mercury6_start_time) << " y"; out.endl();
            out << "# mercury6 stop  : " << JDtoJ(mercury6_stop_time) << " y"; out.endl();
            out << "#--------------------------- "; out.endl();
            out << "# mercury6 start : " << mercury6_start_time << " JD"; out.endl();
            out << "# mercury6 stop  : " << mercury6_stop_time << " JD"; out.endl();
            out << "#--------------------------- "; out.endl();
            out << std::setprecision(10) << std::scientific;
            out.endl();

            current_n_created = temp_q.size()-1;

            out << "# Updating mercury6 start time to " << mercury6_start_time; out.endl();
            out << "# Updating mercury6 stop time to  " << mercury6_stop_time; out.endl();
            if(file_exists(mercury6_param)) {
                remove(mercury6_param.c_str());
                out << "Removed file: " << mercury6_param; out.endl();
            }
            error(print_mercury6_param(mercury6_param, &opt, mercury6_start_time),matrix_created_files);
            out << "Saved file: " << mercury6_param; out.endl();
            out << "# Number of particles loaded     " << current_n_created; out.endl();
            out << "# Number of parent bodies loaded " << 1; out.endl();

            total_n_created += current_n_created;
            used_t.clear();
            for(i=0; i < temp_q.size(); i++) {
                used_t.push_back(mercury6_start_time);
            }

            out << "# Building Small-body initial data: " << mercury6_small; out.endl();
            if(file_exists(mercury6_small)) {
                remove(mercury6_small.c_str());
                out << "Removed file: " << mercury6_small; out.endl();
            }
            error(print_mercury6_small(mercury6_small, &temp_q, &temp_v, &PBE_particle_m, &used_t, temp_sim[0][2], max_particle_number_mercury, &particle_number_mercury),matrix_created_files);
            out << "Saved file: " << mercury6_small; out.endl();
            temp_q.clear(); temp_v.clear(); temp_m_mat.clear();

            
            out << "# Small ini file for mercury created with " << particle_number_mercury; out.endl();
            total_particles_integrated = total_particles_integrated + particle_number_mercury;
            out << "########################################################################"; out.endl();
            out << "# Total_met_collided         = " << total_met_collided; out.endl();
            out << "# total_particles_integrated = " << total_particles_integrated; out.endl();
            out << "# total_close_investigated   = " << total_close_encouters_investigated; out.endl();
            out << "# Total_showers_created      = " << total_showers_created; out.endl();
            out << "# total_X_created            = " << total_X_created << ", " << opt.X - total_X_created << " left"; out.endl();
            out << "# Propagation_counter        = " << propagation_counter; out.endl();
            out << "# Starting new particle propagation"; out.endl();
            out << "########################################################################"; out.endl();
            
            execute_command = "cd " + selfpath + mercury6_folder + "; " + mercury6_program;
            out << "# Calling " << execute_command; out.endl();
            
        out.endl(); out << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" ; out.endl();
        MODULE_RETURN = exec(execute_command.c_str());
        out << MODULE_RETURN;
        out << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" ; out.endl();out.endl(); 
        propagation_counter++;
        
        //exit(-1);

EXEC_TIME_ALL[2].push_back((double)(clock() - tStart)/CLOCKS_PER_SEC);
    tStart = clock();

//"########################################################################";
        //snapshot_date
        if(opt.snapshot_date != 0) {
            out << "# Preparing snapshot"; out.endl();

            snapshot_data_pos.clear();
            snapshot_data_vel.clear();
            snapshot_data_earth.clear();

            if(file_exists(mercury6_element_cfg)) {
                remove(mercury6_element_cfg.c_str());
                out << "Removed file: " << mercury6_element_cfg; out.endl();
            }

            out << "# Building Small-body element extract data file: " << mercury6_element_cfg; out.endl();

            mercury6_element_bodies.clear();
            mercury6_element_bodies.push_back("EARTH");
            for(ui=0; ui < current_n_created; ui++) {
                mercury6_element_bodies.push_back("OBJ" + IntToString(ui));
            }

            error(print_mercury6_element(mercury6_element_cfg,mercury6_element_bodies,1),matrix_created_files);

            execute_command = "cd " + selfpath + mercury6_folder + "; " + element6_program;
            out << "Calling " << execute_command; out.endl();

            out.endl(); out << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" ; out.endl();
            MODULE_RETURN = exec(execute_command.c_str());
            out << MODULE_RETURN;
            out << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" ; out.endl();out.endl(); 

            
            snapshot_offset = (opt.snapshot_date - mercury6_start_time);

            for(ui=0; ui < current_n_created; ui++) {
                temp_obj_data.clear();
                mercury6_temp_obj_str = selfpath + mercury6_folder + "OBJ" + IntToString(ui) + ".aei";
                out << "# Attempting to load: " << mercury6_temp_obj_str; out.endl();
                error(read_mercury6_matrix(&temp_obj_data,mercury6_temp_obj_str),matrix_created_files);


                temp_data = abs_element_v(extract_col(temp_obj_data,0)-snapshot_offset);
                TEMP_ID = std::min_element(temp_data.begin(),temp_data.end());
                TEMP_ID_index = (unsigned int)(TEMP_ID - temp_data.begin());

                temp_vec.clear();
                temp_vec.push_back(propagation_counter);
                temp_vec.push_back(temp_obj_data[TEMP_ID_index][1]);
                temp_vec.push_back(temp_obj_data[TEMP_ID_index][2]);
                temp_vec.push_back(temp_obj_data[TEMP_ID_index][3]);
                snapshot_data_pos.push_back(temp_vec);
                
                temp_vec.clear();
                temp_vec.push_back(propagation_counter);
                temp_vec.push_back(temp_obj_data[TEMP_ID_index][4]);
                temp_vec.push_back(temp_obj_data[TEMP_ID_index][5]);
                temp_vec.push_back(temp_obj_data[TEMP_ID_index][6]);
                snapshot_data_vel.push_back(temp_vec);
            }

            temp_obj_data.clear();
            mercury6_temp_obj_str = selfpath + mercury6_folder + "EARTH" + ".aei";
            out << "# Attempting to load: " << mercury6_temp_obj_str; out.endl();
            error(read_mercury6_matrix(&temp_obj_data,mercury6_temp_obj_str),matrix_created_files);

            temp_data = abs_element_v(extract_col(temp_obj_data,0)-snapshot_offset);
            TEMP_ID = std::min_element(temp_data.begin(),temp_data.end());
            TEMP_ID_index = (unsigned int)(TEMP_ID - temp_data.begin());

            if(TEMP_ID_index > 5) {
                earth_snap_range = 5;
            }
            else {
                earth_snap_range = TEMP_ID_index;
            }

            for(i=-earth_snap_range; i <= earth_snap_range; i++) {
                temp_vec.clear();
                temp_vec.push_back(propagation_counter);
                temp_vec.push_back(temp_obj_data[TEMP_ID_index + i][0]);
                temp_vec.push_back(temp_obj_data[TEMP_ID_index + i][1]);
                temp_vec.push_back(temp_obj_data[TEMP_ID_index + i][2]);
                temp_vec.push_back(temp_obj_data[TEMP_ID_index + i][3]);
                temp_vec.push_back(temp_obj_data[TEMP_ID_index + i][4]);
                temp_vec.push_back(temp_obj_data[TEMP_ID_index + i][5]);
                temp_vec.push_back(temp_obj_data[TEMP_ID_index + i][6]);
                snapshot_data_earth.push_back(temp_vec);
            }


            mercury6_temp_obj_str = selfpath + mercury6_folder + "EARTH" + ".aei";
            if(file_exists(mercury6_temp_obj_str)) {
                remove(mercury6_temp_obj_str.c_str());
                out << "# Removed file: " << mercury6_temp_obj_str; out.endl();
            }

            for(ui=0; ui < current_n_created; ui++) {
                mercury6_temp_obj_str = selfpath + mercury6_folder + "OBJ" + IntToString(ui) + ".aei";
                if(file_exists(mercury6_temp_obj_str)) {
                    remove(mercury6_temp_obj_str.c_str());
                    out << "# Removed file: " << mercury6_temp_obj_str; out.endl();
                }
            }

            save_mat(particle_snapshot_pos, &snapshot_data_pos);
            out << "# Saved data to: " << particle_snapshot_pos; out.endl();
            save_mat(particle_snapshot_vel, &snapshot_data_vel);
            out << "# Saved data to: " << particle_snapshot_vel; out.endl();
            save_mat(particle_snapshot_earth, &snapshot_data_earth);
            out << "# Saved data to: " << particle_snapshot_earth; out.endl();
            
            snapshot_data_earth.clear();
            snapshot_data_pos.clear();
            snapshot_data_vel.clear();

        }
//"########################################################################";

EXEC_TIME_ALL[3].push_back((double)(clock() - tStart)/CLOCKS_PER_SEC);
    tStart = clock();

            //mercury6_status = system(execute_command.c_str());

            //out << "- mercury6 exit status: " << mercury6_status; out.endl();
            //error(mercury6_status,matrix_created_files);
            //exit(-1);
            execute_command = "cd " + selfpath + mercury6_folder + "; " + close6_program;
            out << "# Calling " << execute_command; out.endl();
            
        out.endl(); out << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" ; out.endl();
        MODULE_RETURN = exec(execute_command.c_str());
        out << MODULE_RETURN;
        out << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" ; out.endl();out.endl(); 

EXEC_TIME_ALL[4].push_back((double)(clock() - tStart)/CLOCKS_PER_SEC);


            //mercury6_status = system(execute_command.c_str());
            //out << "- close6 exit status: " << mercury6_status; out.endl();

            encounter_data.clear();
            if(MODULE_RETURN.compare("ERROR") != 0) {
                error(read_mercury6_matrix(&encounter_data,mercury6_earth_clo),matrix_created_files);
            }
            else {
                out << "-------- close6 ERROR ENCOUNTERED, SKIPPING LOOP"; out.endl();
            }

            //error(mercury6_status,matrix_created_files);
            filtered_encounter_data.clear();
//exit(-1);
            total_close_encouters_investigated = total_close_encouters_investigated + encounter_data.size();
            if(encounter_data.size() > 0) {
                for(i=(encounter_data.size() - 1); i >= 0; i--) {
                    if((int)encounter_data[i][1] == -1) {
                        encounter_data.erase(encounter_data.begin()+i);
                        out << "Parent body made close enoucnter with earth, removing from close data"; out.endl();
                    }
                }

                id_vec.clear();
                orig_id_vec.clear();
                out << "Encounters detected, checking hill radii of: " << round(encounter_data.size()/2) << " objects"; out.endl();
                //Check for encounters inside hill radii
                for(ui=0; ui < encounter_data.size(); ui=ui+2) {
                    /*out << std::setprecision(2) << std::fixed;
                    out << "IP OBJ" << ui << " : " << (Rearth + Hatm)/(encounter_data[ui][2]*AU); out.endl();
                    out << "            : " << encounter_data[ui][2] << " AU"; out.endl();
                    out << std::setprecision(10) << std::scientific;*/
                    

                    if(encounter_data[ui][2] <= (earth_hill_r*opt.close_match)) {
                        id_vec.push_back((unsigned int)encounter_data[ui][1]);
                        orig_id_vec.push_back((unsigned int)encounter_data[ui][1]);
                        filtered_encounter_data.push_back(encounter_data[ui]);
                        filtered_encounter_data.back().insert(filtered_encounter_data.back().end(),encounter_data[ui+1].begin(),encounter_data[ui+1].end());
                    }
                }

                encounter_data.clear();
                out.endl(); out <<"# Encounters filtered down to: " << filtered_encounter_data.size() << " objects"; out.endl();
                
                if(filtered_encounter_data.size() > 0) {

                    temp_mat.clear(); temp_vec.clear();

                    error(load_data_matrix(&temp_mat,PBE_m_out),matrix_created_files);
                    temp_mat[0].insert(temp_mat[0].begin(),total_PB_created);
                    save_mat(PBE_all_m_statistics, &temp_mat);
                    temp_mat.clear(); temp_vec.clear();

                    error(load_data_matrix(&temp_mat,PBE_abs_v_out),matrix_created_files);
                    temp_mat[0].insert(temp_mat[0].begin(),total_PB_created);
                    temp_mat[1].insert(temp_mat[1].begin(),total_PB_created);
                    save_mat(PBE_abs_v_statistics, &temp_mat);
                    temp_mat.clear(); temp_vec.clear();

                    error(load_data_matrix(&temp_mat,PBE_sim_out),matrix_created_files);
                    temp_vec.push_back(temp_mat[0][0]);
                    temp_mat.clear();
                    temp_mat.push_back(temp_vec);
                    save_mat(PBE_weight_statistics, &temp_mat);
                    temp_mat.clear(); temp_vec.clear();
                    
                    error(load_data_matrix(&temp_mat,PBE_M_comet_out),matrix_created_files);
                    temp_mat = matrix_transpose(temp_mat);
                    save_mat(PBE_M_statistics, &temp_mat);
                    temp_mat.clear(); temp_vec.clear();
                    
                    temp_vec = time_vec;
                    temp_vec.insert(temp_vec.begin(),total_PB_created);
                    temp_mat.push_back(temp_vec);
                    save_mat(simulation_time_data, &temp_mat);
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

switch(opt.X_type) {
    case 0:
    total_X_created = total_showers_created;
        break;
    case 1:
    total_X_created = total_PB_created;
        break;
}

if(total_X_created >= opt.X) {
    not_done = FALSE;
}

total_met_collided+= filtered_encounter_data.size();
out << "## Total meteoroids encountered incresed by " << filtered_encounter_data.size() << " to " << total_met_collided << "."; out.endl();

    tStart = clock();

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


if(total_X_created < opt.X) {

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

out << "# Building Small-body element extract data file: " << mercury6_element_cfg; out.endl();

mercury6_element_bodies.clear();
mercury6_element_bodies.push_back("EARTH");
mercury6_element_bodies.push_back("PB");
if(opt.stream_calc) {
    for(ui=0; ui < current_n_created; ui++) {
        mercury6_element_bodies.push_back("OBJ" + IntToString(ui));
    }
}
else {
    for(ui=0; ui < id_vec.size(); ui++) {
        mercury6_element_bodies.push_back("OBJ" + IntToString(id_vec[ui]));
    }
}


error(print_mercury6_element(mercury6_element_cfg,mercury6_element_bodies,0),matrix_created_files);

out << "Created file: " << mercury6_element_cfg; out.endl();

execute_command = "cd " + selfpath + mercury6_folder + "; " + element6_program;
out << "Calling " << execute_command; out.endl();

out.endl(); out << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" ; out.endl();
        MODULE_RETURN = exec(execute_command.c_str());
        out << MODULE_RETURN;
        out << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" ; out.endl();out.endl(); 

//mercury6_status = system(execute_command.c_str());
//out << "element6 exit status: " << mercury6_status; out.endl();
//error(mercury6_status,matrix_created_files);

EXEC_TIME_ALL[5].push_back((double)(clock() - tStart)/CLOCKS_PER_SEC);
    tStart = clock();

encounter_data.clear();
earth_encounter_data.clear();
filtered_encounter_data_save.clear();
PB_time_series_data.clear();

D_SH_diff_mat_time_series.clear();
D_D_diff_mat_time_series.clear();
D_rho2_diff_mat_time_series.clear();
D_varrho1_diff_mat_time_series.clear();

D_SH_stream_dissipation_time_series.clear();
D_D_stream_dissipation_time_series.clear();
D_rho2_stream_dissipation_time_series.clear();
D_varrho1_stream_dissipation_time_series.clear();

    out << "Attempting to load: " << mercury6_PB_aei; out.endl();
    error(read_mercury6_matrix(&PB_time_series_data,mercury6_PB_aei),matrix_created_files);

D_SH_diff_mat_time_series.push_back(EMPTY_DOUBLE_VEC);
D_D_diff_mat_time_series.push_back(EMPTY_DOUBLE_VEC);
D_rho2_diff_mat_time_series.push_back(EMPTY_DOUBLE_VEC);
D_varrho1_diff_mat_time_series.push_back(EMPTY_DOUBLE_VEC);
for(uj=0; uj < 3; uj++) {
    D_SH_stream_dissipation_time_series.push_back(EMPTY_DOUBLE_VEC);
    D_D_stream_dissipation_time_series.push_back(EMPTY_DOUBLE_VEC);
    D_rho2_stream_dissipation_time_series.push_back(EMPTY_DOUBLE_VEC);
    D_varrho1_stream_dissipation_time_series.push_back(EMPTY_DOUBLE_VEC);
}

internal_dissipation_vector.clear();
for(uj=0; uj < PB_time_series_data.size(); uj++) {
    internal_dissipation_vector.push_back(0);

    D_SH_diff_mat_time_series[0].push_back(PB_time_series_data[uj][0]);
    D_D_diff_mat_time_series[0].push_back(PB_time_series_data[uj][0]);
    D_rho2_diff_mat_time_series[0].push_back(PB_time_series_data[uj][0]);
    D_varrho1_diff_mat_time_series[0].push_back(PB_time_series_data[uj][0]);

    D_SH_stream_dissipation_time_series[0].push_back(PB_time_series_data[uj][0]);
    D_D_stream_dissipation_time_series[0].push_back(PB_time_series_data[uj][0]);
    D_rho2_stream_dissipation_time_series[0].push_back(PB_time_series_data[uj][0]);
    D_varrho1_stream_dissipation_time_series[0].push_back(PB_time_series_data[uj][0]);

    D_SH_stream_dissipation_time_series[1].push_back(0);
    D_D_stream_dissipation_time_series[1].push_back(0);
    D_rho2_stream_dissipation_time_series[1].push_back(0);
    D_varrho1_stream_dissipation_time_series[1].push_back(0);

    D_SH_stream_dissipation_time_series[2].push_back(0);
    D_D_stream_dissipation_time_series[2].push_back(0);
    D_rho2_stream_dissipation_time_series[2].push_back(0);
    D_varrho1_stream_dissipation_time_series[2].push_back(0);
}

if(opt.stream_calc) {
    out << "# Starting stream dissipation calculation"; out.endl();

    for(ui=0; ui < current_n_created-1; ui++) {
        temp_obj_data.clear();
        mercury6_temp_obj_str = selfpath + mercury6_folder + "OBJ" + IntToString(ui) + ".aei";
        //out << "# Attempting to load: " << mercury6_temp_obj_str; out.endl();
        error(read_mercury6_matrix(&temp_obj_data,mercury6_temp_obj_str),matrix_created_files);
        
        out << "# Calculating row " << ui << " of " << current_n_created-1; out.endl();
        for(uj=ui+1; uj < current_n_created; uj++) {
            temp_obj_data2.clear();
            mercury6_temp_obj_str = selfpath + mercury6_folder + "OBJ" + IntToString(uj) + ".aei";
            //out << "# Attempting to load: " << mercury6_temp_obj_str; out.endl();
            error(read_mercury6_matrix(&temp_obj_data2,mercury6_temp_obj_str),matrix_created_files);
            internal_dissipation_counter++;
            
            for(uk=0; uk < PB_time_series_data.size(); uk++) {

                if((temp_obj_data[uk][2] < 1.0) && (temp_obj_data2[uk][2] < 1.0)) {
                    internal_dissipation_vector[uk] = internal_dissipation_vector[uk] + 1;
                    internal_dissipation_counter = internal_dissipation_vector[uk];

                    TP_temp.a = temp_obj_data[uk][1];
                    TP_temp.e = temp_obj_data[uk][2];
                    TP_temp.i = temp_obj_data[uk][3];
                    TP_temp.omega = temp_obj_data[uk][4];
                    TP_temp.Omega = temp_obj_data[uk][5];

                    TP_temp2.a = temp_obj_data2[uk][1];
                    TP_temp2.e = temp_obj_data2[uk][2];
                    TP_temp2.i = temp_obj_data2[uk][3];
                    TP_temp2.omega = temp_obj_data2[uk][4];
                    TP_temp2.Omega = temp_obj_data2[uk][5];

                    if(internal_dissipation_counter == 1) {
                        XNp1 = d_SH( &(TP_temp), &(TP_temp2));
                        D_SH_stream_dissipation_time_series[1][uk] = XNp1;
                        D_SH_stream_dissipation_time_series[2][uk] = 0.0;
                        
                        XNp1 = d_D( &(TP_temp), &(TP_temp2));
                        D_D_stream_dissipation_time_series[1][uk] = XNp1;
                        D_D_stream_dissipation_time_series[2][uk] = 0.0;

                        XNp1 = rho2( &(TP_temp), &(TP_temp2));
                        D_rho2_stream_dissipation_time_series[1][uk] = XNp1;
                        D_rho2_stream_dissipation_time_series[2][uk] = 0.0;

                        XNp1 = varrho1( &(TP_temp), &(TP_temp2));
                        D_varrho1_stream_dissipation_time_series[1][uk] = XNp1;
                        D_varrho1_stream_dissipation_time_series[2][uk] = 0.0;
                    }
                    else {
                        XNp1 = d_SH( &(TP_temp), &(TP_temp2));
                        mean_stream_dissipation = D_SH_stream_dissipation_time_series[1][uk];
                        D_SH_stream_dissipation_time_series[1][uk] = mean_stream_dissipation*(1.0 - 1.0/((double)internal_dissipation_counter)) + XNp1/((double)internal_dissipation_counter);
                        D_SH_stream_dissipation_time_series[2][uk] = D_SH_stream_dissipation_time_series[2][uk]*(1.0 - 1.0/((double)internal_dissipation_counter)) + (XNp1*XNp1 + ((double)internal_dissipation_counter - 1.0)*mean_stream_dissipation*mean_stream_dissipation)/((double)internal_dissipation_counter) - D_SH_stream_dissipation_time_series[1][uk]*D_SH_stream_dissipation_time_series[1][uk];
                        
                        XNp1 = d_D( &(TP_temp), &(TP_temp2));
                        mean_stream_dissipation = D_D_stream_dissipation_time_series[1][uk];
                        D_D_stream_dissipation_time_series[1][uk] = mean_stream_dissipation*(1.0 - 1.0/((double)internal_dissipation_counter)) + XNp1/((double)internal_dissipation_counter);
                        D_D_stream_dissipation_time_series[2][uk] = D_D_stream_dissipation_time_series[2][uk]*(1.0 - 1.0/((double)internal_dissipation_counter)) + (XNp1*XNp1 + ((double)internal_dissipation_counter - 1.0)*mean_stream_dissipation*mean_stream_dissipation)/((double)internal_dissipation_counter) - D_D_stream_dissipation_time_series[1][uk]*D_D_stream_dissipation_time_series[1][uk];

                        XNp1 = rho2( &(TP_temp), &(TP_temp2));
                        mean_stream_dissipation = D_rho2_stream_dissipation_time_series[1][uk];
                        D_rho2_stream_dissipation_time_series[1][uk] = mean_stream_dissipation*(1.0 - 1.0/((double)internal_dissipation_counter)) + XNp1/((double)internal_dissipation_counter);
                        D_rho2_stream_dissipation_time_series[2][uk] = D_rho2_stream_dissipation_time_series[2][uk]*(1.0 - 1.0/((double)internal_dissipation_counter)) + (XNp1*XNp1 + ((double)internal_dissipation_counter - 1.0)*mean_stream_dissipation*mean_stream_dissipation)/((double)internal_dissipation_counter) - D_rho2_stream_dissipation_time_series[1][uk]*D_rho2_stream_dissipation_time_series[1][uk];

                        XNp1 = varrho1( &(TP_temp), &(TP_temp2));
                        mean_stream_dissipation = D_varrho1_stream_dissipation_time_series[1][uk];
                        D_varrho1_stream_dissipation_time_series[1][uk] = mean_stream_dissipation*(1.0 - 1.0/((double)internal_dissipation_counter)) + XNp1/((double)internal_dissipation_counter);
                        D_varrho1_stream_dissipation_time_series[2][uk] = D_varrho1_stream_dissipation_time_series[2][uk]*(1.0 - 1.0/((double)internal_dissipation_counter)) + (XNp1*XNp1 + ((double)internal_dissipation_counter - 1.0)*mean_stream_dissipation*mean_stream_dissipation)/((double)internal_dissipation_counter) - D_varrho1_stream_dissipation_time_series[1][uk]*D_varrho1_stream_dissipation_time_series[1][uk];
                    }

                }
            }
        }
    }

    for(uk=0; uk < PB_time_series_data.size(); uk++) {
        D_SH_stream_dissipation_time_series[2][uk] = sqrt( (((double)internal_dissipation_counter))/(((double)internal_dissipation_counter) - 1.0)*D_SH_stream_dissipation_time_series[2][uk] );
        D_D_stream_dissipation_time_series[2][uk] = sqrt( (((double)internal_dissipation_counter))/(((double)internal_dissipation_counter) - 1.0)*D_D_stream_dissipation_time_series[2][uk] );
        D_rho2_stream_dissipation_time_series[2][uk] = sqrt( (((double)internal_dissipation_counter))/(((double)internal_dissipation_counter) - 1.0)*D_rho2_stream_dissipation_time_series[2][uk] );
        D_varrho1_stream_dissipation_time_series[2][uk] = sqrt( (((double)internal_dissipation_counter))/(((double)internal_dissipation_counter) - 1.0)*D_varrho1_stream_dissipation_time_series[2][uk] );
    }



}

EXEC_TIME_ALL[6].push_back((double)(clock() - tStart)/CLOCKS_PER_SEC);
    tStart = clock();

for(ui=0; ui < id_vec.size(); ui++) {
    temp_obj_data.clear();
    mercury6_temp_obj_str = selfpath + mercury6_folder + "OBJ" + IntToString(id_vec[ui]) + ".aei";
    out << "# Attempting to load: " << mercury6_temp_obj_str; out.endl();
    error(read_mercury6_matrix(&temp_obj_data,mercury6_temp_obj_str),matrix_created_files);

    for(uj=0; uj < filtered_encounter_data.size(); uj++ ) {
        if(orig_id_vec[uj] == id_vec[ui]) {
            first_id = uj;
            break;
        }
    }

    out << "# Starting PB divergance calculation for object " << id_vec[ui]; out.endl();
    D_SH_diff_mat_time_series.push_back(EMPTY_DOUBLE_VEC);
    D_D_diff_mat_time_series.push_back(EMPTY_DOUBLE_VEC);
    D_rho2_diff_mat_time_series.push_back(EMPTY_DOUBLE_VEC);
    D_varrho1_diff_mat_time_series.push_back(EMPTY_DOUBLE_VEC);
    for(uj=0; uj < temp_obj_data.size(); uj++) {

        TP_temp.a = temp_obj_data[uj][1];
        TP_temp.e = temp_obj_data[uj][2];
        TP_temp.i = temp_obj_data[uj][3];
        TP_temp.omega = temp_obj_data[uj][4];
        TP_temp.Omega = temp_obj_data[uj][5];

        PB_temp.a = PB_time_series_data[uj][1];
        PB_temp.e = PB_time_series_data[uj][2];
        PB_temp.i = PB_time_series_data[uj][3];
        PB_temp.omega = PB_time_series_data[uj][4];
        PB_temp.Omega = PB_time_series_data[uj][5];

        D_SH_diff_mat_time_series[ui+1].push_back(d_SH( &(TP_temp), &(PB_temp)));
        D_D_diff_mat_time_series[ui+1].push_back(d_D( &(TP_temp), &(PB_temp)));
        D_rho2_diff_mat_time_series[ui+1].push_back(rho2( &(TP_temp), &(PB_temp)));
        D_varrho1_diff_mat_time_series[ui+1].push_back(varrho1( &(TP_temp), &(PB_temp)));
    }

    temp_data = abs_element_v(extract_col(temp_obj_data,0)-(filtered_encounter_data[first_id][0]));
    TEMP_ID = std::min_element(temp_data.begin(),temp_data.end());
    TEMP_ID_index = (unsigned int)(TEMP_ID - temp_data.begin());

    if((filtered_encounter_data[first_id][0] - temp_obj_data[TEMP_ID_index][0]) < 0 && TEMP_ID_index >= 1) {
        TEMP_ID_index = TEMP_ID_index - 1;
    }

    temp_vec.clear();
    temp_vec.push_back(filtered_encounter_data[first_id][0]); //time
    temp_vec.push_back(filtered_encounter_data[first_id][2]); //Closest distnace
    temp_vec.push_back(filtered_encounter_data[first_id][9]); //a
    temp_vec.push_back(filtered_encounter_data[first_id][10]);//e
    temp_vec.push_back(filtered_encounter_data[first_id][11]);//i
    temp_vec.push_back(filtered_encounter_data[first_id][12]);//longitude of perihelion
    temp_vec.push_back(filtered_encounter_data[first_id][13]);//longitude of ascending node
    temp_vec.push_back(filtered_encounter_data[first_id][14]);//mean anomaly (or mean longitude if e < 1.e-8)
    temp_vec.push_back(PBE_particle_m[id_vec[ui]]);

    filtered_encounter_data_save.push_back(temp_vec);

    out << "# Creating object encounter data structure"; out.endl();
    encounter_data.push_back(temp_obj_data[TEMP_ID_index]);

    temp_vec.clear();
    temp_vec.push_back(filtered_encounter_data[first_id][0]);//time
    temp_vec.push_back(filtered_encounter_data[first_id][2]);//Closest distnace
    temp_vec.push_back(filtered_encounter_data[first_id][3]);//a
    temp_vec.push_back(filtered_encounter_data[first_id][4]);//e
    temp_vec.push_back(filtered_encounter_data[first_id][5]);//i
    temp_vec.push_back(filtered_encounter_data[first_id][6]);//longitude of perihelion
    temp_vec.push_back(filtered_encounter_data[first_id][7]);//longitude of ascending node
    temp_vec.push_back(filtered_encounter_data[first_id][8]);//mean anomaly (or mean longitude if e < 1.e-8)

    out << "# Creating earth encounter data structure"; out.endl();
    earth_encounter_data.push_back(temp_vec);


}

encounter_stat_data.clear();
encounter_stat_data = PB_time_series_data;
for(ui=0; ui < encounter_stat_data.size(); ui++) {
    encounter_stat_data[ui].insert(encounter_stat_data[ui].begin(),total_PB_created);
}
save_mat(statistics_time_series_PB_kep_out, &encounter_stat_data);
out << "# Saved data to: " << statistics_time_series_PB_kep_out; out.endl();

for(ui=0; ui < current_n_created; ui++) {
    mercury6_temp_obj_str = selfpath + mercury6_folder + "OBJ" + IntToString(ui) + ".aei";
    if(file_exists(mercury6_temp_obj_str)) {
        remove(mercury6_temp_obj_str.c_str());
        out << "# Removed file: " << mercury6_temp_obj_str; out.endl();
    }
}

encounter_stat_data.clear();
encounter_stat_data = filtered_encounter_data_save;
for(ui=0; ui < encounter_stat_data.size(); ui++) {
    encounter_stat_data[ui].insert(encounter_stat_data[ui].begin(),total_PB_created);
}
save_mat(met_encounter_data_out, &encounter_stat_data);
out << "# Saved data to: " << met_encounter_data_out; out.endl();

encounter_stat_data.clear();
encounter_stat_data = encounter_data;
for(ui=0; ui < encounter_stat_data.size(); ui++) {
    encounter_stat_data[ui].insert(encounter_stat_data[ui].begin(),total_PB_created);
}
save_mat(met_before_encounter_data_out, &encounter_stat_data);
out << "# Saved data to: " << met_before_encounter_data_out; out.endl();

encounter_stat_data.clear();
encounter_stat_data = earth_encounter_data;
for(ui=0; ui < encounter_stat_data.size(); ui++) {
    encounter_stat_data[ui].insert(encounter_stat_data[ui].begin(),total_PB_created);
}
save_mat(earth_encounter_data_out, &encounter_stat_data);
out << "# Saved data to: " << earth_encounter_data_out; out.endl();

encounter_stat_data.clear();
encounter_stat_data.push_back(current_body_char);
save_mat(statistics_PB_char_out, &encounter_stat_data);
out << "# Saved data to: " << statistics_PB_char_out; out.endl();

encounter_stat_data.clear();
encounter_stat_data.push_back(current_body_kep);
    encounter_stat_data[0][0] = encounter_stat_data[0][0]/AU; //convert from m to AU (for output consistency)
    encounter_stat_data[0][2] = encounter_stat_data[0][2]/(PI/180); //deg
    encounter_stat_data[0][3] = encounter_stat_data[0][3]/(PI/180);
    encounter_stat_data[0][4] = encounter_stat_data[0][4]/(PI/180);
save_mat(statistics_PB_kep_out, &encounter_stat_data);
out << "# Saved data to: " << statistics_PB_kep_out; out.endl();


encounter_stat_data.clear();
encounter_stat_data = D_SH_diff_mat_time_series; // first entry: PB ID, rest = time series
for(ui=0; ui < encounter_stat_data.size(); ui++) {
    encounter_stat_data[ui].insert(encounter_stat_data[ui].begin(),total_PB_created);
}
save_mat(statistics_function_divergance[0], &encounter_stat_data);
out << "# Saved data to: " << statistics_function_divergance[0]; out.endl();

encounter_stat_data.clear();
encounter_stat_data = D_D_diff_mat_time_series; // first entry: PB ID, rest = time series
for(ui=0; ui < encounter_stat_data.size(); ui++) {
    encounter_stat_data[ui].insert(encounter_stat_data[ui].begin(),total_PB_created);
}
save_mat(statistics_function_divergance[1], &encounter_stat_data);
out << "# Saved data to: " << statistics_function_divergance[1]; out.endl();

encounter_stat_data.clear();
encounter_stat_data = D_rho2_diff_mat_time_series; // first entry: PB ID, rest = time series
for(ui=0; ui < encounter_stat_data.size(); ui++) {
    encounter_stat_data[ui].insert(encounter_stat_data[ui].begin(),total_PB_created);
}
save_mat(statistics_function_divergance[2], &encounter_stat_data);
out << "# Saved data to: " << statistics_function_divergance[2]; out.endl();

encounter_stat_data.clear();
encounter_stat_data = D_varrho1_diff_mat_time_series; // first entry: PB ID, rest = time series
for(ui=0; ui < encounter_stat_data.size(); ui++) {
    encounter_stat_data[ui].insert(encounter_stat_data[ui].begin(),total_PB_created);
}
save_mat(statistics_function_divergance[3], &encounter_stat_data);
out << "# Saved data to: " << statistics_function_divergance[3]; out.endl();


if(opt.stream_calc) {

encounter_stat_data.clear();
encounter_stat_data = D_SH_stream_dissipation_time_series; // first entry: PB ID, rest = time series
for(ui=0; ui < encounter_stat_data.size(); ui++) {
    encounter_stat_data[ui].insert(encounter_stat_data[ui].begin(),total_PB_created);
}
save_mat(stream_dissipation[0], &encounter_stat_data);
out << "# Saved data to: " << stream_dissipation[0]; out.endl();

encounter_stat_data.clear();
encounter_stat_data = D_D_stream_dissipation_time_series; // first entry: PB ID, rest = time series
for(ui=0; ui < encounter_stat_data.size(); ui++) {
    encounter_stat_data[ui].insert(encounter_stat_data[ui].begin(),total_PB_created);
}
save_mat(stream_dissipation[1], &encounter_stat_data);
out << "# Saved data to: " << stream_dissipation[1]; out.endl();

encounter_stat_data.clear();
encounter_stat_data = D_rho2_stream_dissipation_time_series; // first entry: PB ID, rest = time series
for(ui=0; ui < encounter_stat_data.size(); ui++) {
    encounter_stat_data[ui].insert(encounter_stat_data[ui].begin(),total_PB_created);
}
save_mat(stream_dissipation[2], &encounter_stat_data);
out << "# Saved data to: " << stream_dissipation[2]; out.endl();

encounter_stat_data.clear();
encounter_stat_data = D_varrho1_stream_dissipation_time_series; // first entry: PB ID, rest = time series
for(ui=0; ui < encounter_stat_data.size(); ui++) {
    encounter_stat_data[ui].insert(encounter_stat_data[ui].begin(),total_PB_created);
}
save_mat(stream_dissipation[3], &encounter_stat_data);
out << "# Saved data to: " << stream_dissipation[3]; out.endl();

    
}

 encounter_stat_data.clear();

EXEC_TIME_ALL[7].push_back((double)(clock() - tStart)/CLOCKS_PER_SEC);
    tStart = clock();

if(opt.calc_clusters) {

out << "########################################################################"; out.endl();
out << "# STARTING ASSOCIATION                                                 #"; out.endl();
out << "########################################################################"; out.endl();

temp_mat.clear();
temp_mat = earth_encounter_data;

if(file_exists(OAA_input_file)) {
    remove(OAA_input_file.c_str());
    out << "Removed file: " << OAA_input_file; out.endl();
}
save_mat(OAA_input_file, &temp_mat);
out << "Saved input file: " << OAA_input_file; out.endl(); 


OAA_config_vector.clear();
OAA_config_vector.push_back(0);         //Analysis type
OAA_config_vector.push_back(1);         //Input format
OAA_config_vector.push_back(1);         //Relative object row
OAA_config_vector.push_back(1);         //Cluster analysis      
OAA_config_vector.push_back(0);         //Cluster analysis type    
OAA_config_vector.push_back(1);         //Parameter sweep        
OAA_config_vector.push_back(0);         //Parameter sweep adaptive step 
OAA_config_vector.push_back(0.2);       //Adaptive step enable limit (fraction)  
OAA_config_vector.push_back(1.7);       //Adaptive step multiplyer   
OAA_config_vector.push_back(0);         //Parameter sweep termination    
OAA_config_vector.push_back(1);         //Error function       
OAA_config_vector.push_back(1e+03);     //Memory allocation (Mb)   
OAA_config_vector.push_back(1);         //Logfile output
OAA_config_vector.push_back(1);         //D_SH criterion        
OAA_config_vector.push_back(0.1);       //D_SH static value    
OAA_config_vector.push_back(0.0001);     //D_SH step   
OAA_config_vector.push_back(1);         //D_D criterion    
OAA_config_vector.push_back(0.2);       //D_D static value   
OAA_config_vector.push_back(0.0001);     //D_D step   
OAA_config_vector.push_back(1);         //rho2 metric     
OAA_config_vector.push_back(0.2);       //rho2 static value   
OAA_config_vector.push_back(0.0001);    //rho2 step    
OAA_config_vector.push_back(1);         //varrho1 metric      
OAA_config_vector.push_back(0.2);       //varrho1 static value    
OAA_config_vector.push_back(0.0001);    //varrho1 step

temp_mat.clear(); 
temp_mat.push_back(OAA_config_vector);

execute_command = selfpath + settings_folder + OAA_config_simple;
if(file_exists(execute_command)) {
    remove(execute_command.c_str());
    out << "Removed file: " << execute_command; out.endl();
}
save_mat(execute_command, &temp_mat);
out << "Saved cfg file: " << execute_command; out.endl();
temp_mat.clear();

execute_command = "cd " + selfpath + OAA_folder + "; " + OAA_program + " " + OAA_arguments + ";";
out << "Calling " << execute_command; out.endl();

        out.endl(); out << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" ; out.endl();
        MODULE_RETURN = exec(execute_command.c_str());
        out << MODULE_RETURN;
        out << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" ; out.endl();out.endl(); 

//OAA_status = system(execute_command.c_str());
//out << "OAA exit status: " << OAA_status; out.endl();
//error(OAA_status,matrix_created_files);

for(ui=0; ui < FUNCTIONS.size(); ui++) {
    temp_mat.clear(); 
    load_data_matrix(&temp_mat,OAA_OUTPUT_profile_files[ui]);
    for(uj=0; uj < temp_mat.size(); uj++) {
        temp_mat[uj].insert(temp_mat[uj].begin(),total_PB_created);
    }   
    save_mat(statistics_association_profile[ui], &temp_mat);

    temp_mat.clear(); 
    load_data_matrix(&temp_mat,OAA_OUTPUT_error_files[ui]);
    for(uj=0; uj < temp_mat.size(); uj++) {
        temp_mat[uj].insert(temp_mat[uj].begin(),total_PB_created);
    }   
    save_mat(statistics_error_profile[ui], &temp_mat);
}
temp_mat.clear(); 

out << "########################################################################"; out.endl();
}
clean_up(matrix_created_files,1);
clean_up(matrix_created_files,2);
clean_up(matrix_created_files,3);

EXEC_TIME_ALL[8].push_back((double)(clock() - tStart)/CLOCKS_PER_SEC);
    tStart = clock();

//        std::cout << std::endl << "Test pause: ";
//    std::cin >> TEST_INT;

}

switch(opt.X_type) {
    case 0:
    total_X_created = total_showers_created;
        break;
    case 1:
    total_X_created = total_PB_created;
        break;
}

if(total_X_created >= opt.X) {
    MC_not_done = FALSE;

    out << "########################################################################"; out.endl();
    out << "# Goal met of total_X_created = " << total_X_created << ", exiting loop"; out.endl();
    out << "########################################################################"; out.endl();
}
else {
    out << "# Goal NOT met, continuing loop"; out.endl();
}

} while(MC_not_done);


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

    out << "# Moving file " << PBE_all_m_statistics; out.endl();
    execute_command = data_out_folder_string + "particle_ejection_m_stat.txt";
    rename(PBE_all_m_statistics.c_str(), execute_command.c_str());

    out << "# Moving file " << PBE_M_statistics; out.endl();
    execute_command = data_out_folder_string + "PBE_M_stat.txt";
    rename(PBE_M_statistics.c_str(), execute_command.c_str());

    out << "# Moving file " << statistics_all_PB_kep_out; out.endl();
    execute_command = data_out_folder_string +  "PB_kep_data.txt";
    rename(statistics_all_PB_kep_out.c_str(), execute_command.c_str());

    out << "# Moving file " << statistics_time_series_PB_kep_out; out.endl();
    execute_command = data_out_folder_string +  "PB_kep_timeS_data.txt";
    rename(statistics_time_series_PB_kep_out.c_str(), execute_command.c_str());

    for(ui=0; ui < FUNCTIONS.size(); ui++) {
        if(opt.calc_clusters) {
            out << "# Moving file " << statistics_association_profile[ui]; out.endl();
            execute_command = data_out_folder_string + "association_profile_" + FUNCTIONS[ui] + ".txt";
            rename(statistics_association_profile[ui].c_str(), execute_command.c_str());

            out << "# Moving file " << statistics_error_profile[ui]; out.endl();
            execute_command = data_out_folder_string + "error_profile_" + FUNCTIONS[ui] + ".txt";
            rename(statistics_error_profile[ui].c_str(), execute_command.c_str());
        }
        out << "# Moving file " << statistics_function_divergance[ui]; out.endl();
        execute_command = data_out_folder_string + "function_divergance_" + FUNCTIONS[ui] + ".txt";
        rename(statistics_function_divergance[ui].c_str(), execute_command.c_str());

        if(opt.stream_calc) {
            out << "# Moving file " << stream_dissipation[ui]; out.endl();
            execute_command = data_out_folder_string + "stream_dissipation_data_" + FUNCTIONS[ui] + ".txt";
            rename(stream_dissipation[ui].c_str(), execute_command.c_str());
        }
    }

    out << "# Moving file " << PB_kepler_dist; out.endl();
    execute_command = data_out_folder_string + "PB_kep_dist.txt";
    rename(PB_kepler_dist.c_str(), execute_command.c_str());

    out << "# Moving file " << PBE_m_dist_in; out.endl();
    execute_command = data_out_folder_string + "PBE_m_dist_in.txt";
    rename(PBE_m_dist_in.c_str(), execute_command.c_str());

    out << "# Moving file " << simulation_time_data; out.endl();
    execute_command = data_out_folder_string + "time_data.txt";
    rename(simulation_time_data.c_str(), execute_command.c_str());

    out << "# Moving file " << particle_snapshot_vel; out.endl();
    execute_command = data_out_folder_string + "snapshot_vel.txt";
    rename(particle_snapshot_vel.c_str(), execute_command.c_str());

        out << "# Moving file " << particle_snapshot_pos; out.endl();
    execute_command = data_out_folder_string + "snapshot_pos.txt";
    rename(particle_snapshot_pos.c_str(), execute_command.c_str());

        out << "# Moving file " << particle_snapshot_earth; out.endl();
    execute_command = data_out_folder_string + "snapshot_earth.txt";
    rename(particle_snapshot_earth.c_str(), execute_command.c_str());

    out << "# Moving file " << log_file; out.endl();
    execute_command = data_out_folder_string + "outlog.txt";
    rename(log_file.c_str(), execute_command.c_str());


EXEC_TIME_ALL[9].push_back((double)(clock() - time0)/CLOCKS_PER_SEC);
    
    execute_command = data_out_folder_string + "execute_time.txt";
    out << "# Saving time data " << execute_command; out.endl();
    save_mat(execute_command, &EXEC_TIME_ALL);

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

