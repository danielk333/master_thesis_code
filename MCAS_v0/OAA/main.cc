/*
 * main.cc
 *
 *  Created on: Jun 17, 2015
 *      Author: dankas
 */

//Orbit association analysis

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

int main (int argc, char *argv[]) {

/*######################################################################## */
 // DECLARE VARIABLES
/*######################################################################## */

    /* ITTERATORS */ 
    int i,j,k;

    int k1;

    /* RANDOM SEED */
    srand (time(NULL));

    int TEST_INT;

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

    char THE_ONE = '1';

    //Streams for copy of file
    std::ifstream  src; 
    std::ofstream  dst;

    //Logfile out
    std::ofstream log_out;

    //DATABASE
    std::vector<OBSERVATION_record> objs;
    std::vector<OBSERVATION_record> real_objs;
    std::vector<OBSERVATION_record> real_objs_batch;
    std::vector<OBSERVATION_record> real_objs_batch2;
    OBSERVATION_record OBS_temp;

    OBSERVATION_record PB_temp;
    OBSERVATION_record TP_temp;

    //FUNCTIONS TO USE
    std::vector<std::string> FUNCTIONS;
    FUNCTIONS.push_back("SH");
    FUNCTIONS.push_back("D");
    FUNCTIONS.push_back("rho2");
    FUNCTIONS.push_back("varrho1");

    /* Generate programs self_path for for dynamics folders */
std::string selfpath = get_selfpath();
selfpath.erase(selfpath.end()-3,selfpath.end());

    std::string INPUT_db_file;
    if(argc >= 3) {
        INPUT_db_file = argv[2];
    }
    else {
        INPUT_db_file = selfpath + "INPUT_DB.txt";
    }

    std::string OAA_config;
    if(argc >= 2) {
        OAA_config = argv[1];
    }
    else {
        OAA_config = selfpath + "OAA_config.cfg";
    }

    std::vector<std::string> OUTPUT_D_mat_folders;
    std::vector<std::string> OUTPUT_D_mat_files;
    std::vector<std::string> OUTPUT_profile_files;
    std::vector<std::string> OUTPUT_error_files;
    std::vector<std::string> OUTPUT_cluster_files;
    for(i=0; i < FUNCTIONS.size(); i++) {
        OUTPUT_D_mat_folders.push_back(selfpath + "OUTPUT_" + FUNCTIONS[i]);
        OUTPUT_D_mat_files.push_back(selfpath + "OUTPUT_" + FUNCTIONS[i] + "_mat.txt");
        OUTPUT_profile_files.push_back(selfpath + "OUTPUT_" + FUNCTIONS[i] + "_association_profile.txt");
        OUTPUT_error_files.push_back(selfpath + "OUTPUT_" + FUNCTIONS[i] + "_error_profile.txt");
        OUTPUT_cluster_files.push_back(selfpath + "OUTPUT_" + FUNCTIONS[i] + "_clusters.txt");
    }

    std::string log_file = selfpath + "outlog.txt";

    std::string TEMP_file = selfpath + "TEMP_FILE.tmp";

    //###########
    //CREATE LIST OF ALL CREATED FILES
    // - FOR OPENING, CLOSING, COPY,
    // - DELETE, CREATE, and CLEANUP
    //###########
    std::vector<std::string> list_created_files;
    std::vector<std::vector<std::string> > matrix_created_files;

    list_created_files.push_back(log_file);
    list_created_files.push_back(TEMP_file);
    matrix_created_files.push_back(list_created_files); //Output restults 0
    list_created_files.clear();

    //###########
    //CLEAN UP:
    // - IN CASE LAST RUN DIDNT EXIT PROPERLY
    // - TO MAKE ROOM FOR NEW RESULTS
    //###########
    clean_up(matrix_created_files,-1);


    /* LOAD OPTIONS */
    std::string config_path = OAA_config;

    if(argc >= 4) {
        if((*argv[3]) == THE_ONE) {
            error(load_file(&opt,config_path,DATA_OPTIONS_VEC),matrix_created_files);
        }
        else {
            error(load_file(&opt,config_path,DATA_OPTIONS),matrix_created_files);
        }
    }
    else {
        error(load_file(&opt,config_path,DATA_OPTIONS),matrix_created_files);
    }


//                std::cout << std::endl << "Test pause: ";
//    std::cin >> TEST_INT;

    // INIT OUTPUT
    log_out.open(log_file.c_str(), std::ios::app);

    output_stream out(std::cout, log_out);
    out.type = opt.logfile;

    out << "-- Output type " << opt.logfile << " selected."; out.endl(); out.endl();

    out << "-- Input file used: " << INPUT_db_file; out.endl(); out.endl();

std::vector<bool> MET_bool_list;
MET_bool_list.push_back(opt.D_SH);
MET_bool_list.push_back(opt.D_D);
MET_bool_list.push_back(opt.rho2);
MET_bool_list.push_back(opt.varrho1);

std::vector<double> MET_step_list;
MET_step_list.push_back(opt.D_SH_step);
MET_step_list.push_back(opt.D_D_step);
MET_step_list.push_back(opt.rho2_step);
MET_step_list.push_back(opt.varrho1_step);

std::vector<double> MET_calc_list;
MET_calc_list.push_back(opt.D_SH_Crit);
MET_calc_list.push_back(opt.D_D_Crit);
MET_calc_list.push_back(opt.rho2_Crit);
MET_calc_list.push_back(opt.varrho1_Crit);

    //###########
    //DATABASE STORAGE
    //###########

    std::vector<std::vector<double> > INPUT_db;
    std::vector<std::vector<double> > INPUT_db_distance_matrix_part;

    //###########
    //SYSTEM VARIABLES
    //###########

    double mem_batch_estimate;
    //double faction_of_avalible_mem;
    unsigned int calculation_batches;
    unsigned int items_a_batch;
    unsigned int residue_items;

    unsigned int residue_toggle_loop;
//bool residue_toggle;

    //###########
    //TEMPORARY VARIABLES
    //###########

  //double r_vec_kep_state[3], v_vec_kep_state[3], mu_kep_state;
  //double p_kep_state, a_kep_state, ecc_kep_state, incl_kep_state, omega_kep_state, argp_kep_state, nu_kep_state, m_kep_state, arglat_kep_state, truelon_kep_state, lonper_kep_state;

    std::vector<std::vector<double> > time_corrected_vectors_obj,time_corrected_vectors_earth;
    std::vector<double> temp_kepler,temp_state_vec,x_temp,v_temp;
    //unsigned int propagation_counter;
    //unsigned int TEMP_ID_index;
    std::vector<double> temp_obj_pos_a;
    std::vector<double> temp_obj_pos_b;
    double d_crit_used;
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

    std::vector<double> current_body_xv;

    double d_crit;
    double d_crit_old;
    double d_crit_step;
    //unsigned int d_crit_bins;
    std::vector<double> association_vec;
    std::vector<double> d_crit_vec;
    double fraction_associated,fraction_associated_old,fraction_associated_next,fraction_associated_min;
    double error_function_current,error_function_old;

    bool TERMINATION_VAR;

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
    std::vector<std::vector<double> > PB_time_series_data;
    std::vector<std::vector<double> > earth_encounter_data;
    std::vector<std::vector<double> > encounter_stat_data;

    std::vector<std::vector<double> > PB_kep_dist_data;

    std::vector<std::vector<double> > D_SH_diff_mat_time_series;
    std::vector<std::vector<double> > D_D_diff_mat_time_series;
    std::vector<std::vector<double> > D_rho2_diff_mat_time_series;
    std::vector<std::vector<double> > D_varrho1_diff_mat_time_series;

std::vector<double> error_function;
double error_min;
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
    bool MET_ABORT;
    unsigned int first_id;
    unsigned int largest_clust_size;
    std::vector<std::vector<double> > body_family;

    //###########
    //Counters
    //###########
    //unsigned int total_n_created, total_PB_created, total_met_collided, current_n_created, total_showers_created;

    //###########
    //MISC
    //###########
    bool orbit_ok;
    //double body_mass;
    std::vector<double>::iterator TEMP_ID;


    /*######################################################################## */
    // LOAD INPUT DATA 
    /*######################################################################## */

out << "###################################"; out.endl();
out << "###### CHECKING INPUT DATA ######"; out.endl();
out << "###################################"; out.endl();

bool DO_THEY_EXIST;
DO_THEY_EXIST = TRUE;

bool DO_WE_CALC_CRITERION;
DO_WE_CALC_CRITERION = (opt.D_SH || opt.D_D || opt.rho2 || opt.varrho1);


for(i=0; i < OUTPUT_D_mat_folders.size(); i++) {
    DO_THEY_EXIST = DO_THEY_EXIST && folder_exists(OUTPUT_D_mat_folders[i]);
}

if(!DO_WE_CALC_CRITERION) {
    out << "# No distance matricies to be calculated"; out.endl();
}    
else {

    if(DO_THEY_EXIST && opt.analysis_type == 1) {
        out << "# All distance matrices pre-calculated"; out.endl();
    }
    else {
        out << "# Distance matrices missing: Filling the gaps"; out.endl();
        temp_file_path = INPUT_db_file;
        error(load_file(&INPUT_db,temp_file_path,DATA_MATRIX),matrix_created_files);

        out << "# Building filtered database distance matrix calculation"; out.endl();
        for(i=0; i < INPUT_db.size(); i++) { 
            MET_ABORT = FALSE;
            //PUT PROGRAM RESTRAINTS HERE (I.E if it cant handle hyperbolic, remove them)
            if(opt.rho2 || opt.varrho1) {
                MET_ABORT = MET_ABORT || (INPUT_db[i][27] >= 1.0);
                // can add more on the same way
                //MET_ABORT = MET_ABORT || (whatever dont work);
            }


            if(!MET_ABORT) {
                // HERE INPUT FORMAT WILL BE USED
                switch(opt.input_format) {
                    case 0: //MURMHED
                        OBS_temp.a = INPUT_db[i][26];
                        OBS_temp.e = INPUT_db[i][27];
                        OBS_temp.i = INPUT_db[i][30];
                        OBS_temp.omega = INPUT_db[i][31];
                        OBS_temp.Omega = INPUT_db[i][29];

                        OBS_temp.ra = INPUT_db[i][9];
                        OBS_temp.dec = INPUT_db[i][10];
                        OBS_temp.v_g = INPUT_db[i][17];
                        OBS_temp.lambda = INPUT_db[i][36];

                        break;
                    case 1: //MCAS meteors
                        OBS_temp.a = INPUT_db[i][2];
                        OBS_temp.e = INPUT_db[i][3];
                        OBS_temp.i = INPUT_db[i][4];
                        OBS_temp.omega = INPUT_db[i][5] - INPUT_db[i][6];
                        OBS_temp.Omega = INPUT_db[i][6];

                        break;
                }

                real_objs.push_back(OBS_temp);
            }
        }
        out << "# " << real_objs.size() << " Objects built, " << (INPUT_db.size() - real_objs.size()) << " shower members skipped"; out.endl();
        INPUT_db.clear();

        if(opt.analysis_type == 0) {
            out << "# BIG DATA dissabled, will not split work even if too much memory used."; out.endl();
            for(j=0; j < FUNCTIONS.size(); j++) {
                if(!(file_exists(OUTPUT_D_mat_files[j])) && MET_bool_list[j]) {
                    out << "# Calculating " << FUNCTIONS[j] << " distance matrix"; out.endl();
                    D_matix.clear();
                    D_matix = calculate_metric_matrix(&real_objs,FUNCTIONS[j]);
                    out << "# " << FUNCTIONS[j] << " distance matrix calculated"; out.endl();
                    save_mat(OUTPUT_D_mat_files[j], &D_matix);
                    out << "# Saved " << FUNCTIONS[j] << " to file: " << OUTPUT_D_mat_files[j]; out.endl();
                }
            }

        }
        else {

            out << "# Avalible RAM: " << opt.mem_allowed << " Mb"; out.endl();
            mem_batch_estimate = ((double)real_objs.size())*((double)real_objs.size() - (double)1)*0.5*((double)sizeof(double));
            mem_batch_estimate = mem_batch_estimate/((double)(8*1024*1024)); //Mb

            calculation_batches = (unsigned int)ceil(mem_batch_estimate/(opt.mem_allowed));
            items_a_batch = (unsigned int)floor(((double)real_objs.size())/((double)calculation_batches));
            residue_items = real_objs.size() - calculation_batches*items_a_batch;

            if(residue_items > 0) {
                residue_toggle_loop = 1;
            }
            else {
                residue_toggle_loop = 0;
            }

            out << "# Estimating memory usage to " << mem_batch_estimate << " Mb, allowed usage is " << opt.mem_allowed << "Mb"; out.endl();
            out << "# Splitting work into " << calculation_batches << " batches, eatch of " << items_a_batch << " items, leaving " << residue_items << " residue items"; out.endl();
            for(j=0; j < FUNCTIONS.size(); j++) {
                if(!(folder_exists(OUTPUT_D_mat_files[j])) && MET_bool_list[j]) {

                    boost::filesystem::create_directory(OUTPUT_D_mat_files[j]);

                    out << "# Calculating " << FUNCTIONS[j] << " distance matrix"; out.endl();
                    for(k=0; k < calculation_batches+residue_toggle_loop; k++) {

                        temp_batch_name = OUTPUT_D_mat_files[j] + "/" + IntToString((unsigned int)(k*items_a_batch)) + "_" + IntToString((unsigned int)((k+1)*items_a_batch)) + "_" + FUNCTIONS[j] + ".txt";

                        for(k1=0; k1 < (calculation_batches+residue_toggle_loop-k); k1++) {
                            real_objs_batch.clear();
                            real_objs_batch2.clear();
                            INPUT_db_distance_matrix_part.clear();
                            std::cout << "- Doing batch coordinate (" << k << "," << k1 << "):" << std::endl;
                            if(k != calculation_batches+1) {
                                real_objs_batch.insert( real_objs_batch.begin(), real_objs.begin() + k*items_a_batch, real_objs.begin() + (k+1)*items_a_batch);
                            }
                            else {
                                real_objs_batch.insert( real_objs_batch.begin(), real_objs.begin() + k*items_a_batch, real_objs.begin() + k*items_a_batch + residue_items);
                            }
                            
                            if(k1 == 0) {
                                INPUT_db_distance_matrix_part = calculate_metric_matrix(&real_objs_batch,FUNCTIONS[j]);
                                error(save_mat(temp_batch_name, &INPUT_db_distance_matrix_part),matrix_created_files);
                            }
                            else if(k1 == (calculation_batches+1-k) && k != calculation_batches+1) {
                                real_objs_batch2.insert( real_objs_batch2.begin(), real_objs.begin() + k1*items_a_batch, real_objs.begin() + k1*items_a_batch + residue_items);
                                INPUT_db_distance_matrix_part = augment_metric_matrix(&real_objs_batch,&real_objs_batch2,FUNCTIONS[j]);
                                error(save_batch(temp_batch_name, TEMP_file, &INPUT_db_distance_matrix_part,k*items_a_batch),matrix_created_files);
                            }
                            else {
                                real_objs_batch2.insert( real_objs_batch2.begin(), real_objs.begin() + k1*items_a_batch, real_objs.begin() + (k1+1)*items_a_batch);
                                INPUT_db_distance_matrix_part = augment_metric_matrix(&real_objs_batch,&real_objs_batch2,FUNCTIONS[j]);
                                error(save_batch(temp_batch_name, TEMP_file, &INPUT_db_distance_matrix_part,k*items_a_batch),matrix_created_files);
                            }
                        }
                    }
                out << "# " << FUNCTIONS[j] << " distance matrix calculated"; out.endl();
                out << "# Saved " << FUNCTIONS[j] << " to file: " << OUTPUT_D_mat_files[j]; out.endl();
                }
            }
        }
        real_objs.clear();

        out << "# Database distance matrix calculation done, clearing data"; out.endl(); out.endl();   
    }
}


if(!(opt.cluster_analysis) || ((opt.analysis_type == 1) && !(opt.parameter_sweep))) {
    out << "# No cluster analysis to be performed"; out.endl();
}
else {
    for(j=0; j < FUNCTIONS.size(); j++) {
        if(!(file_exists(OUTPUT_profile_files[j]))  && MET_bool_list[j]) {
            out << "# Missing database association profiles: Filling the gaps"; out.endl();
            out << "# -- Creating database association profile " << OUTPUT_profile_files[j]; out.endl();

            switch(opt.analysis_type) {
                case 0:
                    D_matix.clear();
                    error(load_file(&D_matix,OUTPUT_D_mat_files[j],DATA_MATRIX),matrix_created_files);
                    out << "- Loaded D matrix file: " << OUTPUT_D_mat_files[j]; out.endl();
                    number_of_objects_db_load = D_matix.size()+1;
                    break;
                case 1:
                    list_of_batches_d_matrix.clear();
                    list_of_batches_d_matrix = list_files(OUTPUT_D_mat_folders[j]);
                    index_matrix.clear();
                    for(k=0; k < list_of_batches_d_matrix.size(); k++) {
                        list_of_batches_d_matrix[k].erase(0,OUTPUT_D_mat_folders[j].size()+1);
                        out << "- Index range recorded for file " << list_of_batches_d_matrix[k]; out.endl();

                        char_index_1 = list_of_batches_d_matrix[k].find("_", 0);
                        char_index_2 = list_of_batches_d_matrix[k].find("_", char_index_1+1);

                        index_matrix.push_back(EMPTY_UINT_VEC);
                        buffer_temp = list_of_batches_d_matrix[k].substr(0,char_index_1);
                        index_matrix[k].push_back((unsigned int)strtod(buffer_temp.c_str(), NULL));

                        buffer_temp = list_of_batches_d_matrix[k].substr(char_index_1+1,char_index_2);
                        index_matrix[k].push_back((unsigned int)strtod(buffer_temp.c_str(), NULL));
                    }
                    out << "- Completed index matrix:"; out.endl();
                    print_matrix(index_matrix,out);
                    number_of_objects_db_load = index_matrix.back().back()-1;
                    break;
            }
            out << "# Objects loaded from D file: " << number_of_objects_db_load; out.endl();  
            error_function.clear();

            switch(opt.error_function) {
                case 0:
                    error_min = 1;
                    break;
                case 1:
                    error_min = (1.0/number_of_objects_db_load);
                    break;
            }
            out << "# Minimum error calculated to " << error_min; out.endl();  

            d_crit = 0;
            d_crit_step = MET_step_list[j];

            min_d_crit_step = d_crit_step*1e-2;
            max_d_crit_step = d_crit_step*1e2;

            out << "# Running cluster analysis with | step    :" << d_crit_step; out.endl();
            out << "#                               | min step:" << min_d_crit_step; out.endl();
            out << "#                               | max step:" << max_d_crit_step; out.endl();

            
            association_adaptive_step_threshold = opt.adaptive_step_threshold;//0.01;
            step_multiplyer_inc = opt.adaptive_step_mult;

            fraction_associated = 0;
            fraction_associated_next = 0;
    
            association_vec.clear();
            out << "# Starting association loop."; out.endl();
            association_loop_counter = 0;
            fraction_associated_old = 0;
            d_crit_vec.clear();
            TERMINATION_VAR = FALSE;

            d_crit_vec.push_back(d_crit);
            association_vec.push_back(fraction_associated);
            switch(opt.error_function) {
                case 0:
                    error_function.push_back( ((double)number_of_objects_db_load) );
                    break;
                case 1:
                    error_function.push_back( (double)1.0 );
                    break;
            }
            error_function_current = error_function.back();
            error_function_old = error_function_current;

            while(!TERMINATION_VAR) {
                association_loop_counter++;
            
                cluster_matrix.clear();
                fractions.clear();
                
                switch(opt.analysis_type) {
                case 0:
                    cluster_matrix = single_linkage_clustering_matrix(&D_matix, d_crit);
                    break;
                case 1:
                    cluster_matrix = single_linkage_clustering_files(OUTPUT_D_mat_folders[j], list_of_batches_d_matrix, index_matrix, d_crit);
                    break;
                }

                largest_clust_size = 0;
                for(i=0; i < cluster_matrix.size(); i++) {

                    if(cluster_matrix[i].size() > largest_clust_size) {
                        largest_clust_size = cluster_matrix[i].size();
                    }

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

                

                if( ((fraction_associated - fraction_associated_old != 0) || (error_function_current - error_function_old != 0) ) || !opt.adaptive_step ) {
                    d_crit_vec.push_back(d_crit);
                    association_vec.push_back(fraction_associated);
                    switch(opt.error_function) {
                        case 0:
                            error_function.push_back( ((double)number_of_objects_db_load)/((double)largest_clust_size) );
                            break;
                        case 1:
                            error_function.push_back( (double)cluster_matrix.size()/((double)number_of_objects_db_load) );
                            break;
                    }
                    error_function_current = error_function.back();
                }

                // ############# ADAPTIVE STEP #######################
                if(opt.adaptive_step) {

                    fraction_associated_old = fraction_associated;
                    error_function_old = error_function_current;

                    d_crit_old = d_crit;
                    d_crit += d_crit_step;

                    visited_zero_progress_counter = 0;
                    visited_max_progress_counter = 0;
                    association_loop_error_counter = 0;
                    do {
                        association_loop_error_counter++;
                        cluster_matrix_next.clear();
                        fractions_next.clear();

                        
                        switch(opt.analysis_type) {
                            case 0:
                                cluster_matrix_next = single_linkage_clustering_matrix(&D_matix, d_crit);
                                break;
                            case 1:
                                cluster_matrix_next = single_linkage_clustering_files(OUTPUT_D_mat_folders[j], list_of_batches_d_matrix, index_matrix, d_crit);
                                break;
                        }
                        d_crit_used = d_crit;
                        for(i=0; i < cluster_matrix_next.size(); i++) {
                            if((double)cluster_matrix_next[i].size() > 1) {
                                fractions_next.push_back((double)cluster_matrix_next[i].size());
                            }
                        }
                        fractions_next = fractions_next/((double)number_of_objects_db_load);
                        switch(opt.error_function) {
                            case 0:
                                error_function_current = ((double)number_of_objects_db_load)/((double)largest_clust_size);
                                break;
                            case 1:
                                error_function_current = (double)cluster_matrix.size()/((double)number_of_objects_db_load);
                                break;
                        }
                        
                        if(fractions_next.size() > 0) {
                            fraction_associated_next = sum_v(fractions_next);
                        }
                        else {
                            fraction_associated_next = 0;
                        }
                        
                        out << "# Checking next step: " << "frac = " << fraction_associated << ", frac_next = " << fraction_associated_next << ", with D_c = " << d_crit; out.endl();
                        out << "#                   : " << "err  = " << error_function_old << ", err_next  = " << error_function_current; out.endl();
                        if((((fraction_associated_next - fraction_associated) > association_adaptive_step_threshold) ) ) {
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
                        else if((((fraction_associated_next - fraction_associated) == 0) || ( fraction_associated_next == 1 && (error_function_current - error_function_old == 0)) )) {
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

                        switch(opt.parameter_sweep_terminate) {
                            case 0:
                                TERMINATION_VAR = (error_function_current == error_min);
                                break;
                            case 1:
                                TERMINATION_VAR = (fraction_associated == 1);
                                break;
                            case 2:
                                TERMINATION_VAR = (cluster_matrix.size() == 1);
                                break;
                        }
                    
                        if(association_loop_error_counter > 100) {
                            out << "# AN ERROR MIGHT HAVE OCCURED, EXITING LOOP: " << "association_loop_error_counter = " << association_loop_error_counter; out.endl();
                            out << "#                                          : " << "cluster_matrix.size() = " << cluster_matrix.size(); out.endl();
                            out << "# Printing cluster matrix: "; out.endl();
                            print_matrix(cluster_matrix,out);

                            break;
                        }

                    } while(!next_step_ok && !TERMINATION_VAR);

                    if(TERMINATION_VAR) {
                        d_crit_vec.push_back(d_crit_used);
                        association_vec.push_back(fraction_associated_next);
                        error_function.push_back( error_function_current );
                    }

                    out << "Sample associated   : " << fraction_associated*100 << " percent with " << d_crit_old << " threshold."; out.endl();
                    out << "cluster count       : " << cluster_matrix.size(); out.endl();
                    out << "error_function      : " << error_function.back() << ", error_min = " << error_min; out.endl();

                    if(association_loop_error_counter > 10000) {
                        out << "# AN ERROR MIGHT HAVE OCCURED, EXITING LOOP: " << "association_loop_counter = " << association_loop_counter; out.endl();
                        out << "#                                          : " << "cluster_matrix.size() = " << cluster_matrix.size(); out.endl();
                        out << "# Printing cluster matrix: "; out.endl();
                        print_matrix(cluster_matrix,out);

                        break;
                    }
                }
                else {// ############# REGULAR STEP #######################

                    d_crit += d_crit_step;
                    switch(opt.parameter_sweep_terminate) {
                        case 0:
                            TERMINATION_VAR = (error_function_current == error_min);
                            break;
                        case 1:
                            TERMINATION_VAR = (fraction_associated == 1);
                            break;
                        case 2:
                            TERMINATION_VAR = (cluster_matrix.size() == 1);
                            break;
                    }
                }
            }
        
            out << "# Association loop exited."; out.endl();
            out << "# Saving data."; out.endl();

            out << "# Saving err function Dc to file " << OUTPUT_error_files[j]; out.endl();
            temp_mat.clear();
            temp_mat.push_back(d_crit_vec);
            temp_mat.push_back(error_function);
            save_mat(OUTPUT_error_files[j], &temp_mat);
        
            out << "# Saving association curve to file " << OUTPUT_profile_files[j]; out.endl();
            temp_mat.clear(); 
            temp_mat.push_back(d_crit_vec);
            temp_mat.push_back(association_vec);
            save_mat(OUTPUT_profile_files[j], &temp_mat);

        }
        else {
            out << "# Files exist or no calculation scheduled"; out.endl();
        } 
    }
}

if(!opt.parameter_sweep && opt.cluster_analysis) {
    out << "# No parameter sweep is to be performed, using suplied critical values"; out.endl();
    for(j=0; j < FUNCTIONS.size(); j++) {
        if(MET_bool_list[j]) {

            switch(opt.analysis_type) {
                case 0:
                    D_matix.clear();
                    error(load_file(&D_matix,OUTPUT_D_mat_files[j],DATA_MATRIX),matrix_created_files);
                    out << "- Loaded D matrix file: " << OUTPUT_D_mat_files[j]; out.endl();
                    number_of_objects_db_load = D_matix.size()+1;
                    break;
                case 1:
                    list_of_batches_d_matrix.clear();
                    list_of_batches_d_matrix = list_files(OUTPUT_D_mat_folders[j]);
                    index_matrix.clear();
                    for(k=0; k < list_of_batches_d_matrix.size(); k++) {
                        list_of_batches_d_matrix[k].erase(0,OUTPUT_D_mat_folders[j].size()+1);
                        out << "- Index range recorded for file " << list_of_batches_d_matrix[k]; out.endl();

                        char_index_1 = list_of_batches_d_matrix[k].find("_", 0);
                        char_index_2 = list_of_batches_d_matrix[k].find("_", char_index_1+1);

                        index_matrix.push_back(EMPTY_UINT_VEC);
                        buffer_temp = list_of_batches_d_matrix[k].substr(0,char_index_1);
                        index_matrix[k].push_back((unsigned int)strtod(buffer_temp.c_str(), NULL));

                        buffer_temp = list_of_batches_d_matrix[k].substr(char_index_1+1,char_index_2);
                        index_matrix[k].push_back((unsigned int)strtod(buffer_temp.c_str(), NULL));
                    }
                    out << "- Completed index matrix:"; out.endl();
                    print_matrix(index_matrix,out);
                    number_of_objects_db_load = index_matrix.back().back()-1;
                    break;
            }
            out << "# Objects loaded from D file: " << number_of_objects_db_load; out.endl();  
            error_function.clear();

            out << "# Running cluster analysis"; out.endl();
            out << "# Calculating " << FUNCTIONS[j] << " using a critical value of " << MET_calc_list[j]; out.endl();
            cluster_matrix.clear();
            switch(opt.analysis_type) {
                case 0:
                    cluster_matrix = single_linkage_clustering_matrix(&D_matix, MET_calc_list[j]);
                    break;
                case 1:
                    cluster_matrix = single_linkage_clustering_files(OUTPUT_D_mat_folders[j], list_of_batches_d_matrix, index_matrix, MET_calc_list[j]);
                    break;
            }
            
            out << "# Saving cluster matrix."; out.endl();
            temp_mat.clear(); 
            temp_mat = uint_mat_to_double(cluster_matrix);
            save_mat(OUTPUT_cluster_files[j], &temp_mat);
        }
    }
}

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

