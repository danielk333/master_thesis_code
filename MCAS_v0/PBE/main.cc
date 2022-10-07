// Gravity simulator main function caller file

// Library includes


//Define global variables and physical constants
#include "define.hh"

//Functions include
#include "functions.hh"

int main (int argc, char *argv[]) {

    std::string config_path;
    if(argc == 2) {
        config_path = argv[1];
    }
    else {
        config_path = "config.cfg";
    }

    std::cout << "Config path chosen to: " << config_path << std::endl;

    #ifdef __linux__
        std::string folder_1 = "INIT/MASSIVE_BODIES/";
        std::string folder_2 = "INIT/PRODUCER/";
        std::string out_data_folder = "OUT_DATA/";
    #endif /* linux */
    #ifdef _WIN32
        std::string folder_1 = "INIT\\MASSIVE_BODIES\\";
        std::string folder_2 = "INIT\\PRODUCER\\";
        std::string out_data_folder = "OUT_DATA\\";
    #endif /* win */

    std::string selfpath = get_selfpath();
    selfpath.erase(selfpath.end()-3,selfpath.end());

    std::vector<unsigned int> N_substep;
    N_substep.push_back(2);
    N_substep.push_back(4);
    N_substep.push_back(6);
    N_substep.push_back(8);
    N_substep.push_back(10);
    N_substep.push_back(12);
    N_substep.push_back(14);
    N_substep.push_back(16);

    /* INIT */
    OPTIONS opt;
    SIMULATION sim;
    sim.done = 0; sim.simulation_output = -1;

    bool ejection_ok;

    error(console_output(&sim,&opt));
    sim.simulation_output++;

    /* initialize random seed: */
    srand ( time(NULL) );

    unsigned int i, j;

    // ######### LOAD EXTERNAL CONFIG ##########
    std::cout << "Loading config from " << config_path << std::endl;
    error(load_config(config_path,&opt));

    std::cout << std::setprecision(opt.precision);

    // ######### OUTPUT CONFIG ##########
    std::cout << "Generating needed paths" << std::endl;

    std::string DEBUG_data_q_out = selfpath + out_data_folder + "debug_pos_data.txt";
    std::string DEBUG_energy_out = selfpath + out_data_folder + "debug_energy_data.txt";

    std::string data_q_out = selfpath + out_data_folder + "particle_q_data.txt";
    std::string data_v_out = selfpath + out_data_folder + "particle_v_data.txt";
    std::string data_m_out = selfpath + out_data_folder + "particle_m_data.txt";
    std::string data_t_out = selfpath + out_data_folder + "particle_t_data.txt";
    std::string abs_v_out = selfpath + out_data_folder + "particle_abs_v_data.txt";
    std::string kep_out = selfpath + out_data_folder + "particle_kep_data.txt";
    std::string data_M_out = selfpath + out_data_folder + "body_m_data.txt";
    std::string body_kep_out = selfpath + out_data_folder + "body_kep_data.txt";
    std::string sim_out = selfpath + out_data_folder + "sim_data.txt";

    std::string mass_name = "mass.data";
    std::string pos_name = "pos.data";
    std::string vel_name = "vel.data";
    /*FORMAT a_1 a_2 a_3 ... from (1e-7g - 1g)*/
    std::string mass_dist_name = "mass_dist.data";

    /*FORMAT initial pos (3d); inital vel (3d); body radius, sigma, alpha, crit radius */
    std::string body_data_name = "body_data.data";
    std::string temp_file_path;

    std::cout << "Calculating time and ejection variables" << std::endl;
    char *time_char;

    time_t rawtime;
    struct tm *timeinfo;

    time(&rawtime);
    timeinfo = localtime(&rawtime);

    time_char = asctime(timeinfo);

    // ####################################

    // ######### LOAD INITIAL CONDITION DATA ##########

    std::cout << "Allocating variables" << std::endl;
    std::vector<unsigned int> particle_ejection_vector;
    std::vector<double> temp_vec, baryC, baryCv;
    std::vector<std::vector<double> > temp_vec3;
    temp_vec3.push_back(temp_vec);
    temp_vec3.push_back(temp_vec);
    std::vector<std::vector<double> > temp_vec5;
    temp_vec5.push_back(temp_vec);
    unsigned int MAX_EJECT_FAILS;
    unsigned int EJECT_FAILS;
    unsigned int EJECT_FAILS_total;
    std::vector<std::string> folders;
    folders.push_back(folder_1);
    folders.push_back(folder_2);
    unsigned int particle_counter;
    unsigned int single_eject_count;

    double init_time;

    std::vector<double> temp_vec_kep;
    std::vector<double> temp_vec_kep_2;

    std::vector<std::vector<double> > temp_mat;
    std::vector<double> TEMP_VEL;
    std::vector<double> dMdt_vec;
    double total_mass_loss;
    double number_of_met_ejected;
    double particle_per_ejection;
    double met_weight;
    double mass_per_ejection;
    unsigned int ejection_number_n;

    unsigned int particles_this_ejection;
    unsigned int total_particles_ejected;
    double total_simulation_time;
    double mass_buildup;
    std::vector<double> dM_vec;
    double vel_mag;
    std::vector<double> m_vector;
    std::vector<double> m_vector_init;
    std::vector<double> abs_v_vector;
    std::vector<std::vector<double> > m_dist_file;
    std::vector<double> m_dist_vector;
    std::vector<double> m_mid_vector;
    std::vector<std::vector<double> > m_vector_save;
    std::vector<int> type_vec;
    std::vector<int> type_vec_init;

    std::vector<double> peri_hel;
    std::vector<double> delta_peri_hel;
    std::vector<double> time_peri_hel;

    bool DEBUG_DATA_OUT;
    bool DEBUG_KEP_OUT;
    std::vector<double> CURRENT_COORD;
    std::vector<double> CURRENT_COORD_init;

    std::vector<std::vector<double> > body_data;
    double v_g;
    std::vector<std::vector<double> > RANDOM_SPHERE_POINTS;

    std::vector<std::vector<double> > kep_save;
    std::vector<std::vector<double> > body_kep_save;
    std::vector<std::vector<double> > v_vec_save;
    std::vector<std::vector<double> > q_vec_save;
    std::vector<double> m_vec_save;
    std::vector<double> t_vec_save;
    std::vector<double> weight_vector_temp;
    std::vector<std::vector<double> > weight_vector_file;

    std::vector<double> q_vec_sync;
    std::vector<double> p_vec_sync;
    double dmdt_tot;
    std::vector<std::vector<double> > p_vec_init;
    std::vector<std::vector<double> > q_vec_init;

    std::vector<std::vector<double> > p_vec;
    std::vector<std::vector<double> > q_vec;
    std::vector<std::vector<double> > q_vec_close;
    std::vector<double> helio_r_vec;
    std::vector<std::vector<double> > p_vec_close;
    unsigned int FILE_COUNTER;

    std::cout << "Loading files." << std::endl << std::endl;
    FILE_COUNTER = 0;

    //## MASSIVE ##
    temp_file_path = selfpath + mass_name.insert(0,folders[0].c_str());
    if(file_exists(temp_file_path)) {
        FILE_COUNTER++;
        std::cout << "Loading data file from " << temp_file_path << std::endl;
        error(load_data_vector(&m_vector,temp_file_path));
    }


    temp_file_path = selfpath + pos_name.insert(0,folders[0].c_str());
    if(file_exists(temp_file_path)) {
        FILE_COUNTER++;
        std::cout << "Loading data file from " << temp_file_path << std::endl;
        error(load_data_matrix(&q_vec,temp_file_path));
    }
    temp_file_path = selfpath + vel_name.insert(0,folders[0].c_str());
    if(file_exists(temp_file_path)) {
        FILE_COUNTER++;
        std::cout << "Loading data file from " << temp_file_path << std::endl;
        error(load_data_matrix(&p_vec,temp_file_path));
    }
    
    for (i = 0; i < m_vector.size(); ++i) {
        type_vec.push_back(1);
    }
    opt.massive_body_n = m_vector.size();
    for (i = 0; i < opt.massive_body_n; ++i) {
        for(j = 0; j < 3; ++j) {
            p_vec[i][j] = p_vec[i][j]*m_vector[i];
        }
    }
    for (i = 0; i < opt.massive_body_n; ++i) {
        for(j = 0; j < 3; ++j) {
            q_vec[i][j] = q_vec[i][j]*AU;
        }
    }

    
    //## MASSLESS ##
    temp_file_path = selfpath + mass_dist_name.insert(0,folders[1].c_str());
    if(file_exists(temp_file_path)) {
        FILE_COUNTER++;
        std::cout << "Loading data file from " << temp_file_path << std::endl;
        error(load_data_matrix(&m_dist_file,temp_file_path));
    }
    temp_file_path = selfpath + body_data_name.insert(0,folders[1].c_str());
    if(file_exists(temp_file_path)) {
        FILE_COUNTER++;
        std::cout << "Loading data file from " << temp_file_path << std::endl;
        error(load_data_matrix(&body_data,temp_file_path));
    }


    m_vector.push_back(body_data[2][4]);

    temp_vec.clear();
    for(j = 0; j < 3; ++j) {
        temp_vec.push_back(body_data[0][j]);
    }
    q_vec.push_back(temp_vec);

    temp_vec.clear();
    for(j = 0; j < 3; ++j) {
        temp_vec.push_back(body_data[1][j]);
    }
    p_vec.push_back(temp_vec);

    opt.body_n = m_vector.size();
    
    for (i = opt.massive_body_n; i < opt.body_n; ++i) {
        for(j = 0; j < 3; ++j) {
            p_vec[i][j] = p_vec[i][j]*m_vector[i];
        }
    }
    type_vec.push_back(0);

    unsigned int PB_index;

    PB_index = m_vector.size() - 1;

    if(FILE_COUNTER < 5) {error(NO_INIT_EXIST);}

    // ####################################
    std::vector<double> v_v;
    std::vector<double> r_v;

    double mu = Ggrav*(m_vector[0]+m_vector[PB_index]);

    v_v.push_back((p_vec[PB_index])[0]/m_vector[PB_index]);
    v_v.push_back((p_vec[PB_index])[1]/m_vector[PB_index]);
    v_v.push_back((p_vec[PB_index])[2]/m_vector[PB_index]);
    r_v.push_back((q_vec[PB_index])[0]);
    r_v.push_back((q_vec[PB_index])[1]);
    r_v.push_back((q_vec[PB_index])[2]);

    std::vector<double> ecc_temp = cross_product(v_v,cross_product(r_v,v_v));
    std::vector<double> ecc_temp2;

    ecc_temp2.push_back(ecc_temp[0]/mu - r_v[0]/abs_v(r_v));
    ecc_temp2.push_back(ecc_temp[1]/mu - r_v[1]/abs_v(r_v));
    ecc_temp2.push_back(ecc_temp[2]/mu - r_v[2]/abs_v(r_v));
    double ecc = abs_v(ecc_temp2);

    double epsilon = dot_product(v_v,v_v)*0.5 - mu/abs_v(r_v);
    double a = -mu/(2*epsilon);
    if(a < 0) {
        a = -a;
        //do something for hyperbolic???
    }
    double peri_dist = (1-ecc)*a;
    double orbit_T = 2*PI*sqrt(pow(a,3)/mu);
    double M_t;
    double appi_dist = (1+ecc)*a;

    // INPUT IN DAYS -> s
    opt.dt = opt.dt*(3600*24);
    opt.loops_orbit = (unsigned int)(orbit_T/(opt.dt));

    double approx_prev_appihelion = orbit_T*0.5;
    double prev_appi_direction = -1.0;
    unsigned int prev_appi_loops;

    time(&(sim.start_t));

	std::vector<std::vector<double> > integrator_coef;
	
    error(load_integrator_scheme(&integrator_coef, &opt));

    double helio_r, v2,s, dMdt, dM, grav_effect;
    double helio_r_m1,helio_r_m2;
    std::vector<double> helio_body_r;
    std::vector<double> helio_body_p;
    double v_dot_r;
    double N;
    unsigned int n_ind;
    std::vector<unsigned int> ejection_occurance_index;

    m_mid_vector = m_dist_file[0];
    m_dist_vector = m_dist_file[1];

    error(console_output(&sim,&opt));
    sim.simulation_output++;

    std::string model_string;
    switch(opt.ejection_model) {
        case 0:
            model_string = "Ma et al. 2002";
            break;
        case 1:
            model_string = "Hughes et al. 2000";
            break;
        case 2:
            model_string = "Whipple 1951";
            break;
        case 3:
            model_string = "Fixed maximum model (linear increase with heliocentric distance)";
            break;
        case 4:
            model_string = "Single ejection, perihelion";
            break;
    }
    
    opt.loops = opt.orbits*opt.loops_orbit;
    DEBUG_DATA_OUT = true;
    DEBUG_KEP_OUT = false;
    std::cout << "-------------------------------" << std::endl;
    std::cout << "      COMET RADIUS            : " << body_data[2][0] << " m" << std::endl;
    std::cout << "      BULK DENSITY            : " << body_data[2][1] << " kg/m^3" << std::endl;
    std::cout << "      ACTIVITY                : " << body_data[2][2] << " " << std::endl;
    std::cout << "      CRITICAL SUBLIM         : " << body_data[2][3]/AU << " AU" << std::endl;
    std::cout << "      MASS                    : " << body_data[2][4] << " kg" << std::endl;
    std::cout << "      PERIHELION PASSAGE DIFF : " << body_data[2][5] << " days" << std::endl;
    std::cout << "-------------------------------" << std::endl;
    std::cout << "      Ejection model          : " << model_string << std::endl;
    std::cout << "      Loops                   : " << opt.loops << " " << std::endl;
    std::cout << "      Numer of orbits         : " << opt.orbits << " " << std::endl;
    std::cout << "      Loops/orbit             : " << opt.loops_orbit << " " << std::endl;
    std::cout << "      Orbit time              : " << orbit_T/(3600*24*365.25) << " y " << std::endl;
    std::cout << "      Time step               : " << opt.dt/(3600*24) << " d " << std::endl;
    std::cout << "      Semi major axis         : " << a/AU << " AU " << std::endl;
    std::cout << "      Perihelion distnace     : " << peri_dist/AU << " AU " << std::endl;
    std::cout << "      Max Ejections           : " << opt.ejection_number << " " << std::endl;
    std::cout << "-------------------------------" << std::endl;
    unsigned int eject_cout = 0;

    unsigned int orbits_of_mass_loss;
    bool performing_mass_loss;

    double tota_time_sync = body_data[2][5]*(3600.0*24.0);
    double desierd_start = body_data[2][6];
    unsigned int loops_to_sync;
    double time_direction;
    double time_direction_2;
    double sync_step_left;
    double steps_sync;
    std::vector<double> current_kep;
    double current_helio_r;
    double last_helio_r;
    double current_time;

    unsigned int peri_hel_distance_n;

    double body_semi_major_axis;

    double sync_dt = opt.dt*1;

    prev_appi_loops = (unsigned int)round(approx_prev_appihelion/sync_dt);
    loops_to_sync = (unsigned int)round(abs_working(tota_time_sync/sync_dt));

    if(body_data[2][5] >= 0) {
        time_direction = 1.0;
    }
    else {
        time_direction = -1.0;
    }

    std::vector<double> T_TOT;
    std::vector<double> ENG_TOT;
    std::vector<double> ANG_TOT;

    T_TOT.push_back(0);
    ENG_TOT.push_back(E_t(q_vec, p_vec, m_vector, type_vec));
    ANG_TOT.push_back(abs_v(total_angular_momentum_qp_heliocentric(&q_vec,&p_vec,&type_vec)));

    std::cout << "STARTING ENERGY      : " << ENG_TOT[0] << " J " << std::endl;
    std::cout << "STARTING ANG MOMENTUM: " << ANG_TOT[0] << " kg m^2 / s " << std::endl << std::endl;

    steps_sync = loops_to_sync*time_direction*sync_dt;
    std::cout << "      Sync loops    : " << loops_to_sync << " " << std::endl << std::endl;

    sync_step_left = tota_time_sync - steps_sync;

    CURRENT_COORD = temp_heliocentric_qp(&q_vec, &p_vec, &m_vector);
    //restore_coord_qp(&q_vec, &p_vec, &m_vector,CURRENT_COORD);

    q_vec_sync.clear();
    p_vec_sync.clear();

    q_vec_sync = q_vec[PB_index];
    p_vec_sync = p_vec[PB_index];

    temp_vec_kep.clear();
    temp_vec_kep = qp_to_kepler(q_vec[PB_index], p_vec[PB_index],m_vector[PB_index],m_vector[0]);
    body_semi_major_axis = temp_vec_kep[0]/AU;

    total_simulation_time = 0;
    std::cout << "## Current relative simulation time: " << total_simulation_time/(3600*24*365.25) << " y" << std::endl<< std::endl;
if(loops_to_sync > 0) {

    std::cout << "Propagating major bodies by " << steps_sync/(3600.0*24.0) << " days" << std::endl;
    for(i=0;i<loops_to_sync;i++) {
        total_simulation_time = total_simulation_time + time_direction*sync_dt;

        restore_coord_qp(&q_vec, &p_vec, &m_vector,CURRENT_COORD);
        //Leapfrog KIN POT 8split method

        BS_step(&q_vec,&p_vec,&m_vector,&type_vec,N_substep,time_direction*sync_dt,body_data[2][1]);
        /*H_kin(&q_vec,&p_vec,&m_vector,&type_vec,time_direction*sync_dt*integrator_coef[0][0]);    
        H_pot_tp_DIS(&q_vec,&p_vec,&m_vector,&type_vec,time_direction*sync_dt*integrator_coef[1][0],1,body_data[2][1]);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,time_direction*sync_dt*integrator_coef[0][0]);*/

        /*H_kin(&q_vec,&p_vec,&m_vector,&type_vec,time_direction*opt.dt*integrator_coef[0][0]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,time_direction*opt.dt*integrator_coef[1][0],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,time_direction*opt.dt*integrator_coef[0][1]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,time_direction*opt.dt*integrator_coef[1][1],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,time_direction*opt.dt*integrator_coef[0][2]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,time_direction*opt.dt*integrator_coef[1][2],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,time_direction*opt.dt*integrator_coef[0][3]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,time_direction*opt.dt*integrator_coef[1][3],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,time_direction*opt.dt*integrator_coef[0][4]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,time_direction*opt.dt*integrator_coef[1][3],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,time_direction*opt.dt*integrator_coef[0][3]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,time_direction*opt.dt*integrator_coef[1][2],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,time_direction*opt.dt*integrator_coef[0][2]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,time_direction*opt.dt*integrator_coef[1][1],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,time_direction*opt.dt*integrator_coef[0][1]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,time_direction*opt.dt*integrator_coef[1][0],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,time_direction*opt.dt*integrator_coef[0][0]);*/
        CURRENT_COORD = temp_heliocentric_qp(&q_vec, &p_vec, &m_vector);
        if(DEBUG_DATA_OUT) {
            temp_vec.clear();
            temp_mat.clear();
            for(j=0; j < q_vec.size(); j++) {
                temp_vec = q_vec[j];
                temp_vec.insert(temp_vec.begin(),(double)j);
                temp_mat.push_back(temp_vec);
            }
            save_mat(DEBUG_data_q_out.c_str(),&temp_mat,&opt);
        }
        
        if(DEBUG_KEP_OUT) {
            body_kep_save.push_back(qp_to_kepler(q_vec[PB_index], p_vec[PB_index], m_vector[PB_index], m_vector[0]));

            (body_kep_save.back())[0] = (body_kep_save.back())[0]/AU;
            (body_kep_save.back()).push_back((sync_dt)*(i*time_direction)/((double)(3600*24)));
            (body_kep_save.back()).push_back(m_vector[PB_index]);
        }

    }

}
    std::cout << "## Current relative simulation time: " << total_simulation_time/(3600*24*365.25) << " y" << std::endl<< std::endl;

    if(abs_working(sync_step_left) > 60.0) {
        total_simulation_time = total_simulation_time + sync_step_left;
        restore_coord_qp(&q_vec, &p_vec, &m_vector,CURRENT_COORD);
        std::cout << "Propagating major bodies by last step " << sync_step_left/(3600.0*24.0) << " days" << std::endl;

        BS_step(&q_vec,&p_vec,&m_vector,&type_vec,N_substep,sync_step_left,body_data[2][1]);

        /*H_kin(&q_vec,&p_vec,&m_vector,&type_vec,sync_step_left*integrator_coef[0][0]);
        H_pot_tp_DIS(&q_vec,&p_vec,&m_vector,&type_vec,sync_step_left*integrator_coef[1][0],1,body_data[2][1]);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,sync_step_left*integrator_coef[0][0]);*/

        /*H_kin(&q_vec,&p_vec,&m_vector,&type_vec,sync_step_left*integrator_coef[0][0]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,sync_step_left*integrator_coef[1][0],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,sync_step_left*integrator_coef[0][1]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,sync_step_left*integrator_coef[1][1],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,sync_step_left*integrator_coef[0][2]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,sync_step_left*integrator_coef[1][2],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,sync_step_left*integrator_coef[0][3]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,sync_step_left*integrator_coef[1][3],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,sync_step_left*integrator_coef[0][4]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,sync_step_left*integrator_coef[1][3],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,sync_step_left*integrator_coef[0][3]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,sync_step_left*integrator_coef[1][2],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,sync_step_left*integrator_coef[0][2]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,sync_step_left*integrator_coef[1][1],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,sync_step_left*integrator_coef[0][1]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,sync_step_left*integrator_coef[1][0],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,sync_step_left*integrator_coef[0][0]);*/
        CURRENT_COORD = temp_heliocentric_qp(&q_vec, &p_vec, &m_vector);
    }

    std::cout << "Resetting PB to complete sync" << std::endl;
    q_vec[PB_index] = q_vec_sync;
    p_vec[PB_index] = p_vec_sync;
    unsigned int correct_peri_find = 0;
    std::cout << "Checking if at wanted first perihelion" << std::endl;
    current_time = desierd_start + tota_time_sync/(3600.0*24.0);
    std::cout << "## Current relative simulation time: " << total_simulation_time/(3600*24*365.25) << " y" << std::endl;
    std::cout << "## Years from desired start        : " << (tota_time_sync)/(3600.0*24.0*365.25) << " y" << std::endl<< std::endl;

    T_TOT.push_back(total_simulation_time);
    ENG_TOT.push_back(E_t(q_vec, p_vec, m_vector, type_vec));
    ANG_TOT.push_back(abs_v(total_angular_momentum_qp_heliocentric(&q_vec,&p_vec,&type_vec)));

if(abs_working(tota_time_sync) > orbit_T/2.0) {

    if(tota_time_sync < 0) {
        time_direction_2 = 1.0;
    }
    else {
        time_direction_2 = -1.0;
    }
    std::cout << "Desired orbits passed: Propagating to correct perihelion" << std::endl;

    helio_body_r = add_v_v(q_vec[PB_index],multiply_v(-1.0,q_vec[0]));
    helio_body_p = add_v_v(p_vec[PB_index],multiply_v(-1.0,p_vec[0]));
    helio_r = abs_v(helio_body_r);
    v_dot_r = dot_product(helio_body_p,helio_body_r);
    current_helio_r = v_dot_r/abs_working(v_dot_r);
    last_helio_r = current_helio_r;
    while( ((last_helio_r*current_helio_r) > 0 && time_direction_2*current_helio_r < 0) || abs_working(current_time - desierd_start) > orbit_T/4.0/(3600.0*24.0) ) {

        last_helio_r = current_helio_r;
        current_time = current_time + time_direction_2*sync_dt/(3600.0*24.0);
        total_simulation_time = total_simulation_time + time_direction_2*sync_dt;
        correct_peri_find = correct_peri_find + 1;
        restore_coord_qp(&q_vec, &p_vec, &m_vector,CURRENT_COORD);
        //Leapfrog KIN POT 8split method

        BS_step(&q_vec,&p_vec,&m_vector,&type_vec,N_substep,time_direction_2*sync_dt,body_data[2][1]);

        /*H_kin(&q_vec,&p_vec,&m_vector,&type_vec,time_direction_2*sync_dt*integrator_coef[0][0]);    
        H_pot_tp_DIS(&q_vec,&p_vec,&m_vector,&type_vec,time_direction_2*sync_dt*integrator_coef[1][0],1,body_data[2][1]);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,time_direction_2*sync_dt*integrator_coef[0][0]);*/

        /*H_kin(&q_vec,&p_vec,&m_vector,&type_vec,time_direction_2*opt.dt*integrator_coef[0][0]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,time_direction_2*opt.dt*integrator_coef[1][0],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,time_direction_2*opt.dt*integrator_coef[0][1]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,time_direction_2*opt.dt*integrator_coef[1][1],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,time_direction_2*opt.dt*integrator_coef[0][2]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,time_direction_2*opt.dt*integrator_coef[1][2],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,time_direction_2*opt.dt*integrator_coef[0][3]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,time_direction_2*opt.dt*integrator_coef[1][3],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,time_direction_2*opt.dt*integrator_coef[0][4]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,time_direction_2*opt.dt*integrator_coef[1][3],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,time_direction_2*opt.dt*integrator_coef[0][3]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,time_direction_2*opt.dt*integrator_coef[1][2],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,time_direction_2*opt.dt*integrator_coef[0][2]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,time_direction_2*opt.dt*integrator_coef[1][1],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,time_direction_2*opt.dt*integrator_coef[0][1]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,time_direction_2*opt.dt*integrator_coef[1][0],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,time_direction_2*opt.dt*integrator_coef[0][0]);*/
        CURRENT_COORD = temp_heliocentric_qp(&q_vec, &p_vec, &m_vector);

        if(DEBUG_DATA_OUT) {
            temp_vec.clear();
            temp_mat.clear();
            for(j=0; j < q_vec.size(); j++) {
                temp_vec = q_vec[j];
                temp_vec.insert(temp_vec.begin(),(double)j);
                temp_mat.push_back(temp_vec);
            }
            save_mat(DEBUG_data_q_out.c_str(),&temp_mat,&opt);
        }

        helio_body_r = add_v_v(q_vec[PB_index],multiply_v(-1.0,q_vec[0]));
        helio_body_p = add_v_v(p_vec[PB_index],multiply_v(-1.0,p_vec[0]));
        helio_r = abs_v(helio_body_r);
        v_dot_r = dot_product(helio_body_p,helio_body_r);
        current_helio_r = v_dot_r/abs_working(v_dot_r);
    }


    std::cout << "Correct perihelion found after " << correct_peri_find << " loops and " << time_direction_2*correct_peri_find*sync_dt/(3600*24*365.25) << " years" << std::endl;
            std::cout << "last_helio_r      :" << last_helio_r << std::endl;
            std::cout << "current_helio_r   :" << current_helio_r << std::endl;
            std::cout << "helio_r           :" << helio_r/AU << std::endl;
}

    helio_r = sqrt(pow((q_vec[PB_index])[0] - (q_vec[0])[0],2) + pow((q_vec[PB_index])[1] - (q_vec[0])[1],2) + pow((q_vec[PB_index])[2] - (q_vec[0])[2],2) );
    std::cout << "Heliocentric distance: " << helio_r/AU << " AU" << std::endl; 

    T_TOT.push_back(total_simulation_time);
    ENG_TOT.push_back(E_t(q_vec, p_vec, m_vector, type_vec));
    ANG_TOT.push_back(abs_v(total_angular_momentum_qp_heliocentric(&q_vec,&p_vec,&type_vec)));

    std::cout << "## Current relative simulation time: " << total_simulation_time/(3600*24*365.25) << " y" << std::endl<< std::endl;

if(prev_appi_loops > 0) {

    helio_body_r = add_v_v(q_vec[PB_index],multiply_v(-1.0,q_vec[0]));
    helio_body_p = add_v_v(p_vec[PB_index],multiply_v(-1.0,p_vec[0]));
    helio_r = abs_v(helio_body_r);
    v_dot_r = dot_product(helio_body_p,helio_body_r);
    current_helio_r = v_dot_r/abs_working(v_dot_r);

    std::cout << "Propagating bodies to reach appihelion" << std::endl;
    for(i=0;i<prev_appi_loops*2;i++) {
        total_simulation_time = total_simulation_time + prev_appi_direction*sync_dt;

        last_helio_r = current_helio_r;
        restore_coord_qp(&q_vec, &p_vec, &m_vector,CURRENT_COORD);
        //Leapfrog KIN POT 8split method

        BS_step(&q_vec,&p_vec,&m_vector,&type_vec,N_substep,prev_appi_direction*sync_dt,body_data[2][1]);

        /*
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,prev_appi_direction*sync_dt*integrator_coef[0][0]);
        H_pot_tp_DIS(&q_vec,&p_vec,&m_vector,&type_vec,prev_appi_direction*sync_dt*integrator_coef[1][0],1,body_data[2][1]);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,prev_appi_direction*sync_dt*integrator_coef[0][0]);*/
        
        /*H_kin(&q_vec,&p_vec,&m_vector,&type_vec,prev_appi_direction*opt.dt*integrator_coef[0][0]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,prev_appi_direction*opt.dt*integrator_coef[1][0],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,prev_appi_direction*opt.dt*integrator_coef[0][1]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,prev_appi_direction*opt.dt*integrator_coef[1][1],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,prev_appi_direction*opt.dt*integrator_coef[0][2]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,prev_appi_direction*opt.dt*integrator_coef[1][2],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,prev_appi_direction*opt.dt*integrator_coef[0][3]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,prev_appi_direction*opt.dt*integrator_coef[1][3],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,prev_appi_direction*opt.dt*integrator_coef[0][4]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,prev_appi_direction*opt.dt*integrator_coef[1][3],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,prev_appi_direction*opt.dt*integrator_coef[0][3]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,prev_appi_direction*opt.dt*integrator_coef[1][2],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,prev_appi_direction*opt.dt*integrator_coef[0][2]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,prev_appi_direction*opt.dt*integrator_coef[1][1],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,prev_appi_direction*opt.dt*integrator_coef[0][1]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,prev_appi_direction*opt.dt*integrator_coef[1][0],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,prev_appi_direction*opt.dt*integrator_coef[0][0]);*/
        CURRENT_COORD = temp_heliocentric_qp(&q_vec, &p_vec, &m_vector);
        if(DEBUG_DATA_OUT) {
            temp_vec.clear();
            temp_mat.clear();
            for(j=0; j < q_vec.size(); j++) {
                temp_vec = q_vec[j];
                temp_vec.insert(temp_vec.begin(),(double)j);
                temp_mat.push_back(temp_vec);
            }
            save_mat(DEBUG_data_q_out.c_str(),&temp_mat,&opt);
        }

        helio_body_r = add_v_v(q_vec[PB_index],multiply_v(-1.0,q_vec[0]));
        helio_body_p = add_v_v(p_vec[PB_index],multiply_v(-1.0,p_vec[0]));
        helio_r = abs_v(helio_body_r);
        v_dot_r = dot_product(helio_body_p,helio_body_r);
        current_helio_r = v_dot_r/abs_working(v_dot_r);
        
        if(DEBUG_KEP_OUT) {
            body_kep_save.push_back(qp_to_kepler(q_vec[PB_index], p_vec[PB_index], m_vector[PB_index], m_vector[0]));

            (body_kep_save.back())[0] = (body_kep_save.back())[0]/AU;
            (body_kep_save.back()).push_back((sync_dt)*(((double)i)*prev_appi_direction)/((double)(3600*24)));
            (body_kep_save.back()).push_back(m_vector[PB_index]);
        }
        

        if( ((last_helio_r*current_helio_r) < 0 && prev_appi_direction*current_helio_r < 0) && ((double)i)*sync_dt > orbit_T*0.1) {
            std::cout << "Appihelion found, exiting loop" << std::endl;
            std::cout << "last_helio_r      :" << last_helio_r << std::endl;
            std::cout << "current_helio_r   :" << current_helio_r << std::endl;
            std::cout << "current time      :" << ((double)i)*sync_dt/(3600*24*365.25) << " y" << std::endl;
            break;
        }

    }

}

    std::cout << "Appihelion loop exited after " << ((double)i)*prev_appi_direction*sync_dt/(3600.0*24.0) << " days" << std::endl;

    helio_r = sqrt(pow((q_vec[PB_index])[0] - (q_vec[0])[0],2) + pow((q_vec[PB_index])[1] - (q_vec[0])[1],2) + pow((q_vec[PB_index])[2] - (q_vec[0])[2],2) );
    std::cout << "Heliocentric distance: " << helio_r/AU << " AU" << std::endl; 

    T_TOT.push_back(total_simulation_time);
    ENG_TOT.push_back(E_t(q_vec, p_vec, m_vector, type_vec));
    ANG_TOT.push_back(abs_v(total_angular_momentum_qp_heliocentric(&q_vec,&p_vec,&type_vec)));

    std::cout << "## Current relative simulation time: " << total_simulation_time/(3600*24*365.25) << " y" << std::endl << std::endl;

    CURRENT_COORD_init = CURRENT_COORD;
    init_time = total_simulation_time;

    std::cout << "Saving initial state for post sublimation sync" << std::endl;
    q_vec_init = q_vec;
    p_vec_init = p_vec;
    m_vector_init = m_vector;
    type_vec_init = type_vec;
    dmdt_tot = 0;
    orbits_of_mass_loss = 0;
    performing_mass_loss = false;
    peri_hel_distance_n = 0;

    helio_r = sqrt(pow((q_vec[PB_index])[0] - (q_vec[0])[0],2) + pow((q_vec[PB_index])[1] - (q_vec[0])[1],2) + pow((q_vec[PB_index])[2] - (q_vec[0])[2],2) );
    std::cout << "Heliocentric distance: " << helio_r/AU << " AU" << std::endl; 

    std::cout << "MAIN LOOPS ENERGY      : " << ENG_TOT.back() << " J " << std::endl;
    std::cout << "MAIN LOOPS ANG MOMENTUM: " << ANG_TOT.back() << " kg m^2 / s " << std::endl << std::endl;

    std::cout << std::setprecision(6) << std::fixed;
    std::cout << "Loop   |Time     |Mass lost   |Mass loss   |Heliocentric distance " << std::endl;
    for(i=0;i<opt.loops;i++) {
        sim.loops_performed = i;

        /*temp_vec.clear();
        temp_mat.clear();
        for(j=0; j < q_vec.size(); j++) {
            temp_vec = q_vec[j];
            temp_vec.insert(temp_vec.begin(),(double)j);
            temp_mat.push_back(temp_vec);
        }
        save_mat(DEBUG_data_q_out.c_str(),&temp_mat,&opt);*/
        //Progressbar
        if(fmod((double)(sim.loops_performed), (double)(opt.loops/10)) < 1) {
            if(dMdt_vec.size() > 0) {
                dmdt_tot = sum_v(dMdt_vec);
            }
            else {
                dmdt_tot = 0;
            }

            std::cout << i << "|" << (i*opt.dt)/(3600*24*365.25)<< " y |" << dmdt_tot << " kg   |" << dMdt << " kg/s |" << helio_r/AU << " AU " << std::endl;
        }
        
        /* PROGRESS TIME */
        //Backup previus step

        //Leapfrog KIN POT 8split method
        restore_coord_qp(&q_vec, &p_vec, &m_vector,CURRENT_COORD);

        BS_step(&q_vec,&p_vec,&m_vector,&type_vec,N_substep,opt.dt,body_data[2][1]);

        /*H_kin(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[0][0]);
        H_pot_tp_DIS(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[1][0],1,body_data[2][1]);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[0][0]);*/

        /*H_kin(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[0][0]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[1][0],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[0][1]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[1][1],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[0][2]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[1][2],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[0][3]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[1][3],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[0][4]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[1][3],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[0][3]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[1][2],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[0][2]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[1][1],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[0][1]);
        H_pot(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[1][0],1);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[0][0]);*/
        CURRENT_COORD = temp_heliocentric_qp(&q_vec, &p_vec, &m_vector);

        if(DEBUG_KEP_OUT) {
            body_kep_save.push_back(qp_to_kepler(q_vec[PB_index], p_vec[PB_index], m_vector[PB_index], m_vector[0]));

            (body_kep_save.back())[0] = (body_kep_save.back())[0]/AU;
            (body_kep_save.back()).push_back((opt.dt)*(sim.loops_performed)/((double)(3600*24)));
            (body_kep_save.back()).push_back(m_vector[PB_index]);
        }

        helio_r = sqrt(pow((q_vec[PB_index])[0] - (q_vec[0])[0],2) + pow((q_vec[PB_index])[1] - (q_vec[0])[1],2) + pow((q_vec[PB_index])[2] - (q_vec[0])[2],2) );
        if(i > 1 && peri_hel_distance_n < opt.orbits && helio_r/AU < body_semi_major_axis) {
            if(helio_r_m2 > helio_r_m1 && helio_r_m1 < helio_r/AU) {
                peri_hel.push_back(helio_r_m1);
                delta_peri_hel.push_back(5*(helio_r_m2 + helio_r/AU - 2.0*helio_r_m1)*0.5);
                peri_hel_distance_n++;
                time_peri_hel.push_back(((opt.dt)*((double)sim.loops_performed) + total_simulation_time)/((double)(3600*24)));
                std::cout << "PERI DISTANCE FOUND: " << peri_hel.back() << " AU +- 5 diff of " << delta_peri_hel.back() << " AU" << std::endl;
            }
            helio_r_m2 = helio_r_m1;
            helio_r_m1 = helio_r/AU;
        }
        else if(i==0) {
            helio_r_m2 = helio_r/AU;
        }
        else if(i==1) {
            helio_r_m1 = helio_r/AU;
        }

        if(helio_r < body_data[2][3]) {
            performing_mass_loss = true;


    switch(opt.ejection_model) {
        case 0:
            dMdt = pow(body_data[2][0],2)*sol_L/(4*body_H)*(1/pow(helio_r,2) - 1/pow(body_data[2][3],2));
            break;
        case 1:
            dMdt = (g_Hug*PI*pow(body_data[2][0],2)*sol_L/(2*PI*pow(helio_r,2))/body_H_W);
            break;
        case 2:
            dMdt = (PI*pow(body_data[2][0],2)*sol_L/(4*PI*pow(helio_r,2))/body_H_W);
            break;
        case 3:
            dMdt = pow(body_data[2][0],2)*sol_L/(4*body_H)*(1/pow(helio_r,2) - 1/pow(body_data[2][3],2));
            break;
        case 4:
            dMdt = pow(body_data[2][0],2)*sol_L/(4*body_H)*(1/pow(helio_r,2) - 1/pow(body_data[2][3],2));
            break;
    }

            
            dMdt_vec.push_back(dMdt);
            dM_vec.push_back(dMdt*opt.dt);

            p_vec[PB_index] = multiply_v( 1 - (dMdt*(opt.dt)/(m_vector[PB_index]) ),p_vec[PB_index]);
            m_vector[PB_index] = m_vector[PB_index] - dMdt*(opt.dt);
            //body_data[2][0] = cbrt(m_vector.back()*3/(4*PI*body_data[2][1]));

            temp_vec.clear();
            temp_vec.push_back(m_vector[PB_index]);
            temp_vec.push_back((opt.dt)*((double)sim.loops_performed)/((double)(3600*24)));
            m_vector_save.push_back(temp_vec);
        }
        else if(dM_vec.size() > 0 && performing_mass_loss) {
            orbits_of_mass_loss++;
            std::cout << "Done ejecting mass for orbit " << orbits_of_mass_loss << std::endl;
            performing_mass_loss = false;
            if(orbits_of_mass_loss == opt.orbits) {
                break;
            }
            
        }
        else {
            performing_mass_loss = false;
        }

    }
    std::cout << std::setprecision(10) << std::scientific;
    total_mass_loss = sum_v(dM_vec);
    particle_per_ejection = (double)floor(((double)opt.particles)/((double)opt.ejection_number));
    mass_per_ejection = total_mass_loss/((double)opt.ejection_number);
    met_weight = (total_mass_loss/dot_product(m_dist_vector,m_mid_vector))/((double)opt.particles);
    
    std::cout << std::endl << "--- Additional output " << std::endl;
    std::cout << "Particle vs meteoroids weight : " << met_weight << std::endl;
    std::cout << "Total mass loss               : " << total_mass_loss << std::endl;    
    std::cout << "Max ejection number           : " << opt.ejection_number << std::endl;
    std::cout << "Mass per ejection             : " << mass_per_ejection << std::endl << std::endl;

            std::cout << std::setprecision(6) << std::fixed;
            std::cout << "# -------------------------------------" << std::endl;
            std::cout << "# PB kepler elements before mass loss: " << std::endl;
            std::cout << "| a      : " << temp_vec_kep[0]/AU << std::endl;
            std::cout << "| e      : " << temp_vec_kep[1] << std::endl;
            std::cout << "| i      : " << temp_vec_kep[2]/PI*180.0 << std::endl;
            std::cout << "| omega  : " << temp_vec_kep[3]/PI*180.0 << std::endl;
            std::cout << "| Omega  : " << temp_vec_kep[4]/PI*180.0 << std::endl << std::endl;
            temp_vec_kep_2 = qp_to_kepler(q_vec[PB_index], p_vec[PB_index],m_vector[PB_index],m_vector[0]);
            std::cout << "# PB kepler elements after mass loss: " << std::endl;
            std::cout << "| a      : " << temp_vec_kep_2[0]/AU << std::endl;
            std::cout << "| e      : " << temp_vec_kep_2[1] << std::endl;
            std::cout << "| i      : " << temp_vec_kep_2[2]/PI*180.0 << std::endl;
            std::cout << "| omega  : " << temp_vec_kep_2[3]/PI*180.0 << std::endl;
            std::cout << "| Omega  : " << temp_vec_kep_2[4]/PI*180.0 << std::endl << std::endl;
            std::cout << std::setprecision(10) << std::scientific;

    q_vec.clear();
    p_vec.clear();
    m_vector.clear();
    type_vec.clear();

    q_vec = q_vec_init;
    p_vec = p_vec_init;
    m_vector = m_vector_init;
    type_vec = type_vec_init;

    std::cout << "Clearing memory" << std::endl;
    kep_save.clear();

    q_vec_init.clear();
    p_vec_init.clear();
    m_vector_init.clear();
    type_vec_init.clear();
    temp_vec.clear();

    std::cout << "### Starting system sync and particle ejection" << std::endl;

    single_eject_count = 0;
    particle_counter = 0;
    performing_mass_loss = false;
    

switch(opt.ejection_model) {
                case 0:
                    MAX_EJECT_FAILS = particle_per_ejection*10;
                    break;
                case 1:
                    MAX_EJECT_FAILS = particle_per_ejection*10;
                    break;
                case 2:
                    MAX_EJECT_FAILS = particle_per_ejection*10;
                    break;
                case 3:
                    MAX_EJECT_FAILS = particle_per_ejection*10;
                    break;
                case 4:
                    MAX_EJECT_FAILS = 1;
                    break;
            }

    mass_buildup = 0;
    ejection_number_n = 0;
    EJECT_FAILS_total = 0;
    total_particles_ejected = 0;
    orbits_of_mass_loss = 0;
    dMdt = 0;
    particle_ejection_vector.clear();
    for(i=0;i<opt.loops;i++) {


        T_TOT.push_back(total_simulation_time);
        ENG_TOT.push_back(E_t(q_vec, p_vec, m_vector, type_vec));
        ANG_TOT.push_back(abs_v(total_angular_momentum_qp_heliocentric(&q_vec,&p_vec,&type_vec)));

        if(fmod((double)(i), (double)(opt.loops/30)) < 1) {
            std::cout << i << "|" << (i*opt.dt)/(3600*24*365.25) << " y |" << particle_counter << " particles | " <<  dMdt << " kg/s |" << helio_r/AU << " AU " << "| dE: " << (ENG_TOT[0] - ENG_TOT.back())/ENG_TOT[0] << std::endl;
        }
        total_simulation_time = total_simulation_time + opt.dt;

        restore_coord_qp(&q_vec, &p_vec, &m_vector,CURRENT_COORD_init);

        BS_step(&q_vec,&p_vec,&m_vector,&type_vec,N_substep,opt.dt,body_data[2][1]);

        /*H_kin(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[0][0]);
        H_pot_tp_DIS(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[1][0],1,body_data[2][1]);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[0][0]);*/

        /*
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[0][0]);
        H_pot_tp_DIS(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[1][0],1,body_data[2][1]);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[0][1]);
        H_pot_tp_DIS(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[1][1],1,body_data[2][1]);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[0][2]);
        H_pot_tp_DIS(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[1][2],1,body_data[2][1]);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[0][3]);
        H_pot_tp_DIS(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[1][3],1,body_data[2][1]);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[0][4]);
        H_pot_tp_DIS(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[1][3],1,body_data[2][1]);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[0][3]);
        H_pot_tp_DIS(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[1][2],1,body_data[2][1]);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[0][2]);
        H_pot_tp_DIS(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[1][1],1,body_data[2][1]);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[0][1]);
        H_pot_tp_DIS(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[1][0],1,body_data[2][1]);
        H_kin(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[0][0]);*/

        CURRENT_COORD_init = temp_heliocentric_qp(&q_vec, &p_vec, &m_vector);
        if(DEBUG_DATA_OUT) {
            temp_mat.clear();
            for(j=0; j < q_vec.size(); j++) {
                temp_vec.clear();
                temp_vec = q_vec[j];
                temp_vec.insert(temp_vec.begin(),(double)j);
                temp_mat.push_back(temp_vec);
            }
            save_mat(DEBUG_data_q_out.c_str(),&temp_mat,&opt);
        }

        dMdt = 0;
        helio_r = sqrt(pow((q_vec[PB_index])[0] - (q_vec[0])[0],2) + pow((q_vec[PB_index])[1] - (q_vec[0])[1],2) + pow((q_vec[PB_index])[2] - (q_vec[0])[2],2) );
        if(helio_r < body_data[2][3]) {
            performing_mass_loss = true;

            switch(opt.ejection_model) {
                case 0:
                    dMdt = pow(body_data[2][0],2)*sol_L/(4*body_H)*(1/pow(helio_r,2) - 1/pow(body_data[2][3],2));
                    break;
                case 1:
                    dMdt = (g_Hug*PI*pow(body_data[2][0],2)*sol_L/(2*PI*pow(helio_r,2))/body_H_W);
                    break;
                case 2:
                    dMdt = (PI*pow(body_data[2][0],2)*sol_L/(4*PI*pow(helio_r,2))/body_H_W);
                    break;
                case 3:
                    dMdt = pow(body_data[2][0],2)*sol_L/(4*body_H)*(1/pow(helio_r,2) - 1/pow(body_data[2][3],2));
                    break;
                case 4:
                    dMdt = pow(body_data[2][0],2)*sol_L/(4*body_H)*(1/pow(helio_r,2) - 1/pow(body_data[2][3],2));
                    break;
            }
            p_vec[PB_index] = multiply_v( 1 - (dMdt*(opt.dt)/(m_vector[PB_index]) ),p_vec[PB_index]);
            m_vector[PB_index] = m_vector[PB_index] - dMdt*(opt.dt);
            //body_data[2][0] = cbrt(m_vector.back()*3/(4*PI*body_data[2][1]));
        }
        else if(performing_mass_loss) {
            orbits_of_mass_loss++;
            std::cout << "Done ejecting mass for orbit " << orbits_of_mass_loss << std::endl;
            performing_mass_loss = false;
            if(orbits_of_mass_loss == opt.orbits) {
                break;
            }
        }

        mass_buildup = mass_buildup + dMdt*(opt.dt);
        
        if(mass_buildup >= mass_per_ejection) {

            grav_effect = 2*Ggrav*(m_vector[PB_index])/body_data[2][0];

            switch(opt.ejection_model) {
                case 0:
                    v2 = water_W*dMdt/(2*PI*body_data[2][2]*body_data[2][1]*body_data[2][0]);
                    break;
                case 1:
                    v_g = sqrt(8*k_B*T_g_Hug/(PI*mu_Hug));
                    v2 = 2*dMdt/(Psi_Hug*PI*body_data[2][0])*(theta_Hug*v_g)*xi_Hug;
                    break;
                case 2:
                    v_g = sqrt(8*k_B*T_g0/(PI*mu_Wip))/pow(helio_r/AU,tau_Wip);
                    v2 = dMdt*K_drag*v_g/(2*PI*body_data[2][0]);
                    break;
                case 3:
                    v2 = helio_r*(-opt.model_max_vel/(body_data[2][0] - peri_dist)) - body_data[2][0]*(-opt.model_max_vel/(body_data[2][0] - peri_dist));
                    break;
                case 4:
                    v2 = opt.model_max_vel;
                    break;
            }
            EJECT_FAILS = 0;

            particles_this_ejection = (unsigned int)floor((mass_buildup/dot_product(m_dist_vector,m_mid_vector))/met_weight);
            if(total_particles_ejected + particles_this_ejection > opt.particles) {
                particles_this_ejection = opt.particles - total_particles_ejected;
            }
            if(particles_this_ejection > 0) {
                particle_ejection_vector.push_back(particles_this_ejection);
                q_vec_save.clear();
                RANDOM_SPHERE_POINTS.clear();
                RANDOM_SPHERE_POINTS = generate_random_sphere_normals(particles_this_ejection);
                q_vec_save = multiply_mat(body_data[2][0],RANDOM_SPHERE_POINTS);

                v_vec_save.clear();
                kep_save.clear();
                t_vec_save.clear();
                m_vec_save.clear();
                abs_v_vector.clear();

                for(j=0; j < particles_this_ejection; j++) {
                    ejection_ok = FALSE;

                    do {
                        n_ind = draw_from_dist(m_dist_vector);
                        s = cbrt(m_mid_vector[n_ind]*3/(4*PI*body_data[2][1]));
                        switch(opt.ejection_model) {
                            case 0:
                                vel_mag = v2/s - grav_effect;
                                if(vel_mag > 0) {
                                    ejection_ok = TRUE;
                                    vel_mag = sqrt(vel_mag);
                                    EJECT_FAILS = 0;
                                }
                                else {
                                    EJECT_FAILS++;
                                    EJECT_FAILS_total++;
                                }
                                break;
                            case 1:
                                vel_mag = v2/(s*body_data[2][1])*(3.0/4.0) - grav_effect;
                                if(vel_mag > 0) {
                                    ejection_ok = TRUE;
                                    vel_mag = sqrt(vel_mag);
                                    EJECT_FAILS = 0;
                                }
                                else {
                                    EJECT_FAILS++;
                                    EJECT_FAILS_total++;
                                }
                                break;
                            case 2:
                                vel_mag = v2/(s*body_data[2][1])*(3.0/4.0) - grav_effect;
                                if(vel_mag > 0) {
                                    ejection_ok = TRUE;
                                    vel_mag = sqrt(vel_mag);
                                    EJECT_FAILS = 0;
                                }
                                else {
                                    EJECT_FAILS++;
                                    EJECT_FAILS_total++;
                                }
                                break;
                            case 3:
                                vel_mag = v2;
                                ejection_ok = TRUE;
                                break;
                            case 4:
                                if(single_eject_count == 0 && helio_r/AU < (peri_hel[orbits_of_mass_loss] + delta_peri_hel[orbits_of_mass_loss])) {
                                    vel_mag = v2;
                                    ejection_ok = TRUE;
                                }
                                else {
                                    EJECT_FAILS++;
                                    EJECT_FAILS_total++;
                                }

                                
                                break;
                        }

                    } while(!ejection_ok && EJECT_FAILS < MAX_EJECT_FAILS);
                    
                    if(EJECT_FAILS < MAX_EJECT_FAILS) {

                        TEMP_VEL = multiply_v(vel_mag,RANDOM_SPHERE_POINTS[j]);

                        abs_v_vector.push_back(abs_v(TEMP_VEL));
                        v_vec_save.push_back(TEMP_VEL);

                        m_vec_save.push_back(m_mid_vector[n_ind]);
                        t_vec_save.push_back(((double)i)*opt.dt);
                        kep_save.push_back(xv_to_kepler(add_v_v(q_vec[PB_index],q_vec_save[j]), add_v_v(multiply_v(1/(m_vector[PB_index]),p_vec[PB_index]),v_vec_save[j]), m_vector[0]));
                        kep_save[j][0] = kep_save[j][0]/AU;

                        temp_vec3[0].push_back(helio_r/AU);
                        temp_vec3[1].push_back(abs_v_vector[j]);

                        
                    }

                }

                if(EJECT_FAILS < MAX_EJECT_FAILS) {
                    ejection_number_n++;
                    single_eject_count++;
                    total_particles_ejected = total_particles_ejected + particles_this_ejection;

                    for(j=0; j < particles_this_ejection; j++) {
                        q_vec.push_back(add_v_v(q_vec[PB_index],q_vec_save[j]));
                        p_vec.push_back(multiply_v(m_vec_save[j],add_v_v(multiply_v(1.0/m_vector[PB_index],p_vec[PB_index]),v_vec_save[j])));
                        m_vector.push_back(m_vec_save[j]);
                        type_vec.push_back(0);

                        particle_counter++;
                    }

                    temp_vec5[0] = t_vec_save;
                    error(save_mat(kep_out.c_str(), &kep_save, &opt));
                    error(save_mat(data_t_out.c_str(), &temp_vec5, &opt));
                    
                    mass_buildup = 0;
                }
            }

        }

    }
    number_of_met_ejected = total_particles_ejected*met_weight;

    std::cout << "Number of meteroids ejected   : " << number_of_met_ejected << std::endl;
    std::cout << "Ejection number               : " << ejection_number_n << std::endl;
    std::cout << "Particle number               : " << total_particles_ejected << std::endl;
    std::cout << "Ejections failed              : " << EJECT_FAILS_total << std::endl;

    error(save_mat(abs_v_out.c_str(), &temp_vec3, &opt));
    error(save_mat(data_M_out.c_str(), &m_vector_save, &opt));
    if(DEBUG_KEP_OUT) {
        error(save_mat(body_kep_out.c_str(), &body_kep_save, &opt));
    }
    error(save_mat(kep_out.c_str(), &kep_save, &opt));

    std::cout << "Particles added to (q,p) data       : " << particle_counter << std::endl;
    std::cout << "System sync complete to T_0 + " << total_simulation_time/(3600*24*365.25) << " years" << std::endl;

    std::cout << "END ENERGY       : " << ENG_TOT.back() << " J " << std::endl;
    std::cout << "END ANG MOMENTUM : " << ANG_TOT.back() << " kg m^2 / s " << std::endl;
    std::cout << "ENERGY DIFF      : " << (ENG_TOT[0] - ENG_TOT.back())/ENG_TOT[0] << " " << std::endl;
    std::cout << "ANG MOMENTUM DIFF: " << (ANG_TOT[0] - ANG_TOT.back())/ANG_TOT[0] << " " << std::endl << std::endl;

    temp_mat.clear();
    temp_mat.push_back(T_TOT);
    temp_mat.push_back(ENG_TOT);
    temp_mat.push_back(ANG_TOT);
    error(save_mat(DEBUG_energy_out.c_str(), &temp_mat, &opt));

    std::cout << std::setprecision(6) << std::fixed;
    std::cout << "# -------------------------------------" << std::endl;
    std::cout << "# PB kepler elements before ejections: " << std::endl;
    std::cout << "| a      : " << temp_vec_kep[0]/AU << std::endl;
    std::cout << "| e      : " << temp_vec_kep[1] << std::endl;
    std::cout << "| i      : " << temp_vec_kep[2]/PI*180.0 << std::endl;
    std::cout << "| omega  : " << temp_vec_kep[3]/PI*180.0 << std::endl;
    std::cout << "| Omega  : " << temp_vec_kep[4]/PI*180.0 << std::endl << std::endl;
    temp_vec_kep.clear();
    temp_vec_kep = qp_to_kepler(q_vec[PB_index], p_vec[PB_index],m_vector[PB_index],m_vector[0]);
    std::cout << "# PB kepler elements after ejections: " << std::endl;
    std::cout << "| a      : " << temp_vec_kep[0]/AU << std::endl;
    std::cout << "| e      : " << temp_vec_kep[1] << std::endl;
    std::cout << "| i      : " << temp_vec_kep[2]/PI*180.0 << std::endl;
    std::cout << "| omega  : " << temp_vec_kep[3]/PI*180.0 << std::endl;
    std::cout << "| Omega  : " << temp_vec_kep[4]/PI*180.0 << std::endl << std::endl;

    std::cout << "Writing simulation data" << std::endl;

    weight_vector_temp.clear();
    weight_vector_temp.push_back(met_weight); // particle representation of this many mets 0
    weight_vector_temp.push_back(total_simulation_time); // total simulations time         1
    weight_vector_temp.push_back(body_data[2][1]); // bulk density                         2
    weight_vector_temp.push_back(m_vector[PB_index]); // comet mass                        3
    weight_vector_temp.push_back(init_time); // time sync                                  4
    weight_vector_temp.push_back((double)opt.orbits);//orbits                              5
    for(i=0; i < opt.orbits; i++) {
        weight_vector_temp.push_back(time_peri_hel[i]);//perihelion_passages days relative sim start
    }
    weight_vector_file.clear();
    weight_vector_file.push_back(weight_vector_temp);

    error(save_mat(sim_out.c_str(), &weight_vector_file, &opt));

    remove(data_m_out.c_str());
    remove(data_q_out.c_str());
    remove(data_v_out.c_str());

    
    std::cout << "Converting momentum to velocity " << std::endl;
    for(i=0; i < p_vec.size(); i++) {
        p_vec[i] = multiply_v(1.0/m_vector[i],p_vec[i]);
    }

    std::cout << "Saving synced files" << std::endl;
    temp_vec5.clear();
    temp_vec5.push_back(temp_vec);
    temp_vec5[0] = m_vector;
    error(save_mat(data_m_out.c_str(),&temp_vec5,&opt));
    error(save_mat(data_q_out.c_str(),&q_vec,&opt));
    error(save_mat(data_v_out.c_str(),&p_vec,&opt));

    time(&(sim.end_t));
    sim.t_elapsed = difftime(sim.end_t,sim.start_t);

    sim.simulation_output++;
    error(console_output(&sim,&opt));
    std::cout << "Partices created              : " << total_particles_ejected << std::endl;

    return SIM_OK;
}

void error(int SIM_STATUS) {

    if(SIM_STATUS != SIM_OK) {
        std::cout << std::endl << "ERROR ENCOUTERD: " << SIM_STATUS << std::endl << std::endl;
        exit(SIM_STATUS);
    }
}


