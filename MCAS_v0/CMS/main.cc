// Gravity simulator main function caller file

// Library includes
#include <string>

//Define global variables and physical constants
#include "define.hh"

//Functions include
#include "functions.hh"

int main (int argc, char *argv[]) {


    std::string the_one = "1";
    /* INIT */
    OPTIONS opt;
    SIMULATION sim;
    sim.done = 0; sim.simulation_output = -1; sim.diag_counter = 0; sim.save_counter = 0;

    error(console_output(&sim,&opt));
    sim.simulation_output++;

    std::string selfpath = get_selfpath();
    selfpath.erase(selfpath.end()-3,selfpath.end());

    std::string config_path;
    if(argc >= 2) {
        config_path = argv[1];
    }
    else {
        config_path = selfpath + "config.cfg";
    }

    std::vector<std::string> folders;

#ifdef __linux__
    std::string folder_1 = "INIT/MASSIVE_BODIES/";
    std::string folder_2 = "INIT/MASSLESS_BODIES/";
    std::string out_data_folder = "OUT_DATA/";
#endif /* linux */
#ifdef _WIN32
    std::string folder_1 = "INIT\\MASSIVE_BODIES\\";
    std::string folder_2 = "INIT\\MASSLESS_BODIES\\";
#endif /* win */
    folders.push_back(folder_1);
    folders.push_back(folder_2);

    /* initialize random seed: */
    srand ( time(NULL) );

    unsigned int i, j;

    // ######### LOAD EXTERNAL CONFIG ##########
    std::cout << "Loading config from " << config_path << std::endl;
    error(load_config(config_path,&opt));

    std::cout << std::setprecision(opt.precision);

    // ######### OUTPUT CONFIG ##########
    char *time_char;

    time_t rawtime;
    struct tm *timeinfo;

    time(&rawtime);
    timeinfo = localtime(&rawtime);

    time_char = asctime(timeinfo);
    std::string data_out;
    std::string diag_data_out;
    std::string settings_out;
    std::string kep_out_name;
    std::string pos_out_name;

    data_out = "body_data.txt";
    diag_data_out = "energy_diag.txt";
    settings_out = "sim_settings.txt";
    kep_out_name = "kep_out.txt";
    pos_out_name = "pos_out.txt";

    data_out = selfpath + out_data_folder + data_out;
    diag_data_out = selfpath + out_data_folder + diag_data_out;
    settings_out = selfpath + out_data_folder + settings_out;
    kep_out_name = selfpath + out_data_folder + kep_out_name;
    pos_out_name = selfpath + out_data_folder + pos_out_name;
    // ####################################

    std::vector<double> KEP_TEMP;
    std::vector<std::vector<double> > KEP_MAT_TEMP;

    // ######### LOAD INITIAL CONDITION DATA ##########

    std::vector<double> temp_vec;

    std::vector<double> m_vector;
    std::vector<std::vector<double> > prop_mat;
    std::vector<double> E_vec;
    std::vector<double> T_vec;
    std::vector<int> type_vec;

    std::vector<std::vector<std::vector<double> > > p_vec_save;
    std::vector<std::vector<std::vector<double> > > q_vec_save;
    std::vector<std::vector<double> > m_vec_save;
    std::vector<std::vector<int> > type_vec_save;
    std::vector<double> t_vec_save;

    std::vector<std::vector<double> > p_vec;
    std::vector<std::vector<double> > q_vec;

    std::vector<std::vector<double> > P_vec;
    std::vector<std::vector<double> > Q_vec;

    std::vector<std::vector<double> > p_vec_temp;
    std::vector<std::vector<double> > q_vec_temp;

    std::vector<std::vector<double> > r_mat;

    std::string mass_name =  "mass.data"; //for major bodies, 1 row, each col one obj mass
    std::string prop_name =  "prop.data"; // first row: mass, second row: density
    std::string pos_name = "pos.data"; //for all bodies n rows, 3 cols
    std::string vel_name = "vel.data"; // for all bodies n rows, 3 cols
    std::string temp_file_path;

    unsigned int FILE_COUNTER = 0;


    //## MASSIVE ##
    temp_file_path = selfpath + folders[0] + mass_name;
    if(file_exists(temp_file_path)) {
        FILE_COUNTER++;
        std::cout << "Loading data file from " << temp_file_path << std::endl;
        error(load_data_vector(&m_vector,temp_file_path));
    }
    temp_file_path = selfpath + folders[0] + pos_name;
    if(file_exists(temp_file_path)) {
        FILE_COUNTER++;
        std::cout << "Loading data file from " << temp_file_path << std::endl;
        error(load_data_matrix(&q_vec,temp_file_path));
    }
    temp_file_path = selfpath + folders[0] + vel_name;
    if(file_exists(temp_file_path)) {
        FILE_COUNTER++;
        std::cout << "Loading data file from " << temp_file_path << std::endl;
        error(load_data_matrix(&p_vec,temp_file_path));
    }

    for (i = 0; i < m_vector.size(); i++) {
        m_vector[i] = m_vector[i]/pow(G_gauss,2)*Msol;
        type_vec.push_back(1);
    }
    
    //std::cout << "m: " << m_vector[0] << ", q: " << abs_v(q_vec[3])/AU << ", v: " << abs_v(p_vec[3])*1e-3 << std::endl;

    if(m_vector.size() != q_vec.size() || q_vec.size() != p_vec.size() || m_vector.size() != p_vec.size()) {
        std::cout << "Warning: massive body initial condition sizes non cosistent." << std::endl;
        std::cout << "Position size: " << q_vec.size() << std::endl;
        std::cout << "Velocity size: " << p_vec.size() << std::endl;
        std::cout << "Mass size    : " << m_vector.size() << std::endl;
    }

    opt.massive_body_n = m_vector.size();

    if(FILE_COUNTER < 3) {error(NO_INIT_EXIST);}

    for (i = 0; i < opt.massive_body_n; ++i) {
        for(j = 0; j < 3; ++j) {
            p_vec[i][j] = p_vec[i][j]*m_vector[i];
        }
    }

    //## MASSLESS ##
    temp_file_path = selfpath +  folders[1] + prop_name;
    if(file_exists(temp_file_path)) {
        FILE_COUNTER++;
        std::cout << "Loading data file from " << temp_file_path << std::endl;
        error(load_data_matrix(&prop_mat,temp_file_path));
    }
    if(prop_mat.size() > 0) {
        for(i=0; i < prop_mat[0].size(); i++) {
            m_vector.push_back(prop_mat[0][i]);
        }
    }
    opt.body_n = m_vector.size();
    for(i = type_vec.size(); i < opt.body_n; ++i) {
        type_vec.push_back(0);
    }
    temp_file_path = selfpath +  folders[1] + pos_name;
    if(file_exists(temp_file_path)) {
        FILE_COUNTER++;
        std::cout << "Loading data file from " << temp_file_path << std::endl;
        error(load_data_matrix(&q_vec,temp_file_path));
    }
    for (i = 0; i < opt.body_n; ++i) {
        for(j = 0; j < 3; ++j) {
            q_vec[i][j] = q_vec[i][j]*AU;
        }
    }
    temp_file_path = selfpath +  folders[1] + vel_name;
    if(file_exists(temp_file_path)) {
        FILE_COUNTER++;
        std::cout << "Loading data file from " << temp_file_path << std::endl;
        error(load_data_matrix(&p_vec,temp_file_path));
    }
    temp_file_path = selfpath +  folders[1] + prop_name;
    if(file_exists(temp_file_path)) {
        FILE_COUNTER++;
        std::cout << "Loading data file from " << temp_file_path << std::endl;
        error(load_data_matrix(&prop_mat,temp_file_path));
    }
    for (i = opt.massive_body_n; i < opt.body_n; ++i) {
        for(j = 0; j < 3; ++j) {
            p_vec[i][j] = p_vec[i][j]*m_vector[i];
        }
    }
    unsigned int MAX_TEST_PARTICLES = 10;
    if(opt.body_n - opt.massive_body_n > MAX_TEST_PARTICLES) {
        q_vec.erase(q_vec.begin() + opt.massive_body_n + MAX_TEST_PARTICLES, q_vec.end());
        p_vec.erase(p_vec.begin() + opt.massive_body_n + MAX_TEST_PARTICLES, p_vec.end());
        m_vector.erase(m_vector.begin() + opt.massive_body_n + MAX_TEST_PARTICLES, m_vector.end());
        type_vec.erase(type_vec.begin() + opt.massive_body_n + MAX_TEST_PARTICLES, type_vec.end());
        opt.body_n = m_vector.size();
    }

    // ####################################
    error(calculate_memory_management(&opt));

    sim.save_number = (unsigned int)(opt.loops/opt.save_interval)+1;
    sim.next_clear = opt.memory_clear_interval;

    time(&(sim.start_t));

    if (opt.energy_diag == 1) {
        E_vec.push_back(E_t(&q_vec,&p_vec,&m_vector,&type_vec));
        T_vec.push_back(0);
        sim.diag_counter++;
    }
	

	std::vector<std::vector<double> > integrator_coef;
	
    error(load_integrator_scheme(&integrator_coef, &opt));
    
    std::vector<unsigned int> N_substep;
    N_substep.push_back(2);
    N_substep.push_back(4);
    N_substep.push_back(6);
    N_substep.push_back(8);
    N_substep.push_back(10);
    N_substep.push_back(12);
    N_substep.push_back(14);
    N_substep.push_back(16);

    error(console_output(&sim,&opt));
    sim.simulation_output++;

    std::cout << "Loop      |Elapsed time |Estimated time left |Simulation time " << std::endl;

    for(i=0; i< opt.loops; i++) {
    	sim.loops_performed = i;
        

        //Progressbar
        if((unsigned int)((opt.loops/10)*sim.done) == i) {
            error(console_output(&sim,&opt));
            //std::cout << "Loop: " << i << std::endl;
        }
        
        /* PROGRESS TIME */
        //Backup previus step
        q_vec_temp = q_vec;
        p_vec_temp = p_vec;

        switch(opt.integrator) {
        case 0: //Leapfrog KIN POT method
        	H_kin(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[0][0]);
        	H_pot(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[1][0],1);
        	H_kin(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[0][0]);

        	break;
        case 1: //Leapfrog KIN POT 8split method
            H_kin(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[0][0]);
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
            H_kin(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[0][0]);

            break;
        case 2: //Leapfrog KEP INT method
            H_pot(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[0][0],1);
            H_kep(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[1][0],opt.sun_collision,opt.ejection_criteria);
            H_pot(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[0][0],1);

            break;
        case 3: //DH method
            H_dri(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[0][0]);
            H_pot(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[1][0],0);
            H_kep(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[2][0],opt.sun_collision,opt.ejection_criteria);
            H_pot(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[1][0],0);
            H_dri(&q_vec,&p_vec,&m_vector,&type_vec,opt.dt*integrator_coef[0][0]);
            
            break;
        case 4: //BS-Method
            BS_step(&q_vec,&p_vec,&m_vector,&type_vec,N_substep,opt.dt);
            break;
        }

        //Check for collisions
        check_coll(&q_vec, &p_vec, &q_vec_temp, &p_vec_temp, &m_vector, &type_vec);

        /*if(i == 3) {
            break;
        }*/

        /* SAVE DATA FOR OUTPUT LATER */
        if((unsigned int)((double)sim.save_counter*(double)opt.save_interval) == i) {
            /*switch(opt.integrator) {
            case 2:
                jacobi_inv(&q_vec, &p_vec, &m_vector);
                break;
            case 3:
                canonical_heliocentric_inv(&q_vec, &p_vec, &m_vector);
                break;
            }

            error(add_save_data_to_memory(&opt,&q_vec_save, &p_vec_save, &q_vec, &p_vec, &m_vector));
            */
            q_vec_save.push_back(q_vec);
            p_vec_save.push_back(p_vec);
            m_vec_save.push_back(m_vector);
            type_vec_save.push_back(type_vec);
            t_vec_save.push_back((i+1)*opt.dt/(3600.0*24.0));

            for(j=1; j < q_vec.size(); j++) {
                KEP_TEMP.clear();
                KEP_TEMP = qp_to_kepler(q_vec[j], p_vec[j], m_vector[j], m_vector[0]);
                KEP_TEMP.insert(KEP_TEMP.begin(),(opt.dt)*((double)i)/(3600*24));
                KEP_TEMP.insert(KEP_TEMP.begin(),(double)j);
                KEP_MAT_TEMP.clear();
                KEP_MAT_TEMP.push_back(KEP_TEMP);
                save_mat(kep_out_name,&KEP_MAT_TEMP);

                KEP_TEMP.clear();
                KEP_TEMP = q_vec[j];
                KEP_TEMP.insert(KEP_TEMP.begin(),(opt.dt)*((double)i)/(3600*24));
                KEP_TEMP.insert(KEP_TEMP.begin(),(double)j);
                KEP_MAT_TEMP.clear();
                KEP_MAT_TEMP.push_back(KEP_TEMP);
                save_mat(pos_out_name,&KEP_MAT_TEMP);
            }


            if (opt.energy_diag == 1) {
                E_vec.push_back(E_t(&q_vec,&p_vec,&m_vector,&type_vec));
                T_vec.push_back((i+1)*opt.dt/(3600*24*365));
                sim.diag_counter++;
            }
            /*
            switch(opt.integrator) {
            case 2:
                jacobi(&q_vec, &p_vec, &m_vector);
                break;
            case 3:
                canonical_heliocentric(&q_vec, &p_vec, &m_vector);
                break;
            }*/
            sim.save_counter++;
        }

        /* CLEAR MEMORY AND SAVE TO FILE */
        if(i == sim.next_clear) {
            sim.next_clear += opt.memory_clear_interval;
            error(save(data_out.c_str(), &q_vec_save, &p_vec_save, &m_vec_save, &type_vec_save, &t_vec_save, &opt, &sim));

            if (opt.energy_diag == 1) {
                error(save_diagnostics(diag_data_out.c_str(),E_vec,T_vec));
                E_vec.clear();
                T_vec.clear();
            }

            q_vec_save.clear();
            p_vec_save.clear();
            m_vec_save.clear();
            type_vec_save.clear();
            t_vec_save.clear();
        }

    }
    /*
    switch(opt.integrator) {
    case 2:
        jacobi_inv(&q_vec, &p_vec, &m_vector);
        break;
    case 3:
        canonical_heliocentric_inv(&q_vec, &p_vec, &m_vector);
        break;
    }*/

    //error(add_save_data_to_memory(&opt,&q_vec_save, &p_vec_save, &q_vec, &p_vec, &m_vector));
            q_vec_save.push_back(q_vec);
            p_vec_save.push_back(p_vec);
            m_vec_save.push_back(m_vector);
            type_vec_save.push_back(type_vec);
            t_vec_save.push_back((i+1)*opt.dt/(3600.0*24.0));

    if (opt.energy_diag == 1) {
        E_vec.push_back(E_t(&q_vec,&p_vec,&m_vector,&type_vec));
        T_vec.push_back((i+1)*opt.dt/(3600*24*365));
        sim.diag_counter++;
    }
    sim.save_counter++;
    
    error(save(data_out.c_str(), &q_vec_save, &p_vec_save, &m_vec_save, &type_vec_save, &t_vec_save, &opt, &sim));
    error(save_diagnostics(diag_data_out.c_str(),E_vec,T_vec));
    error(save_settings(settings_out.c_str(),&opt, &sim));
    q_vec_save.clear();
    p_vec_save.clear();

    sim.removed_bodies = 0;
    for(i=0; i<type_vec.size(); i++) {
        if(type_vec[i] == -1) {sim.removed_bodies++;}
    }

    time(&(sim.end_t));
    sim.t_elapsed = difftime(sim.end_t,sim.start_t);

    sim.simulation_output++;
    error(console_output(&sim,&opt));
    
    return SIM_OK;
}

void error(int SIM_STATUS) {

    if(SIM_STATUS != SIM_OK) {
        std::cout << std::endl << "ERROR ENCOUTERD: " << SIM_STATUS << std::endl << std::endl;
        exit(SIM_STATUS);
    }
}


