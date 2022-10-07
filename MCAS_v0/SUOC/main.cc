// Gravity simulator main function caller file

// Library includes


//Define global variables and physical constants
#include "define.hh"

//Functions include
#include "functions.hh"

int main (int argc, char *argv[]) {

    unsigned int clones;
    std::string in_format;
    std::string out_format;
    std::string number;
    if(argc > 1) {
        in_format = argv[1];
    }
    else {
        in_format = "kep";
    }
    if(argc > 2) {
        number = argv[2];
        clones = (unsigned int)round(strtod(number.c_str(),NULL));
    }
    else {
        clones = 1000;
    }

    if(argc > 3) {
        out_format = argv[3];
    }
    else {
        out_format = "kep";
    }

    std::string selfpath = get_selfpath();
    selfpath.erase(selfpath.end()-4,selfpath.end());

    /* initialize random seed: */
    srand ( time(NULL) );

    unsigned int i, j;

    // ######### OUTPUT CONFIG ##########
    std::cout << "Generating needed paths and allocating memory... " << std::endl;
    std::string orb_in = selfpath + "orb_in.data";
    std::string SUOC_config = selfpath + "SUOC_config.cfg";

    std::vector<std::vector<double> > CONFIG_EXTRA;
    std::cout << "Loading extra config data from file: " << SUOC_config << std::endl;
    load_data_matrix(&CONFIG_EXTRA,SUOC_config);
    // x y z
    // sigma_x sigma_y sigma_z
    // vx vy vz
    // sigma_vx sigma_vy sigma_vz
    //3d normal dist, (x,y,z) = mu_vec
    //                (sigma_vx,sigma_vy,sigma_vz) = sigma_vec
    // rho = 0? then 3 separate normal dists

    //std::string data_q_out = selfpath + "particle_q_data.txt";
    //std::string data_v_out = selfpath + "particle_v_data.txt";
    std::string orb_out = selfpath + "orb_clones.txt";
    char *time_char;

    srand (time(NULL));

    time_t rawtime;
    struct tm *timeinfo;

    time(&rawtime);
    timeinfo = localtime(&rawtime);

    time_char = asctime(timeinfo);

    // ####################################

    std::vector<double> temp_vec;
    std::vector<std::vector<double> > temp_mat;

    // ####################################
    std::cout << "Loading body data from file: " << orb_in << std::endl;
    std::vector<std::vector<double> > orb_data;
    load_data_matrix(&orb_data,orb_in);
    unsigned int dist_res = (unsigned int)CONFIG_EXTRA[0][1];
    double sigma_range = CONFIG_EXTRA[0][0];


    // ####################################
    // ####################################
        std::vector<double> a_dist;
    std::vector<double> e_dist;
    std::vector<double> i_dist;
    std::vector<double> omega_dist;
    std::vector<double> Omega_dist;
    std::vector<double> tp_dist;

    std::vector<double> a_dist_mid;
    std::vector<double> e_dist_mid;
    std::vector<double> i_dist_mid;
    std::vector<double> omega_dist_mid;
    std::vector<double> Omega_dist_mid;
    std::vector<double> tp_dist_mid;

    std::vector<std::vector<double> > kep_out;

        std::vector<double> x_dist;
    std::vector<double> y_dist;
    std::vector<double> z_dist;

    std::vector<double> vx_dist;
    std::vector<double> vy_dist;
    std::vector<double> vz_dist;

    std::vector<double> x_dist_mid;
    std::vector<double> y_dist_mid;
    std::vector<double> z_dist_mid;

    std::vector<double> vx_dist_mid;
    std::vector<double> vy_dist_mid;
    std::vector<double> vz_dist_mid;

    std::vector<std::vector<double> > q_out;
    std::vector<std::vector<double> > v_out;

    // ####################################
    // ####################################


    if(in_format.compare("cart") == 0) {

    x_dist = normal_distribution_density_vector(orb_data[0][0],orb_data[1][0],sigma_range,dist_res);
    y_dist = normal_distribution_density_vector(orb_data[0][1],orb_data[1][1],sigma_range,dist_res);
    z_dist = normal_distribution_density_vector(orb_data[0][2],orb_data[1][2],sigma_range,dist_res);

    vx_dist = normal_distribution_density_vector(orb_data[2][0],orb_data[3][0],sigma_range,dist_res);
    vy_dist = normal_distribution_density_vector(orb_data[2][1],orb_data[3][1],sigma_range,dist_res);
    vz_dist = normal_distribution_density_vector(orb_data[2][2],orb_data[3][2],sigma_range,dist_res);

    x_dist_mid = normal_distribution_bin_mid_vector(orb_data[0][0],orb_data[1][0],sigma_range,dist_res);
    y_dist_mid = normal_distribution_bin_mid_vector(orb_data[0][1],orb_data[1][1],sigma_range,dist_res);
    z_dist_mid = normal_distribution_bin_mid_vector(orb_data[0][2],orb_data[1][2],sigma_range,dist_res);

    vx_dist_mid = normal_distribution_bin_mid_vector(orb_data[2][0],orb_data[3][0],sigma_range,dist_res);
    vy_dist_mid = normal_distribution_bin_mid_vector(orb_data[2][1],orb_data[3][1],sigma_range,dist_res);
    vz_dist_mid = normal_distribution_bin_mid_vector(orb_data[2][2],orb_data[3][2],sigma_range,dist_res);

    for(i=0; i < clones; i++) {
        temp_vec.clear();
        temp_vec.push_back(x_dist_mid[draw_from_dist(x_dist)]);
        temp_vec.push_back(y_dist_mid[draw_from_dist(y_dist)]);
        temp_vec.push_back(z_dist_mid[draw_from_dist(z_dist)]);
        q_out.push_back(temp_vec);

        temp_vec.clear();
        temp_vec.push_back(vx_dist_mid[draw_from_dist(vx_dist)]);
        temp_vec.push_back(vy_dist_mid[draw_from_dist(vy_dist)]);
        temp_vec.push_back(vz_dist_mid[draw_from_dist(vz_dist)]);
        v_out.push_back(temp_vec);
    }
    }
    else if(in_format.compare("kep") == 0) {


    // ####################################

    a_dist = normal_distribution_density_vector(orb_data[0][0],orb_data[1][0],sigma_range,dist_res);
    e_dist = normal_distribution_density_vector(orb_data[0][1],orb_data[1][1],sigma_range,dist_res);
    i_dist = normal_distribution_density_vector(orb_data[0][2],orb_data[1][2],sigma_range,dist_res);
    omega_dist = normal_distribution_density_vector(orb_data[0][3],orb_data[1][3],sigma_range,dist_res);
    Omega_dist = normal_distribution_density_vector(orb_data[0][4],orb_data[1][4],sigma_range,dist_res);
    tp_dist = normal_distribution_density_vector(orb_data[0][5],orb_data[1][5],sigma_range,dist_res);

    a_dist_mid = normal_distribution_bin_mid_vector(orb_data[0][0],orb_data[1][0],sigma_range,dist_res);
    e_dist_mid = normal_distribution_bin_mid_vector(orb_data[0][1],orb_data[1][1],sigma_range,dist_res);
    i_dist_mid = normal_distribution_bin_mid_vector(orb_data[0][2],orb_data[1][2],sigma_range,dist_res);
    omega_dist_mid = normal_distribution_bin_mid_vector(orb_data[0][3],orb_data[1][3],sigma_range,dist_res);
    Omega_dist_mid = normal_distribution_bin_mid_vector(orb_data[0][4],orb_data[1][4],sigma_range,dist_res);
    tp_dist_mid = normal_distribution_bin_mid_vector(orb_data[0][5],orb_data[1][5],sigma_range,dist_res);

    for(i=0; i < clones; i++) {
        temp_vec.clear();
        temp_vec.push_back(a_dist_mid[draw_from_dist(a_dist)]);
        temp_vec.push_back(e_dist_mid[draw_from_dist(e_dist)]);
        temp_vec.push_back(i_dist_mid[draw_from_dist(i_dist)]);
        temp_vec.push_back(omega_dist_mid[draw_from_dist(omega_dist)]);
        temp_vec.push_back(Omega_dist_mid[draw_from_dist(Omega_dist)]);
        temp_vec.push_back(tp_dist_mid[draw_from_dist(tp_dist)]);
        kep_out.push_back(temp_vec);
    }
    }
    

    if(out_format.compare("kep") == 0) {
        //a e i argperi Omega 
        if(in_format.compare("kep") == 0) {
            save_mat(orb_out, &kep_out);
        }
        else if(in_format.compare("cart") == 0) {
            temp_mat.clear();
            for(i=0; i < clones; i++) {
                temp_vec.clear();
                temp_vec = xv_to_kepler(q_out[i], v_out[i], Msol);
                temp_vec.erase(temp_vec.begin() + 5, temp_vec.end());
                temp_mat.push_back(temp_vec);
            }

            save_mat(orb_out, &temp_mat);
        }

    }
    
    return CALC_OK;
}

void error(int SIM_STATUS) {

    if(SIM_STATUS != CALC_OK) {
        std::cout << std::endl << "ERROR ENCOUTERD: " << SIM_STATUS << std::endl << std::endl;
        exit(SIM_STATUS);
    }
}


