#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>


#include "define.hh"
#include "functions.hh"


int load_config(std::string file, OPTIONS *opt) {

    char buffer[256];
    std::vector<double> vals;
    std::string s_value;
    int i, j;

    std::ifstream config;
    config.open(file.c_str(), std::ios::in);

    if(config.fail()) {
        return CONFIG_LOAD_FAILED;
    }

    std::cout << "---------- CONFIGURATION SETTINGS -----------" << std::endl;

    do {
        for(j=0; j < 256; j++) {
            buffer[j] = 32;
        }
        config.getline(&(buffer[0]),256);

        if(buffer[0] != 35 && buffer[0] != 37) {
            i=0;
            while(buffer[i] != 61) {
                i++;
            }

            std::cout << &(buffer[0]) << std::endl;
            vals.push_back(strtod(&(buffer[i+1]), NULL));

        }

    } while(buffer[0] != 37);

    std::cout << "---------- ---------------------- -----------" << std::endl << std::endl;

    config.close();

    i=0;
    (*opt).loops = (unsigned int)vals[i++]; //must be above 1e+02
    (*opt).dt = vals[i++];
    (*opt).sun_collision = vals[i++];
    (*opt).hill_sphere = (unsigned int)vals[i++];
    (*opt).ejection_criteria = vals[i++];
    (*opt).load_old = vals[i++];

    (*opt).integrator = (unsigned int)vals[i++];
    
    (*opt).output_type = (unsigned int)vals[i++];
    (*opt).save_interval = (unsigned int)vals[i++];
    (*opt).coord = (unsigned int)vals[i++];
    (*opt).energy_diag = (unsigned int)vals[i++];
    (*opt).precision = (unsigned int)vals[i++];   

    (*opt).memory_alloc = vals[i++];

    return SIM_OK;
}

int load_previus_settings(char *in, OPTIONS *opt, SIMULATION *sim) {

    char buffer[256];
    std::vector<double> vals;
    std::string s_value;
    int i, j;

    std::ifstream config;
    std::string in_file_name = in;
    config.open(in_file_name.c_str(), std::ios::in);

    if(config.fail()) {
        return CONFIG_LOAD_FAILED;
    }

    std::cout << "---------- CONFIGURATION SETTINGS -----------" << std::endl;

    do {
        for(j=0; j < 256; j++) {
            buffer[j] = 32;
        }
        config.getline(&(buffer[0]),256);

        vals.push_back(strtod(buffer, NULL));

    } while(buffer[0] != 37);

    std::cout << "---------- ---------------------- -----------" << std::endl << std::endl;

    config.close();

    i=0;

    (*opt).loops = (unsigned int)vals[i++];
    (*opt).precision = (unsigned int)vals[i++];
    (*opt).memory_alloc = vals[i++];
    (*opt).output_type = (unsigned int)vals[i++];
    (*opt).save_interval = (unsigned int)vals[i++];
    (*opt).coord = (unsigned int)vals[i++];
    (*opt).energy_diag = (unsigned int)vals[i++];
    (*opt).sun_collision = vals[i++];
    (*opt).hill_sphere = (unsigned int)vals[i++];
    (*opt).ejection_criteria = vals[i++];
    (*opt).integrator = (unsigned int)vals[i++];
    (*opt).dt = vals[i++];

    (*opt).body_n = (unsigned int)vals[i++];
    (*opt).massive_body_n = (unsigned int)vals[i++];
    (*opt).memory_clear_interval = (unsigned int)vals[i++];

    (*sim).loops_performed = (unsigned int)vals[i++];

    (*sim).save_number = (unsigned int)vals[i++];
    (*sim).next_clear = (unsigned int)vals[i++];
    (*sim).diag_counter = (unsigned int)vals[i++];
    (*sim).save_counter = (unsigned int)vals[i++];
    
    (*sim).t_elapsed = vals[i++];

    (*sim).removed_bodies = (unsigned int)vals[i++];

    return SIM_OK;
}

int load_data_vector(std::vector<double> *data,std::string vector_file) {

    std::string buffer;
    std::string buffer_temp;
    std::string s_value;
    std::vector<size_t> space_pos;
    double d_value;
    unsigned int j;

    std::ifstream data_vector;
    data_vector.open(vector_file.c_str(), std::ios::in);

    if(data_vector.fail()) {
        return VECTOR_LOAD_FAILED;
    }

    space_pos.push_back(0);

    while(data_vector.good()) {
        buffer.push_back(data_vector.get());
    }

    while(space_pos.back() < buffer.size()) {
        space_pos.push_back(buffer.find(" ", space_pos.back()+1));
    }
    space_pos.back() = buffer.size();

    for(j=0; j < space_pos.size()-1; j++) {
        buffer_temp = buffer.substr(space_pos[j],space_pos[j+1]);
        d_value = strtod(buffer_temp.c_str(), NULL);
        (*data).push_back(d_value);
    }

    data_vector.close();

    return SIM_OK;
}

int load_data_matrix(std::vector<std::vector<double> > *data,std::string matrix_file) {

    std::string buffer;
    std::string buffer_temp;
    std::string s_value;
    std::vector<unsigned int> space_pos;
    std::vector<double> temp_row;
    unsigned int j;

    std::ifstream data_matrix;
    data_matrix.open(matrix_file.c_str(), std::ios::in);

    if(data_matrix.fail()) {
        return MATRIX_LOAD_FAILED;
    }

    while(!data_matrix.eof()) {
        space_pos.clear(); space_pos.push_back(0);
        temp_row.clear();
        buffer.clear();
        getline(data_matrix,buffer);

        while(space_pos.back() < buffer.size()) {
            space_pos.push_back(buffer.find(" ", space_pos.back()+1));
        }
        space_pos.back() = buffer.size();

        for(j=0; j < space_pos.size()-1; j++) {
            buffer_temp = buffer.substr(space_pos[j],space_pos[j+1]);
            temp_row.push_back(strtod(buffer_temp.c_str(), NULL));
        }

        if (temp_row.size() > 0) {
            (*data).push_back(temp_row);
        }
        
    }

    data_matrix.close();

    return SIM_OK;
}