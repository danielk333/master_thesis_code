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
    (*opt).ejection_model = (unsigned int)vals[i++];
    (*opt).dt = (double)vals[i++];
    (*opt).orbits = (unsigned int)vals[i++];
    (*opt).particles = (unsigned int)vals[i++];
    (*opt).ejection_number = (unsigned int)vals[i++];
    (*opt).precision = (unsigned int)vals[i++];
    (*opt).output_type = (unsigned int)vals[i++];
    (*opt).model_max_vel = (double)vals[i++];

    return 0;
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

        //std::cout << j << ": " << d_value << std::endl;

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
