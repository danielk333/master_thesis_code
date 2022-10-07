#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <typeinfo>

#include "define.hh"
#include "functions.hh"

int load_file(const void *data, std::string file, int type) {
    if(type != DATA_MATRIX_SILENT) {
        std::cout << "Loading data from: " << file << std::endl;
        std::cout << "Load type: " << type << std::endl;
    }

    int STATUS_RETURN = CALC_OK;
    if(!file_exists(file)) {
        return FILE_DOES_NOT_EXIST;
    }

    switch(type) {
        case DATA_VECTOR: 
            STATUS_RETURN =  load_data_vector((std::vector<double>*)data,file);
            break;
        case DATA_MATRIX: 
            STATUS_RETURN =  load_data_matrix((std::vector<std::vector<double> >*)data,file);
            break;  
        case DATA_MATRIX_SILENT: 
            STATUS_RETURN =  load_data_matrix((std::vector<std::vector<double> >*)data,file);
            break;  
        case DATA_OPTIONS:
            STATUS_RETURN =  load_config((OPTIONS*)data,file);
            break;
        case DATA_OPTIONS_VEC:
            STATUS_RETURN =  load_config_vec((OPTIONS*)data,file);
            break;
        default:
            STATUS_RETURN =  FAILD_TO_DETECT_DATA_LOAD_TYPE;
    }
    if(type != DATA_MATRIX_SILENT) {
        std::cout << "Load complete with return " << STATUS_RETURN << std::endl << std::endl;
    }

    return STATUS_RETURN;
}

int load_config(OPTIONS *opt, std::string file) {

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
    (*opt).analysis_type = (unsigned int)vals[i++];
    (*opt).input_format = (unsigned int)vals[i++];
    (*opt).relative_obj_row = (unsigned int)vals[i++];
    (*opt).cluster_analysis = (bool)vals[i++];
    (*opt).cluster_merge_function = (unsigned int)vals[i++];
    (*opt).parameter_sweep = (bool)vals[i++];
    (*opt).adaptive_step = (bool)vals[i++];
    (*opt).adaptive_step_threshold = (double)vals[i++];
    (*opt).adaptive_step_mult = (double)vals[i++];
    (*opt).parameter_sweep_terminate = (unsigned int)vals[i++];
    (*opt).error_function = (unsigned int)vals[i++];
    (*opt).mem_allowed = (double)vals[i++];
    (*opt).logfile = (unsigned int)vals[i++];
    (*opt).D_SH = (bool)vals[i++];
    (*opt).D_SH_Crit = (double)vals[i++];
    (*opt).D_SH_step = (double)vals[i++];
    (*opt).D_D = (bool)vals[i++];
    (*opt).D_D_Crit = (double)vals[i++];
    (*opt).D_D_step = (double)vals[i++];
    (*opt).rho2 = (bool)vals[i++];
    (*opt).rho2_Crit = (double)vals[i++];
    (*opt).rho2_step = (double)vals[i++];
    (*opt).varrho1 = (bool)vals[i++];
    (*opt).varrho1_Crit = (double)vals[i++];
    (*opt).varrho1_step = (double)vals[i++];

    return CALC_OK;
}


int load_config_vec(OPTIONS *opt, std::string file) {

    std::string buffer;
    std::string buffer_temp;
    std::string s_value;
    std::vector<size_t> space_pos;
    double d_value;
    unsigned int j,i;
    std::vector<double> vals;

    std::ifstream data_vector;
    data_vector.open(file.c_str(), std::ios::in);

    if(data_vector.fail()) {
        return CONFIG_LOAD_FAILED;
    }

    space_pos.push_back(0);

    while(data_vector.good()) {
        buffer.push_back(data_vector.get());
    }

    while(space_pos.back() < buffer.size()) {
        space_pos.push_back(buffer.find(" ", space_pos.back()+1));
    }
    space_pos.back() = buffer.size();
    std::cout << "---------- CONFIGURATION SETTINGS -----------" << std::endl;
    for(j=0; j < space_pos.size()-1; j++) {
        buffer_temp = buffer.substr(space_pos[j],space_pos[j+1]);
        d_value = strtod(buffer_temp.c_str(), NULL);

        std::cout << d_value << " ";

        vals.push_back(d_value);
    }
    std::cout << std::endl;
    data_vector.close();

    std::cout << "---------- ---------------------- -----------" << std::endl << std::endl;

    i=0;
    (*opt).analysis_type = (unsigned int)vals[i++];
    (*opt).input_format = (unsigned int)vals[i++];
    (*opt).relative_obj_row = (unsigned int)vals[i++];
    (*opt).cluster_analysis = (bool)vals[i++];
    (*opt).cluster_merge_function = (unsigned int)vals[i++];
    (*opt).parameter_sweep = (bool)vals[i++];
    (*opt).adaptive_step = (bool)vals[i++];
    (*opt).adaptive_step_threshold = (double)vals[i++];
    (*opt).adaptive_step_mult = (double)vals[i++];
    (*opt).parameter_sweep_terminate = (unsigned int)vals[i++];
    (*opt).error_function = (unsigned int)vals[i++];
    (*opt).mem_allowed = (double)vals[i++];
    (*opt).logfile = (unsigned int)vals[i++];
    (*opt).D_SH = (bool)vals[i++];
    (*opt).D_SH_Crit = (double)vals[i++];
    (*opt).D_SH_step = (double)vals[i++];
    (*opt).D_D = (bool)vals[i++];
    (*opt).D_D_Crit = (double)vals[i++];
    (*opt).D_D_step = (double)vals[i++];
    (*opt).rho2 = (bool)vals[i++];
    (*opt).rho2_Crit = (double)vals[i++];
    (*opt).rho2_step = (double)vals[i++];
    (*opt).varrho1 = (bool)vals[i++];
    (*opt).varrho1_Crit = (double)vals[i++];
    (*opt).varrho1_step = (double)vals[i++];

    return CALC_OK;
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

    return CALC_OK;
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

    return CALC_OK;
}


int fetch_rows_from_matrix_file(std::vector<std::vector<double> > *data, std::string matrix_file, unsigned int n, unsigned int m) {

    std::string buffer;
    std::string buffer_temp;
    std::vector<unsigned int> space_pos;
    std::vector<double> temp_row;
    unsigned int j;

    std::ifstream data_matrix;
    data_matrix.open(matrix_file.c_str(), std::ios::in);

    if(data_matrix.fail()) {
        return MATRIX_LOAD_FAILED;
    }
    for(j=0; j < n; j++) {
        getline(data_matrix,buffer);
    }
    for(j=n; j < m; j++) {
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

        if(temp_row.size() > 0) {
            (*data).push_back(temp_row);
        }
    }

    data_matrix.close();

    return CALC_OK;
}
