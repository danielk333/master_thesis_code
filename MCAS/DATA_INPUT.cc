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
        case DATA_FAMILY: 
            STATUS_RETURN =  read_family_matrix((std::vector<std::vector<double> >*)data,file);
            break;
        case DATA_MERCURY: 
            STATUS_RETURN =  read_mercury6_matrix((std::vector<std::vector<double> >*)data,file);
            break;
        case DATA_OPTIONS:
            STATUS_RETURN =  load_config((OPTIONS*)data,file);
            break;
        case DATA_ORBDIST:
            STATUS_RETURN =  read_orbdist_matrix((std::vector<std::vector<double> >*)data,file);
            break;
        case DATA_ORBDIST_SS:
            STATUS_RETURN =  read_orbdist_ss_matrix((std::vector<std::vector<double> >*)data,file);
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
    (*opt).integrator = (unsigned int)vals[i++];
    (*opt).family = (unsigned int)vals[i++];
    (*opt).adaptive_step = (bool)vals[i++];
    (*opt).d_crit_bins = (unsigned int)vals[i++];
    (*opt).mem_allowed = (double)vals[i++];
    (*opt).showers = (unsigned int)vals[i++];
    (*opt).logfile = (unsigned int)vals[i++];

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

int read_mercury6_matrix(std::vector<std::vector<double> > *data,std::string matrix_file) {

    std::string buffer;
    std::string buffer_temp;
    std::string s_value;
    std::vector<unsigned int> space_pos;
    std::vector<double> temp_row;
    unsigned int j,temp_pos1,temp_pos2;

    std::ifstream data_matrix;
    data_matrix.open(matrix_file.c_str(), std::ios::in);

    if(data_matrix.fail()) {
        return mercury6_MATRIX_LOAD_FAILED;
    }
            getline(data_matrix,buffer); //5th line contains data
            getline(data_matrix,buffer);
            getline(data_matrix,buffer);
            getline(data_matrix,buffer);

    while(!data_matrix.eof()) {
        space_pos.clear(); space_pos.push_back(0);
        temp_row.clear();
        buffer.clear();
        getline(data_matrix,buffer);

        temp_pos1 = buffer.find(" ", 0);
        while(space_pos.back() < buffer.size()) {
            temp_pos2 = buffer.find(" ", temp_pos1+1);
            if(temp_pos2 != (temp_pos1+1)) {
                space_pos.push_back(temp_pos2);
            }
            temp_pos1 = temp_pos2;
        }
        space_pos.back() = buffer.size();

        for(j=0; j < space_pos.size()-1; j++) {
            buffer_temp = buffer.substr(space_pos[j],space_pos[j+1]);
            temp_row.push_back(strtod(buffer_temp.c_str(), NULL));
        }
        if(buffer.size() == 93) {
            buffer_temp = buffer.substr(23,26);
            temp_row[1] = strtod(buffer_temp.c_str(), NULL);
        }

        if (temp_row.size() > 0) {
            (*data).push_back(temp_row);
        }
        
    }

    data_matrix.close();

    return CALC_OK;
}

int read_family_matrix(std::vector<std::vector<double> > *data,std::string matrix_file) {

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

if(buffer.size() > 0) {
        while(space_pos.back() < buffer.size()) {
            space_pos.push_back(buffer.find(" ", space_pos.back()+1));
        }
        space_pos.back() = buffer.size();

        if(buffer.substr(space_pos[0],space_pos[1]).compare("U") == 0 ) {
            temp_row.push_back(0);
        }
        else if(buffer.substr(space_pos[0],space_pos[1]).compare("N") == 0 ) {
            temp_row.push_back(1);
        }
        else {
            return DISTRIBUTION_NOT_RECOGNIZED;
        }

        for(j=1; j < space_pos.size()-1; j++) {
            buffer_temp = buffer.substr(space_pos[j],space_pos[j+1]);
            temp_row.push_back(strtod(buffer_temp.c_str(), NULL));
        }

        if (temp_row.size() > 0) {
            (*data).push_back(temp_row);
        }
}
        
    }

    data_matrix.close();

    return CALC_OK;
}


int read_orbdist_matrix(std::vector<std::vector<double> > *data,std::string matrix_file) {

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

if(buffer.size() > 0 && buffer.substr(0,10).compare("synth_name") != 0) {
        while(space_pos.back() < buffer.size()) {
            space_pos.push_back(buffer.find("	", space_pos.back()+1));
        }
        space_pos.back() = buffer.size();

        for(j=1; j < space_pos.size()-1; j++) {
            buffer_temp = buffer.substr(space_pos[j],space_pos[j+1]);
            temp_row.push_back(strtod(buffer_temp.c_str(), NULL));
        }

        if (temp_row.size() > 0) {
            (*data).push_back(temp_row);
        }
}
        
    }

    data_matrix.close();

    return CALC_OK;
}

int read_orbdist_ss_matrix(std::vector<std::vector<double> > *data,std::string matrix_file) {

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

if(buffer.size() > 0 && buffer.substr(0,2).compare("!!") != 0) {
        while(space_pos.back() < buffer.size()) {
            space_pos.push_back(buffer.find(" ", space_pos.back()+1));
        }
        space_pos.back() = buffer.size();

        for(j=2; j < space_pos.size()-2; j++) {
            buffer_temp = buffer.substr(space_pos[j],space_pos[j+1]);
            temp_row.push_back(strtod(buffer_temp.c_str(), NULL));
        }

        if (temp_row.size() > 0) {
            (*data).push_back(temp_row);
        }
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
