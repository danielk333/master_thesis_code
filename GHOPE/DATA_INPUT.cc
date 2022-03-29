#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>


#include "resources.hh"
#include "functions.hh"


void load_config(std::string file, OPTIONS *opt) {

    char buffer[256];
    std::vector<double> vals;
    std::string s_value;
    int i, j, k, l;
    bool only_space;

    std::ifstream config;
    config.open(file.c_str(), std::ios::in);

    if(config.fail()) {
        throw("Faild to open configuration file");
    }

    //std::cout << "---------- CONFIGURATION SETTINGS -----------" << std::endl;

    do {
        for(j=0; j < 256; j++) {
            buffer[j] = 32;
        }
        //Get config line
        config.getline(&(buffer[0]),256);

        //Check if commented line (35 is #) or if end of config (37 is %)
        if(buffer[0] != 35 && buffer[0] != 37) {
            i=0;
            //Find equals sign to take value (61 is =)
            while(buffer[i] != 61) {
                i++;
            }

            if(i > 254) {
                throw("Config load - Could not find equals sign on non-comment line");
            }

            //Find the line break that comes after that (10 is \n and 13 is CR)
            k=i+1;
            while((buffer[k] != 10 || buffer[k] != 13) && k < 255) {
                k++;
            }
            //Check if there is anything else than spaces between = and \n
            only_space = true;
            for(l=i+1; l < k; l++) {
                only_space = only_space & (buffer[l] == 32 || buffer[l] == 0);
            }

            //If there is a option save it in vals
            //Else save the max numerical limit for uint as a identifier of a empy config option
            if(!only_space) {
                vals.push_back(strtod(&(buffer[i+1]), NULL));
            }
            else {
                vals.push_back((double)std::numeric_limits<unsigned int>::max());
            }

            //std::cout << &(buffer[0]) << " | " << vals.back() << std::endl;
            
        }

    } while(buffer[0] != 37);

    //std::cout << "---------- ---------------------- -----------" << std::endl << std::endl;

    config.close();

    // ###########
    // # HERE HANDLE FAULTY CONFIG
    // ###########

    i=0;
    (*opt).steps = (unsigned int)vals[i++];
    (*opt).dt = vals[i++];
    (*opt).start_time = vals[i++];
    (*opt).end_time = vals[i++];
    (*opt).sun_collision_limit = vals[i++];
    (*opt).close_enc_limit = vals[i++];
    (*opt).ejection_limit = vals[i++];
    (*opt).load_old = ((unsigned int)vals[i++]) != 0;
    (*opt).central_body = (unsigned int)vals[i++];
    (*opt).PR_effect = ((unsigned int)vals[i++]) != 0;

    (*opt).integrator = (unsigned int)vals[i++];
    (*opt).ATOL = vals[i++];
    (*opt).RTOL = vals[i++];
    
    (*opt).save_interval_itter = (unsigned int)vals[i++];
    (*opt).save_inteval_time = vals[i++];
    (*opt).save_close = ((unsigned int)vals[i++]) != 0;
    (*opt).save_remove = ((unsigned int)vals[i++]) != 0;
    (*opt).save_log = ((unsigned int)vals[i++]) != 0;
    (*opt).save_summary = ((unsigned int)vals[i++]) != 0;
    (*opt).snapshot_list = ((unsigned int)vals[i++]) != 0;
    (*opt).coord = (unsigned int)vals[i++];
    (*opt).save_energy_diag = ((unsigned int)vals[i++]) != 0;
    (*opt).save_time_diag = ((unsigned int)vals[i++]) != 0;
    (*opt).precision = (unsigned int)vals[i++];

    (*opt).max_RAM = vals[i++];
    (*opt).delim_int = (int)vals[i++];
    (*opt).verbose_level = (unsigned int)vals[i++];
}

void load_data_matrix(std::vector<std::vector<double> > *data,std::string matrix_file, std::string delim) {

    std::string buffer;
    std::string buffer_temp;
    std::string s_value;
    std::vector<unsigned int> space_pos;
    std::vector<double> temp_row;
    unsigned int j;

    std::ifstream data_matrix;
    data_matrix.open(matrix_file.c_str(), std::ios::in);

    if(data_matrix.fail()) {
        throw("Faild to open file to read matrix");
    }

    while(!data_matrix.eof()) {
        space_pos.clear(); space_pos.push_back(0);
        temp_row.clear();
        buffer.clear();
        getline(data_matrix,buffer);

        while(space_pos.back() < buffer.size()) {
            space_pos.push_back(buffer.find(delim, space_pos.back()+1));
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
}
