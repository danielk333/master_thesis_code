#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "define.hh"
#include "functions.hh"


int save_mat(std::string out_file_name, std::vector<std::vector<double> > *M) {

        unsigned int i;
        unsigned int j;

        std::ofstream out_pos(out_file_name.c_str(), std::ios::app);
        out_pos << std::setprecision(10) << std::scientific;
        if(out_pos.fail()) {
            return DATA_SAVE_FAILED;
        }

        for(j = 0; j < (*M).size(); j++) {
        	for(i = 0; i < (*M)[j].size(); i++) {
        		out_pos << (((*M)[j])[i]);
                if(i < ((*M)[j].size() - 1)) {
                    out_pos << " ";
                }
                else {
                    out_pos << "\r\n";
                }
        	}
        }

        out_pos.close();

        return CALC_OK;
}

int save_batch(std::string file_name, std::string temp_file_name, std::vector<std::vector<double> > *M, unsigned int row_skip) {
        unsigned int i;
        unsigned int j;
    std::string buffer;
    std::string buffer_temp;

    std::ifstream file;
    std::ofstream temp_file;

    std::ofstream file_write;
    std::ifstream temp_file_read;

    file.open(file_name.c_str(), std::ios::in);
    if(file.fail()) {
        return MATRIX_LOAD_FAILED;
    }

    temp_file.open(temp_file_name.c_str(), std::ofstream::out | std::ofstream::trunc);
    temp_file << std::setprecision(10) << std::scientific;
    if(temp_file.fail()) {
        return DATA_SAVE_FAILED;
    }

    for(j=0; j < row_skip; j++) {
        buffer.clear();
        getline(file,buffer);
        if(buffer.size() > 0) {
            temp_file << buffer;
        }
    }

    for(j = 0; j < (*M).size(); j++) {
        buffer.clear();
        getline(file,buffer);

        if(buffer.size() > 0) {
            buffer.erase(buffer.end()-1,buffer.end());

            temp_file << buffer << " ";
            for(i = 0; i < (*M)[j].size(); i++) {
                temp_file << (((*M)[j])[i]);
                if(i < (*M)[j].size()-1) {
                    temp_file << " ";
                }
                else {
                    temp_file << "\r\n";
                }
            }
        }
    }
    
    while(!file.eof()) {
        buffer.clear();
        getline(file,buffer);
        if(buffer.size() > 0) {
            temp_file << buffer;
        }
    }

    temp_file.close();
    file.close();

    file_write.open(file_name.c_str(), std::ofstream::out | std::ofstream::trunc);
    if(file_write.fail()) {
        return DATA_SAVE_FAILED;
    }

    temp_file_read.open(temp_file_name.c_str(),  std::ios::in);
    if(temp_file_read.fail()) {
        return MATRIX_LOAD_FAILED;
    }

    while(!temp_file_read.eof()) {
        buffer.clear();
        getline(temp_file_read,buffer);
        if(buffer.size() > 0) {
            file_write << buffer;
        }
    }
    temp_file_read.close();
    file_write.close();

    return CALC_OK;
}
