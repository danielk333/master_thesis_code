/*
 * main.cc
 *
 *  Created on: Feb 17, 2016
 *      Author: dankas
 */

//Simulation Merger (SM)

#include <time.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string>
#include <vector>
#include <unistd.h>
#include <limits.h>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <typeinfo>
#include <cfloat>
#include <cstdio>
#include <memory>

#include <utility> //pair

#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include "boost/progress.hpp"

#include <sys/types.h>
#include <unistd.h>

void remove_folder(const std::string& name);
std::vector<std::string> list_files(const std::string& folder);
int load_data_matrix(std::vector<std::vector<double> > *data,std::string matrix_file);
int save_mat(std::string out_file_name, std::vector<std::vector<double> > *M);

int main (int argc, char *argv[]) {
	unsigned int i,ui,j,k;

	std::cout << "#######################################" << std::endl;
	std::cout << "# --Simulation Merger program --" << std::endl;
	std::cout << "# Daniel Kastinen, Feb 2016" << std::endl;
	std::cout << "#######################################" << std::endl << std::endl;

    std::string path_a;
    std::string path_b;
    std::string path_c;
    if(argc == 4) {
        path_a = argv[1];
        path_b = argv[2];
        path_c = argv[3];

        std::cout << "- Merging: " << path_a << std::endl;
        std::cout << "- and    : " << path_a << std::endl;
        std::cout << "- INTO   : " << path_c << std::endl;
    }
    else {
        std::cout << "NOT CORRECT AMMOUNT OF ARGUMENTS" << std::endl;
    }
    //FUNCTIONS TO USE
    std::vector<std::string> FUNCTIONS;
    FUNCTIONS.push_back("SH");
    FUNCTIONS.push_back("D");
    FUNCTIONS.push_back("rho2");
    FUNCTIONS.push_back("varrho1");

    std::vector<std::string> statistics_association_profile;
    std::vector<std::string> statistics_error_profile;
    std::vector<std::string> statistics_function_divergance;

    std::vector<std::string> file_list;
    unsigned int id_n,no_id_n;
    id_n = 0;
    no_id_n = 0; 
    //id

    std::string PB_kep_timeS_data_out = "PB_kep_timeS_data.txt";;//check
    file_list.push_back(PB_kep_timeS_data_out);id_n++;
    std::string time_data_out = "time_data.txt";;//check
    file_list.push_back(time_data_out);id_n++;
    std::string earth_encounter_data_out = "earth_encounter_data.txt";//check
    file_list.push_back(earth_encounter_data_out);id_n++;
    std::string statistics_PB_n_out = "ejector_n_data.txt";//check
    file_list.push_back(statistics_PB_n_out);id_n++;
    for(ui=0; ui < FUNCTIONS.size(); ui++) {
        statistics_association_profile.push_back( "association_profile_" + FUNCTIONS[ui] + ".txt");//check
        file_list.push_back(statistics_association_profile[ui]);id_n++;
        statistics_error_profile.push_back("error_profile_" + FUNCTIONS[ui] + ".txt");//check
        file_list.push_back(statistics_error_profile[ui]);id_n++;
        statistics_function_divergance.push_back("function_divergance_" + FUNCTIONS[ui] + ".txt");//check
        file_list.push_back(statistics_function_divergance[ui]);id_n++;
    }
    std::string met_before_encounter_data_out = "met_before_encounter_data.txt";//check
    file_list.push_back(met_before_encounter_data_out);id_n++;
    std::string met_encounter_data_out = "met_encounter_data.txt"; //check
    file_list.push_back(met_encounter_data_out);id_n++;
    std::string PBE_all_m_statistics = "particle_ejection_m_stat.txt";//check
    file_list.push_back(PBE_all_m_statistics);id_n++;
    std::string PBE_abs_v_statistics = "particle_ejection_v_stat.txt";//check
    file_list.push_back(PBE_abs_v_statistics);id_n++;
    std::string statistics_all_PB_kep_out = "PB_kep_data.txt";//check
    file_list.push_back(statistics_all_PB_kep_out);id_n++;


    //no id

    std::string statistics_PB_char_out = "ejector_char_data.txt";//check
    file_list.push_back(statistics_PB_char_out);no_id_n++;
    std::string statistics_PB_kep_out = "ejector_kep_data.txt";//check
    file_list.push_back(statistics_PB_kep_out);no_id_n++;
    std::string PBE_weight_statistics = "particle_weight_stat.txt";//check
    file_list.push_back(PBE_weight_statistics);no_id_n++;
    std::string PB_kepler_dist = "PB_kep_dist.txt";//check
    file_list.push_back(PB_kepler_dist);no_id_n++;
    std::string PBE_M_statistics = "PBE_M_stat.txt";//check
    file_list.push_back(PBE_M_statistics);no_id_n++;


    //std::string log_file = "outlog.txt";//check
    //file_list.push_back(log_file); 

    std::vector<std::string> list_a;
    std::vector<std::string> list_b;

    list_a = list_files(path_a);
    list_b = list_files(path_b);

    if(list_a.size() != list_b.size()) {
    	std::cout << "|-- Warning extra files detected --|" << std::endl;
    }
    if(no_id_n + id_n != file_list.size() || (file_list.size() != list_a.size() || file_list.size() != list_b.size() ) ) {
      std::cout << "|-- Warning file count anomalus --|" << std::endl;
    }
    std::string temp_path_x;
    std::string temp_path_c;
    boost::filesystem::path dir(path_c);
    if(boost::filesystem::create_directory(dir)) {
        std::cout << "# Folder create success: " << path_c << std::endl << std::endl;
    }
    int RET_status;
    //use statistics_all_PB_kep_out to make id list
    std::vector<std::vector<double> > TEMP_DATA;
    std::vector<unsigned int> ID_LIST_1;
    std::vector<unsigned int> ID_LIST_2;
    std::vector<unsigned int> ID_LIST_2_mod;
    TEMP_DATA.clear();
    RET_status = load_data_matrix(&TEMP_DATA,path_a + statistics_all_PB_kep_out);
    std::cout << "# ID list for : " << path_a << std::endl;
    for(i=0; i < TEMP_DATA.size(); i++) {
      ID_LIST_1.push_back(TEMP_DATA[i][0]);
      std::cout << ID_LIST_1[i] << " ";
    }
    std::cout << std::endl << std::endl;

    TEMP_DATA.clear();
    RET_status = load_data_matrix(&TEMP_DATA,path_b + statistics_all_PB_kep_out);
    std::cout << "# ID list for : " << path_b << std::endl;
    for(i=0; i < TEMP_DATA.size(); i++) {
      ID_LIST_2.push_back(TEMP_DATA[i][0]);
      ID_LIST_2_mod.push_back(TEMP_DATA[i][0]);
      std::cout << ID_LIST_2[i] << " ";
    }
    TEMP_DATA.clear();
    std::cout << std::endl << std::endl;

    unsigned int MAX_ID_a = ID_LIST_1.back();
    std::cout << "# New ID list for : " << path_c << std::endl;
    for(i=0; i < ID_LIST_1.size(); i++) {
       std::cout << ID_LIST_1[i] << " ";
    }
    for(i=0; i < ID_LIST_2_mod.size(); i++) {
      ID_LIST_2_mod[i] = ID_LIST_2_mod[i] + MAX_ID_a;
      std::cout << ID_LIST_2_mod[i] << " ";
    }
    std::cout << std::endl << std::endl;

    for(i=0; i < file_list.size(); i++) {
      temp_path_x = path_a + file_list[i];
      temp_path_c = path_c + file_list[i];
      std::cout << "# Moving file " << temp_path_x << " to " << temp_path_c << std::endl  << std::endl;
      rename(temp_path_x.c_str(), temp_path_c.c_str());
    }

    std::ofstream out;
    std::ifstream input;

    for(i=0; i < id_n; i++) {
      temp_path_x = path_b + file_list[i];
      temp_path_c = path_c + file_list[i];
      std::cout << "# Merging files " << temp_path_x << " and " << temp_path_c << std::endl;
      
      TEMP_DATA.clear();
      RET_status = load_data_matrix(&TEMP_DATA,temp_path_x);

      for(j=0; j < ID_LIST_2.size(); j++) {
        for(k=0; k < TEMP_DATA.size(); k++) {
          if(TEMP_DATA[k][0] == ID_LIST_2[j]) {
            TEMP_DATA[k][0] = ID_LIST_2_mod[j];
          }
        }
      }

      RET_status = save_mat(temp_path_c,&TEMP_DATA);
    }

    for(i=id_n; i < (file_list.size()-1); i++) {
      temp_path_x = path_b + file_list[i];
      temp_path_c = path_c + file_list[i];
      std::cout << "# Merging files " << temp_path_x << " and " << temp_path_c << std::endl;
      
      TEMP_DATA.clear();
      RET_status = load_data_matrix(&TEMP_DATA,temp_path_x);
      RET_status = save_mat(temp_path_c,&TEMP_DATA);
    }

    
    temp_path_x = path_b + file_list.back();
    temp_path_c = path_c + file_list.back();
    std::cout << "# Merging files " << temp_path_x << " and " << temp_path_c << std::endl;

    std::string buffer;
    std::ifstream in_file;
    in_file.open(temp_path_x.c_str(), std::ios::in);

    std::ofstream out_file;
    out_file.open(temp_path_c.c_str(), std::ios::app);
    if(in_file.fail() || out_file.fail()) {
        return -1;
    }

    while(!in_file.eof()) {
      buffer.clear();
      getline(in_file,buffer);
      out_file << buffer;
    }
    in_file.close();
    out_file.close();

    std::cout << "# Merging done" << std::endl;
    remove_folder(path_a);
    std::cout << "Removed folder: " << path_b << std::endl;
    remove_folder(path_b);
    std::cout << "Removed folder: " << path_b << std::endl;

    return 0;
}

void remove_folder(const std::string& name) {
  namespace fs = boost::filesystem;
  fs::remove_all(name);
}

std::vector<std::string> list_files(const std::string& folder) {
  std::vector<std::string> ret;
  namespace fs = boost::filesystem;

  boost::progress_timer t( std::clog );

  fs::path full_path( fs::initial_path<fs::path>() );

  full_path = fs::system_complete( fs::path( folder ) );

  unsigned long file_count = 0;
  unsigned long dir_count = 0;
  unsigned long other_count = 0;
  unsigned long err_count = 0;

  if ( !fs::exists( full_path ) ) {
    std::cout << "\nNot found: " << full_path.string() << std::endl;
    return ret;
  }

  if ( fs::is_directory( full_path ) ) {
    std::cout << "\nIn directory: " << full_path.string() << "\n\n";
    fs::directory_iterator end_iter;
    for ( fs::directory_iterator dir_itr( full_path ); dir_itr != end_iter; ++dir_itr ) {
      try {
        if ( fs::is_directory( dir_itr->status() ) ) {
          ++dir_count;
          std::cout << dir_itr->path().filename() << " [directory]\n";
        }
        else if ( fs::is_regular_file( dir_itr->status() ) ) {
          ++file_count;
          ret.push_back(dir_itr->path().string());
          std::cout << dir_itr->path().filename() << "\n";
        }
        else {
          ++other_count;
          std::cout << dir_itr->path().filename() << " [other]\n";
        }

      }
      catch ( const std::exception & ex ) {
        ++err_count;
        std::cout << dir_itr->path().filename() << " " << ex.what() << std::endl;
      }
    }
    std::cout << "\n" << file_count << " files\n"
              << dir_count << " directories\n"
              << other_count << " others\n"
              << err_count << " errors\n";
  }
  else { // must be a file
    std::cout << "\nFound: " << full_path.string() << "\n";    
  }
  return ret;
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
        return -1;
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

    return 0;
}

int save_mat(std::string out_file_name, std::vector<std::vector<double> > *M) {

        unsigned int i;
        unsigned int j;

        std::ofstream out_pos(out_file_name.c_str(), std::ios::app);
        out_pos << std::setprecision(20) << std::scientific;
        if(out_pos.fail()) {
            return -1;
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

        return 0;
}
