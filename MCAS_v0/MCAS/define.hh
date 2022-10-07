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

#include "OVERLOADING.hh"
#include "TEMPLATES.hh"

#ifndef DEFINE_CONSTS_HH
#define DEFINE_CONSTS_HH

#define MCAS_version 1.0
#define PI 3.14159265
#define G 6.67300e-11
#define G_gauss 0.017202098950000
#define mu_o 1.32712440018e11
#define Msol 1.98855e30
#define Mearth 5.97219e24
#define c0 2.99792458e+08
#define AU 1.4960e+11
#define my0 4*3.14159265*1e-07

#define Rearth 6371e3
#define Hatm 120e3

#define TRUE 1
#define FALSE 0

#define MAX_THREADS 5

#endif

#ifndef DEFINE_STRUCTS_HH
#define DEFINE_STRUCTS_HH

typedef std::vector<double>::const_iterator double_vector_iterator;

struct double_vec_ordering_struct {
    bool operator ()(std::pair<size_t, double_vector_iterator> const& a, std::pair<size_t, double_vector_iterator> const& b) {
        return *(a.second) < *(b.second);
    }
};

struct OPTIONS {
    unsigned int MODE;
    unsigned int integrator;
    unsigned int IC_gen;
    double start_date;
    double integration_time;
    double end_date;
    unsigned int histogram_n;
    unsigned int family;
    unsigned int dist_type;
    bool calc_clusters;
    bool stream_calc;
    unsigned int mass_dist;
    unsigned int mass_dist_res;
    double mass_min;
    double mass_max;
    double mem_allowed;
    double close_match;
    double snapshot_date;
    unsigned int mass_format;
    unsigned int logfile;
    unsigned int X_type;
    unsigned int X;
    unsigned int merc_int;
    double merc_step;
    unsigned int merc_pr;
    double merc_eject;
    double merc_close_int;
    unsigned int mercury6_max;
    double SUOC_sigma;
    unsigned int SUOC_hist_res;
};


#endif


#ifndef DEFINE_SIMVARS_HH
#define DEFINE_SIMVARS_HH

enum DATA_TYPES {
  DATA_VECTOR  			     		=1,
  DATA_MATRIX  			     		,
  DATA_MATRIX_SILENT        ,
  DATA_FAMILY  			     		,
  DATA_MERCURY  			     	,
  DATA_OPTIONS  			     	,
  DATA_ORBDIST						,
  DATA_ORBDIST_SS					,
};

enum CALCULATION_STATUS {
  CALC_OK      			     		=0,
  CONFIG_LOAD_FAILED         		=-1,
  FAILD_TO_DETECT_DATA_LOAD_TYPE 	=-3,
  FILE_DOES_NOT_EXIST				=-4,
  DISTRIBUTION_NOT_RECOGNIZED		=-5,
  mercury6_big_DATA_SAVE_FAILED		=-6,
  mercury6_small_DATA_SAVE_FAILED =-7,
  NO_INIT_EXIST              		=-9,
  VECTOR_LOAD_FAILED         		=-10,
  MATRIX_LOAD_FAILED   	     		=-11,
  CONSOLE_OUTPUT_FAILED      		=-19,
  DATA_SAVE_FAILED           		=-20,
  DATA_READ_FAILED           		=-21,
  mercury6_MATRIX_LOAD_FAILED		=-30,
  NOT_SUPPORTED_MASS_DIST       =-31
};

void clean_up(std::vector<std::vector<std::string> > matrix_created_files, int type);
void error(int SIM_STATUS, std::vector<std::vector<std::string> > matrix_created_files);

#endif



