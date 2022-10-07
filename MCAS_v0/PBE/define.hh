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
#include <typeinfo>
#include <cfloat>
#include <sys/types.h>
#include <ostream>


#ifndef DEFINE_CONSTS_HH
#define DEFINE_CONSTS_HH

#define PBE_version 0.1
#define PI 3.14159265
#define Ggrav 6.67300e-11
#define mu_o 1.32712440018e11
#define c0 2.99792458e+08
#define AU 1.4960e+11
#define my0 4*3.14159265*1e-07
#define k_B 1.38064852e-23

#define water_W 580
#define sol_L 4e26
#define body_H 2e6
#define sol_R 696342e3

#define T_g0 300
#define mu_Wip 20*1.661e-24*1e-3
#define tau_Wip 0.25
#define body_H_W 1.88e6
#define K_drag 26.0/9.0

#define g_Hug 0.1
#define xi_Hug 2.0
#define theta_Hug 2.0
#define T_g_Hug 230
#define mu_Hug 3.0e-26
#define tau_Hug 0.038
#define Psi_Hug 0.4

#define TRUE 1
#define FALSE 0

#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include "boost/progress.hpp"
//#include "boost/array.hpp"
//#include "boost/operators.hpp"
//#include "boost/numeric/odeint.hpp"


//#include "boost_point_type.hpp"

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
    unsigned int ejection_model;
    unsigned int loops;
    unsigned int loops_orbit;
    unsigned int particles;
    unsigned int orbits;
    unsigned int precision;
    unsigned int output_type;
    double memory_alloc;
    unsigned int massive_body_n;
    unsigned int ejection_number;
    double model_max_vel;

    double dt;
    unsigned int body_n;

};

struct SIMULATION {
	unsigned int done;
	unsigned int loops_performed;
	int simulation_output;

	unsigned int save_number;
	unsigned int next_clear;
    unsigned int diag_counter;
    unsigned int save_counter;
    
    double t_estimate;
    double t_elapsed;
    time_t start_t, end_t;
};

#endif


#ifndef DEFINE_SIMVARS_HH
#define DEFINE_SIMVARS_HH

enum SIMULATION_STATUS {
	SIM_OK						=0,
	CONFIG_LOAD_FAILED		    =-1,
	SCHW_CONFIG_LOAD_FAILED     =-2,
    NO_INIT_EXIST               =-9,
	VECTOR_LOAD_FAILED          =-10,
	MATRIX_LOAD_FAILED			=-11,
	NOT_ENOUGH_MEMORY			=-15,
    FAILED_INTEGRATOR_SETTINGS  =-16,
	CONSOLE_OUTPUT_FAILED		=-19,
	DATA_SAVE_FAILED        	=-20,
    DIAG_SAVE_FAILED            =-21,
    SETTINGS_SAVE_FAILED        =-22
};

void error(int err);

#endif



