#include <time.h>

#ifndef DEFINE_CONSTS_HH
#define DEFINE_CONSTS_HH

#define sol_L 4e26

#define CMS_version 0.9
#define PI 3.14159265
#define G 6.67300e-11
#define G_gauss 0.017202098950000
#define Msol 1.98855e30
#define mu_o 1.32712440018e11
#define c0 2.99792458e+08
#define AU 1.4960e+11
#define my0 4*3.14159265*1e-07

#endif

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

#include "ast2body.h"

#include "OVERLOADING.hh"
#include "TEMPLATES.hh"

#ifndef DEFINE_STRUCTS_HH
#define DEFINE_STRUCTS_HH

typedef std::vector<double>::const_iterator double_vector_iterator;

struct double_vec_ordering_struct {
    bool operator ()(std::pair<size_t, double_vector_iterator> const& a, std::pair<size_t, double_vector_iterator> const& b) {
        return *(a.second) < *(b.second);
    }
};


struct OPTIONS {
    unsigned int loops; //must be above 1e+02
    unsigned int precision;
    double memory_alloc;
    unsigned int output_type;
    unsigned int save_interval;
    unsigned int coord;
    unsigned int energy_diag;
    double sun_collision;
    unsigned int hill_sphere;
    double ejection_criteria;
    unsigned int integrator;
    double dt;

    unsigned int body_n;
    unsigned int massive_body_n;
    unsigned int memory_clear_interval;

    unsigned int load_old;

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

    unsigned int removed_bodies;
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



