#ifndef DEFINE_CONSTS_HH
#define DEFINE_CONSTS_HH

#define TRUE 1
#define FALSE 0

#define GHOPE_version 0.1
#define PI 3.14159265358979323846264338327950288419716939937510
#define G 6.67300e-11
#define Msol 1.98855e30
#define c0 2.99792458e+08
#define AU 1.4960e+11
#define L0 3.846E+26
#define s2day 1.15740740740741e-05
#define s2year 3.16880878140289e-08
#define day2s 86400.0
#define year2s 31557600.0

#endif

//*********************************************************

#include "DATA_TYPE.hh"
#include "OVERLOADING.hh"

#include <cstddef>
#include <time.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <vector>
#include <unistd.h>
#include <limits.h>
#include <limits>
#include <fstream>
#include <algorithm>
#include <typeinfo>
#include <cfloat>
#include <sys/types.h>
#include <ostream>
#include <algorithm>
#include <cmath>


#ifndef BS_DEF_HH
#define BS_DEF_HH
#include "BS/nr3.h"
#include "BS/stepper.h"
#include "BS/stepperbs.h"
#include "BS/odeint.h"
#endif

#include "AST/ast2body.h"

//#include "OVERLOADING.hh"
//#include "TEMPLATES.hh"

//*********************************************************

#ifndef OUT_OVERLOAD
#define OUT_OVERLOAD

//template <typename T> std::string Num2Str ( T Number );
template <typename T>
    std::string Num2Str ( T Number ) {
        std::ostringstream string_stream;
        string_stream << Number;
        return string_stream.str();
    }

class output_stream {
    public:
    std::ostream& out1_;
    std::ostream& out2_;
    unsigned int type;

    output_stream(std::ostream& out1) : out1_(out1),out2_(out1) {}
    output_stream(std::ostream& out1,std::ostream& out2) : out1_(out1), out2_(out2) {}

    void endl();
};


template <typename T>
output_stream& operator<<(output_stream& O, T const& S) {
    switch(O.type) {
        case 0:
            O.out1_ << S << std::endl;
            break;
        case 1:
            O.out1_ << S << std::endl;
            O.out2_ << S << '\n';
            break;
        default:
            std::cout << "Could not read output type" << std::endl;
    }

    return O;
}

#endif

//*********************************************************

#ifndef DEFINE_STRUCTS_HH
#define DEFINE_STRUCTS_HH

typedef std::vector<double>::const_iterator double_vector_iterator;

struct double_vec_ordering_struct {
    bool operator ()(std::pair<size_t, double_vector_iterator> const& a, std::pair<size_t, double_vector_iterator> const& b) {
        return *(a.second) < *(b.second);
    }
};


// Options for the simulation
struct OPTIONS {
    //Loaded by the user
    unsigned int steps;
    double dt;
    double start_time;
    double end_time;
    double sun_collision_limit;
    double close_enc_limit;
    double ejection_limit;
    bool load_old;
    unsigned int central_body;
    bool PR_effect;

    unsigned int integrator;
    double ATOL;
    double RTOL;
    
    unsigned int save_interval_itter;
    double save_inteval_time;
    bool save_close;
    bool save_remove;
    bool save_log;
    bool save_summary;
    bool snapshot_list;
    unsigned int coord;
    bool save_energy_diag;
    bool save_time_diag;
    unsigned int precision;
    
    double max_RAM;
    int delim_int;
    unsigned int verbose_level;
};

// Simulation statistics and performance data
struct SIMULATION {
	unsigned int loops_performed;
    unsigned int steps_performed;
    unsigned int nok;
    unsigned int nbad;

    unsigned int body_n;
    unsigned int interacting_body_n;
    unsigned int noninteracting_body_n;
    unsigned int central_body_index;

	unsigned int save_number;
    unsigned int snapshot_counter;
    unsigned int save_buildup;
    unsigned int current_save_buildup;
    unsigned int save_counter;
    unsigned int save_interval;
    unsigned int last_align;
    unsigned int console_print_c;
    
    double t_estimate;
    double t_elapsed;
    double runtime;
    clock_t start_t, end_t;

    unsigned int removed_bodies;
};

#endif

//*********************************************************


