#include <time.h>

#ifndef DEFINE_CONSTS_HH
#define DEFINE_CONSTS_HH

#define OrbClone_version 0.1
#define PI 3.14159265
#define G 6.67300e-11
#define Msol 1.9891e30
#define c0 2.99792458e+08
#define AU 1.4960e+11


#endif


#ifndef DEFINE_SIMVARS_HH
#define DEFINE_SIMVARS_HH

enum SIMULATION_STATUS {
	CALC_OK						=0,
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



