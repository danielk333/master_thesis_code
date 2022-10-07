#include <time.h>

#ifndef DEFINE_CONSTS_HH
#define DEFINE_CONSTS_HH

#define PI 3.14159265
#define G 6.67300e-11
#define c0 2.99792458e+08
#define AU 1.4960e+11
#define my0 4*3.14159265*1e-07

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


#endif



