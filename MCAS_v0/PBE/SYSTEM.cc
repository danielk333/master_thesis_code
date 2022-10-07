#include <math.h>
#include <iostream>
#include <sys/stat.h>

#include "define.hh"
#include "functions.hh"


bool file_exists(const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

std::string get_selfpath() {
    char buff[PATH_MAX];
    ssize_t len = ::readlink("/proc/self/exe", buff, sizeof(buff)-1);
    if (len != -1) {
      buff[len] = '\0';
      return std::string(buff);
    }
    /* handle error condition */
}

int console_output(SIMULATION *sim, OPTIONS *opt) {
	unsigned int j;
	switch((*sim).simulation_output) {
		case -1:
            std::cout << std::endl << std::endl;
			std::cout << "---- Parent Body Ejector (ṔBE) program running " << PBE_version << " ----" << std::endl;
			std::cout << "Author: Daniel Kastinen (dankas-1@student.ltu.se)" << std::endl;
			break;
		case 0:
			std::cout << "Simulation starting with: " << std::endl;
		    std::cout << "                          " << (*opt).massive_body_n << " Massive objects" << std::endl;

    		std::cout << "#############################################################" << std::endl << std::endl;
			break;
		case 1:
        	(*sim).done++;
            time(&((*sim).end_t));
            (*sim).t_elapsed = difftime((*sim).end_t,(*sim).start_t);
            (*sim).t_estimate = ((*sim).t_elapsed)*((double)1-(double)(*sim).loops_performed/(double)(*opt).loops);
            std::cout << std::endl << "##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~## " << std::endl;
            std::cout << "Loop: " << (*sim).loops_performed << std::endl;
            std::cout << "Elapsed simulation time: " << (*sim).t_elapsed/3600 << " h" << std::endl;
            std::cout << "Estimated time left: " << (*sim).t_estimate << " s" << std::endl;
            std::cout << "                   : " << (*sim).t_estimate/60 << " min" << std::endl;
            std::cout << "                   : " << (*sim).t_estimate/3600 << " h" << std::endl << std::endl;
            std::cout << "Simulating | ------------------------- |" << std::endl << "           | ";
            for(j = 0; j < (unsigned int)((*sim).done/4); j++) {
              	std::cout << "-";
            }
			break;
		case 2:

	std::cout << "- |" << std::endl << "SIMULATION DONE" << std::endl;
    std::cout << "#############################################################" << std::endl << std::endl;
    
    std::cout << "Loops performed          : " << (*sim).loops_performed << std::endl;
    std::cout << "Simulation time          : " << (*sim).t_elapsed/3600 << " h" << std::endl;
    std::cout << "Time/itteration          : " << (*sim).t_elapsed/(*opt).loops << " s/loop" << std::endl;
    std::cout << "Total time simulated     : " << (*opt).loops*(*opt).dt << " s" << std::endl;
    std::cout << "                         : " << (*opt).loops*(*opt).dt/(3600*24*365) << " earth years" << std::endl << std::endl;

			break;
		default:
		return CONSOLE_OUTPUT_FAILED;
	}

	return SIM_OK;
}
