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
    else {
    	/* handle error condition */
    	return "";
    }
}

int calculate_memory_management(OPTIONS *opt) {
	double loop_mem;
	double const_mem;

	loop_mem = (double)(((double)sizeof(double))*6*(*opt).body_n)/((double)(*opt).save_interval);
	const_mem = (double)(sizeof(double)*6*(*opt).body_n*3);
	if (const_mem > (*opt).memory_alloc*8388608) {
		return NOT_ENOUGH_MEMORY;
	}
	
	(*opt).memory_clear_interval = (unsigned int)floor((((double)(*opt).memory_alloc)*8388608 - const_mem)/loop_mem);

	std::cout << "Memory clear interval calculated to every " << (*opt).memory_clear_interval << "'th itteration." << std::endl;

	return SIM_OK;
}

int console_output(SIMULATION *sim, OPTIONS *opt) {
	unsigned int j;
    double unit;
    std::string unit_str;
	switch((*sim).simulation_output) {
		case -1:
			std::cout << "Welcome to the Celestial Mechanics Simulator (CMS) " << CMS_version << std::endl;
			std::cout << "Programmed by Daniel Kastinen (dankas-1@student.ltu.se)" << std::endl;
			std::cout << "Program article can be found at www.this_does_not_exist_yet.com" << std::endl << std::endl;
			break;
		case 0:
			std::cout << "Simulation starting with: " << std::endl;
		    std::cout << (*opt).massive_body_n << " Massive objects" << std::endl;
		    std::cout << (*opt).body_n - (*opt).massive_body_n << " Massless objects" << std::endl;

    		std::cout << "#############################################################" << std::endl << std::endl;
			break;
		case 1:
        	(*sim).done++;
            time(&((*sim).end_t));
            (*sim).t_elapsed = difftime((*sim).end_t,(*sim).start_t);
            (*sim).t_estimate = ((*sim).t_elapsed)/(double)(*sim).loops_performed*((double)(*opt).loops - (double)(*sim).loops_performed);
            std::cout << std::setprecision(5) << std::fixed;
            std::cout << (*sim).loops_performed << "     |";
            unit = 1.0;
            if((*sim).t_elapsed > 600) {
                unit = 60;
                unit_str = "m";
            }
            else if((*sim).t_elapsed/60 > 600) {
                unit = 3600;
                unit_str = "h";
            }
            else if((*sim).t_elapsed/3600 > 600) {
                unit = 3600*24;
                unit_str = "d";
            }
            std::cout << (*sim).t_elapsed/unit << " " << unit_str << "      |";
            if((*sim).t_estimate > 600) {
                unit = 60;
                unit_str = "m";
            }
            else if((*sim).t_estimate/60 > 600) {
                unit = 3600;
                unit_str = "h";
            }
            else if((*sim).t_estimate/3600 > 600) {
                unit = 3600*24;
                unit_str = "d";
            }
            std::cout << (*sim).t_estimate/unit << " " << unit_str << "             |";
            std::cout << ((*sim).loops_performed*(*opt).dt)/(3600*24*365.25) << " y" << std::endl;
			break;
		case 2:

	std::cout << "- |" << std::endl << "SIMULATION DONE" << std::endl;
    std::cout << "#############################################################" << std::endl << std::endl;

    std::cout << "Massive bodies simulated : " << (*opt).massive_body_n << "" << std::endl;
    std::cout << "Massless bodies simulated: " << (*opt).body_n - (*opt).massive_body_n << "" << std::endl;
    std::cout << "                         : " << (*sim).removed_bodies << " Objects ejected" << std::endl << std::endl;
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
