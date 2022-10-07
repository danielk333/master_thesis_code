/*
 * MODELS.cc
 *
 *  Created on: Jul 9, 2015
 *      Author: dankas
 */

#include <math.h>

#include "define.hh"
#include "functions.hh"

double grun_cum_flux(double m) {
	std::vector<double> c;
	c.push_back(3.156e7);
	c.push_back(2.2e3);
	c.push_back(15);
	c.push_back(1.3e-9);
	c.push_back(1e11);
	c.push_back(1e27);
	c.push_back(1.3e-16);
	c.push_back(1e6);

	double F;
	F = c[0]*(pow(c[1]*pow(m,0.306) + c[2],-4.38) + c[3]*pow(m + c[4]*pow(m,2) + c[5]*pow(m,4),-0.36) + c[6]*pow(m + c[7]*pow(m,2),-0.85));

	return F;
}

std::vector<double> grun_flux(double min_mass_detectable, unsigned int grun_model_resolution) {
	// grun_cum_flux(m) is the cumulative flux of particles with mass m or larger
	// in particles/m^2/year to one side of a randomly tumbling plate which
	// is assumed as stationary with respect to the Earth surface.
	//
	//So grun_flux(min_mass_detectable,grun_model_resolution) = (grun_cum_flux(m_i) - grun_cum_flux(m_(i+1)))/((min_mass_detectable - 1e-18)/grun_model_resolution))

	std::vector<double> ret;
	unsigned int i;
	double model_start_val = 1e-18;
	double delta_m = (min_mass_detectable - model_start_val)/grun_model_resolution;

	for(i=0; i < grun_model_resolution; i++) {
		ret.push_back((grun_cum_flux(model_start_val + i*delta_m) - grun_cum_flux(model_start_val + (i+1)*delta_m))/delta_m);
	}

	return ret;
}


