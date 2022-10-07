/*
 * INTEGRATOR_SETTINGS.cc
 *
 *  Created on: May 20, 2015
 *      Author: dankas
 */

#include <vector>

#include "define.hh"
#include "functions.hh"


int load_integrator_scheme(std::vector<std::vector<double> > *integrator_coef, OPTIONS *opt) {

	std::vector<double> integrator_coef_temp;

		integrator_coef_temp.push_back(0.5);

		(*integrator_coef).push_back(integrator_coef_temp);
		integrator_coef_temp.clear();

		integrator_coef_temp.push_back(1);

		(*integrator_coef).push_back(integrator_coef_temp);


		/*integrator_coef_temp.push_back(0.03809449742241219545697532230863756534060);
		integrator_coef_temp.push_back(0.1452987161169137492940200726606637497442);
		integrator_coef_temp.push_back(0.2076276957255412507162056113249882065158);
		integrator_coef_temp.push_back(0.4359097036515261592231548624010651844006);
		integrator_coef_temp.push_back(-0.6538612258327867093807117373907094120024);

		(*integrator_coef).push_back(integrator_coef_temp);
		integrator_coef_temp.clear();

		integrator_coef_temp.push_back(0.09585888083707521061077150377145884776921);
		integrator_coef_temp.push_back(0.2044461531429987806805077839164344779763);
		integrator_coef_temp.push_back(0.2170703479789911017143385924306336714532);
		integrator_coef_temp.push_back(0.01737538195906509300561788011852699719871);
		(*integrator_coef).push_back(integrator_coef_temp);*/
	

	return SIM_OK;
}
