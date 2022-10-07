/*
 * CLASS_FUNCTIONS.cc
 *
 *  Created on: Jul 2, 2015
 *      Author: dankas
 */
#include <vector>
#include <math.h>
#include <sstream>
 
#include "define.hh"
#include "functions.hh"
#include "class.hh"

void output_stream::endl() {
	switch(type) {
		case 0:
			out1_ << std::endl;
			break;
		case 1:
			out1_ << std::endl;
			out2_ << '\n';
			break;
		default:
			std::cout << "Could not read output type" << std::endl;
	}
}
