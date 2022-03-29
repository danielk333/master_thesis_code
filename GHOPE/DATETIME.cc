#include <vector>
#include <math.h>

#include "resources.hh"
#include "functions.hh"


double JDtoJ(double JD) {
	return 2000.0 + (JD - 2451545.0)/365.25;
}