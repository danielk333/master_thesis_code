#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <iostream>

#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#define AU 1.4960e+11

#include "SpiceUsr.h"

/*
   Define the maximum length for any string, 80
   characters plus one null terminator.
*/
#define STRLEN 81

int save_mat(const char *out, std::vector<std::vector<double> > *M) {

  unsigned int i;
  unsigned int j;
  std::string out_file_name = out;

  std::ofstream out_pos(out_file_name.c_str(), std::ios::app);

  if(out_pos.fail()) {
    return -1;
  }

  for(j = 0; j < (*M).size(); j++) {
    for(i = 0; i < (*M)[j].size(); i++) {
      out_pos << (((*M)[j])[i]);
      if(i < ((*M)[j].size() - 1)) {
        out_pos << " ";
      }
      else {
        out_pos << "\r\n";
      }
    }
  }

  out_pos.close();

  return 0;
}


int main( int argc, char **argv ) {
  unsigned int i,inp;
  if(argc < 6) {
    printf ("NOT ENOUGH INPUT ARGUMENTS");
    exit(-10);
  }


  /*  Declare the needed variables: */

  double     Time;

  /*
  Set a flag to start/stop and continue the
  inquiry loop.
  */

  SpiceDouble            state[6];
  SpiceDouble            ltime;

  SpiceDouble       et;
  SpiceChar    pictur[37] = "Wkd Mon DD HR:MN:SC YYYY ::UTC";
  SpiceInt          lenout;
  SpiceChar         outstring[37];

  std::vector<double> temp_vec;
  std::vector<std::vector<double> > temp_mat;
  std::string output_file;
  std::string output_folder;
  std::string targ;
  std::string type;
  std::string opt;

  std::vector<std::string> OPTS;
  OPTS.push_back("Selfpath      :");
  OPTS.push_back("Target        :");
  OPTS.push_back("Type          :");
  OPTS.push_back("Time          :");
  OPTS.push_back("Output folder :");
  OPTS.push_back("Options       :");
  for(i = 6; i < argc; i++) {
    OPTS.push_back("KERNEL FILE   :");
  }
  /*
      The RETURN mode signals an error then returns to the
      caller. Just what we need. REPORT mode performs almost
      the same function as RETURN, however RETURN mode
      sets the return_c() value to TRUE and so the program does
      not execute those CSPICE routines that check the return_c()
      value. Consider REPORT mode useful for debugging.
  */
  erract_c ( "SET", STRLEN, "RETURN" );
  for(i = 0; i < argc; i++) {
    std::cout << OPTS[i] << " " << argv[i] << std::endl;
  }

  inp=1;
  targ = argv[inp];inp++;
  type = argv[inp];inp++;

  if(type.compare("BARYCENTER") == 0) {
    targ = targ + " " + type;
  }

  Time = atof(argv[inp]);inp++;
  output_folder = argv[inp];inp++;
  opt = argv[inp];inp++;
  et = Time;

  /*
  Load the data we need for state evaluation.
  */

  for(i = inp; i < argc; i++) {
    furnsh_c ( argv[i] );
  }

  spkezr_c ( targ.c_str(), Time, "ECLIPJ2000", opt.c_str(), "SUN", state, &ltime );

  if ( ! failed_c() ) {
    timout_c ( et, pictur, 37, outstring );
    printf ( "%s\n", outstring);
    printf ( "R : %17.5f %17.5f %17.5f\n", 
            state[0] , state[1], state[2] );
    printf ( "V : %17.5f %17.5f %17.5f\n",
            state[3] , state[4], state[5] );
    printf ( "LT: %f\n", ltime );

    output_file = output_folder + "pos.data";
    temp_vec.push_back((state[0])*1e3/AU);
    temp_vec.push_back((state[1])*1e3/AU);
    temp_vec.push_back((state[2])*1e3/AU);
    temp_mat.push_back(temp_vec);
    save_mat(output_file.c_str(), &temp_mat);

    temp_vec.clear();
    output_file = output_folder + "vel.data";
    temp_vec.push_back((state[3])*1e3);
    temp_vec.push_back((state[4])*1e3);
    temp_vec.push_back((state[5])*1e3);
    temp_mat[0] = temp_vec;
    save_mat(output_file.c_str(), &temp_mat);
  }
  else {
    /*
    Problem. Something went wrong. Reset the error
    subsystem for another pass.
    */
    reset_c();
  }

  return 0;
}
