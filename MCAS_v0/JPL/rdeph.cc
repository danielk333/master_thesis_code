/**==========================================================================**/
/**                                                                          **/
/**  SOURCE FILE: rdeph.c                                                    **/
/**                                                                          **/
/**      Purpose: This program interpolates a state from data in a JPL       **/
/**               ephemeris file for a given Julian date. It is designed     **/
/**               to be called from a Tcl/Tk script, so it has a very        **/
/**               unfriendly user interface.                                 **/
/**                                                                          **/
/**               For convenience, here is the mapping between the command   **/
/**               line inputs and the global variables found in read.tcl:    **/
/**                                                                          **/
/**                        argv[1]  <-->  EphName                            **/
/**                        argv[2]  <-->  JD                                 **/
/**                        argv[3]  <-->  TN                                 **/
/**                        argv[4]  <-->  FmtOpt                             **/
/**                        argv[5]  <-->  Cnt                                **/
/**                                                                          **/
/**   Programmer: David Hoffman/EG5                                          **/
/**               NASA, Johnson Space Center                                 **/
/**               Houston, TX 77058                                          **/
/**               e-mail: david.a.hoffman1@jsc.nasa.gov                      **/
/**                                                                          **/
/**==========================================================================**/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ephem_read.hh"

#include "functions.hh"

#ifndef TYPES_DEFINED
#include "ephem_types.hh"
#endif
                
int main (int argc, char *argv[]) {

  stateType  State;
  char       ephemFileName[12] , *tgtName;
  double     Position[3] , Time;
  int        i , Target;

  /* Convert command line arguments to numeric values.........................*/
  
  Time   = atof(argv[2]);
  Target = atoi(argv[3]);

  /* Initialize the ephemeris.................................................*/

  Initialize_Ephemeris(argv[1]);

  /* Compute the desired ephemeris data.......................................*/

  if ( !strcmp(argv[4],"PosVel") ) 
      Interpolate_State( Time , Target , &State );
  else
      Interpolate_Position( Time , Target , Position );

  /* Print the answer.........................................................*/

  if ( !strcmp(argv[4],"PosVel") ) {
       /* Print seperator.....................................................*/
       
       for ( i=0 ; i<21 ; i++ ) putchar('-');
       printf("  Case %2s  ",argv[5]);
       for ( i=0 ; i<21 ; i++ ) putchar('-');
       
       /* Print planet name...................................................*/
       
       switch (Target) {
          case  0: tgtName = "Mercury";
                   break;
          case  1: tgtName = "Venus";
                   break;
          case  2: tgtName = "Earth";
                   break;
          case  3: tgtName = "Mars";
                   break;
          case  4: tgtName = "Jupiter";
                   break;
          case  5: tgtName = "Saturn";
                   break;
          case  6: tgtName = "Uranus";
                   break;
          case  7: tgtName = "Neptune";
                   break;
          case  8: tgtName = "Pluto";
                   break;
          case  9: tgtName = "Moon";
                   break;
          case 10: tgtName = "Sun";
                   break;
          default: tgtName = "Mercury";
                   break;
       }
       
       printf("\n\n  Target:  %s",tgtName);
       printf("\n      JD:  %8.2f",Time);

       /* Print position......................................................*/
       
       printf("\n\n  Position (km):     [1] =  % 22.15e",State.Position[0]);
       printf("\n                     [2] =  % 22.15e",State.Position[1]);
       printf("\n                     [3] =  % 22.15e",State.Position[2]);

       /* Print velocity......................................................*/

       printf("\n\n  Velocity (km/sec): [1] =  % 22.15e",State.Velocity[0]);
       printf("\n                     [2] =  % 22.15e",State.Velocity[1]);
       printf("\n                     [3] =  % 22.15e\n\n\n",State.Velocity[2]);

       std::vector<double> temp_vec;
       std::vector<std::vector<double> > temp_mat;

       std::string output_folder;
       output_folder = argv[6];

       std::string output_file;

       output_file = output_folder + "pos.data";
       temp_vec.push_back((State.Position[0])*1e3/AU);
       temp_vec.push_back((State.Position[1])*1e3/AU);
       temp_vec.push_back((State.Position[2])*1e3/AU);
       temp_mat.push_back(temp_vec);
       save_mat(output_file.c_str(), &temp_mat);

       temp_vec.clear();

       output_file = output_folder + "vel.data";
       temp_vec.push_back((State.Velocity[0])*1e3);
       temp_vec.push_back((State.Velocity[1])*1e3);
       temp_vec.push_back((State.Velocity[2])*1e3);
       temp_mat[0] = temp_vec;
       save_mat(output_file.c_str(), &temp_mat);
      }
  else {
       /* Print seperator.....................................................*/
       
       for ( i=0 ; i<21 ; i++ ) putchar('-');
       printf("  Case %2s  ",argv[5]);
       for ( i=0 ; i<21 ; i++ ) putchar('-');

       /* Print planet name...................................................*/
       
       switch (Target) {
          case  1: tgtName = "Mercury";
                   break;
          case  2: tgtName = "Venus";
                   break;
          case  3: tgtName = "Earth";
                   break;
          case  4: tgtName = "Mars";
                   break;
          case  5: tgtName = "Jupiter";
                   break;
          case  6: tgtName = "Saturn";
                   break;
          case  7: tgtName = "Uranus";
                   break;
          case  8: tgtName = "Neptune";
                   break;
          case  9: tgtName = "Pluto";
                   break;
          case 10: tgtName = "Moon";
                   break;
          case 11: tgtName = "Sun";
                   break;
          default: tgtName = "Mercury";
                   break;
       }
       
       printf("\n\n  Target:  %s",tgtName);
       printf("\n      JD:  %8.2f",Time);

       /* Print position......................................................*/
       
       printf("\n\n  Position (km):     [1] =  % 22.15e",Position[0]);
       printf("\n                     [2] =  % 22.15e",Position[1]);
       printf("\n                     [3] =  % 22.15e\n\n\n",Position[2]);
     }

  /* Exit normally............................................................*/
  
  exit(0);
}
