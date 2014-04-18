#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <stdio.h>
#include <math.h>
#include <cstring>
#include "XsecTransform.h" //load Jacobian
#include "knockout2D.h"
#include "3DeeGenLibrary.h" 

using namespace std;


int main(int argc, char *argv[]){
  time_t Tstart=time(0);

  if(argc < 8) {
    printf("===============Generating Flourine (p,pn) knockout data======================\n");
    printf("          Only for F(p,2p)O knockout [A(a,cd)b]\n");
    printf("Usage: ./3DeeGee_angk_angNN.o MA Z JA JB k BE dAng\n");
    printf("      MA = Mass number of isotop \n");
    printf("      Z  = charge number of isotop \n");
    printf("      JA = spin of isotop\n");
    printf("      JB = spin of Residual \n");
    printf("       k = Fermi-momentum of orbital proton \n");   
    printf("      BE = binding energy of orbital proton \n");
    printf("    dAng = step of angle \n\n");
    exit(0);
  }

  //##################### variables
  int MA   = atoi(argv[1]);
  float Ti;
  switch (MA) {
  case 16: Ti = 200.0; break;
  case 23: Ti = 289.44; break;
  case 25: Ti = 276.779; break;  
  }
  int Z = atoi(argv[2]);
  float JA = atof(argv[3]);
  float JB = atof(argv[4]);
  float k = atof(argv[5]);
  float BE = atof(argv[6]);
  int  AngStep = atoi(argv[7]);
  
  int N = 0;
  int L = 0;
  float J = 0.5;

  int totCount = 6*360*360/AngStep/AngStep;

  float *jaco = new float[2]; // jacobian for x-sec transform
  float *output = new float[8]; // knockout output 

  char filename[50];
  sprintf(filename, "paraOut_%2d%s_JA%3.1f_JB%3.1f_k%d_Sp%4.1f.dat",  MA, symbolZ(Z), JA, JB, (int)k, BE);
  
//#############################  display input condition
  printf("===========================\n");
  printf(" %d%s(p,2p)%d%s \n",MA, symbolZ(Z), MA-1, symbolZ(Z-1));
  printf("JA = %3.1f,  JB = %3.1f\n", JA, JB);
  printf(" k = %6.2f MeV/c\n", k);
  printf("Be = %10.4f MeV\n", BE);
  printf("Ang Step = %2d deg, total loops = %10d \n", AngStep, totCount);
  printf(" output: %s \n", filename);  
  printf("------------------------------------------------------\n");

  //#######################rewrite output file
  FILE * paraOut;
  paraOut = fopen (filename, "w");
  // file header
  fprintf(paraOut, "#A(a,cd)B = %2dF(p,2p)%2dO, JA=%3.1f  JB=%3.1f\n", MA, MA-1, JA, JB);
  fprintf(paraOut, "#BE=%5.1f  Ti=%9.3f k=%9.3f\n", BE, Ti, k);
  fprintf(paraOut, "#%119s", ""); 
  for (int ID = 1; ID<=6 ; ID++) fprintf(paraOut, "%12s%12s", "DWIA", "A00n0") ; fprintf(paraOut, "\n");
  fprintf(paraOut, "#%11d", 1); for (int i = 2; i <= 10+2*6 ; i ++) fprintf(paraOut, "%12d", i); fprintf(paraOut, "\n");
  fprintf(paraOut, "%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s", 
            "angk", "angNN", "T1", "theta1", "T2", "theta2", "T_c","theta_c", "T_d", "theta_d");
  for (int ID= 1; ID <=6 ; ID++){
    orbit(ID, N, L, J);
    char NLJ[9];
    sprintf(NLJ, "%1d%1s%1d/2", N, symbolL(L), (int)(2*J));
    for (int i = 1; i <= 2; i++){
      fprintf(paraOut,"%12s", NLJ);  
    }
  }
  fprintf(paraOut, "\n");
    

  //########################### start looping
  int count = 0;
  for (float angk=0; angk <360; angk+=AngStep){
    for(float angNN=0; angNN <=360; angNN+=AngStep){

      // use knockout2D.h to calculate Tc thetac and thetad;
      output = Knockout2D(MA, Z, Ti, k, angk, angNN, BE);

      // make_infile
      printf("%6d[%4.1f%%]| angk:%6.2f angNN:%6.2f| T_c:%9.3f, theta_c:%9.3f, theta_d:%9.3f\n",
             count, count*100./totCount, angk, angNN, output[0], output[1], output[3]);
      
      if ( isnan(output[0])) {
        printf(" skipped due to impossible kinematic. \n");
        continue;
      }

      // save parameters + readout
      fprintf(paraOut,"%12.3f%9.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f",
              angk, angNN, output[4],output[5],output[6],output[7], output[0], output[1], output[2], output[3]);       

      for (int ID=1; ID <=6; ID ++){
        count ++;
        
        //orbit
        orbit(ID, N, L, J);

        // make infile
        make_infile(MA, Z, JA, JB,  Ti, N, L, J, BE, output[0], output[1], output[3]);

        // Calulate the thetaINV_c and Jacobian for X-sec
        // jaco = Jacobian(938.272, output[0], output[1], -0.6465223); // for 23F in SHARAQ04

        // run 3Dee code
        system("./threedee infile");
      
        // read_outfile : DWIA, PWIA,  Ay
        read_outfile(58);
     
        // save parameters + readout 
        fprintf(paraOut,"%12.6f%12.6f",DWIA ,A00n0); 
      
        //Delete outfile
        //remove("outfile");

      }
      fprintf(paraOut,"\n");
    }
    fprintf(paraOut,"\n");
  }
  fprintf(paraOut,"##################################\n");

  fclose(paraOut);
  // append infile to  output file
  char command[100];
  sprintf(command,"cat %s >> %s","infile", filename);
  system(command);

  //########################## display result
  time_t Tend=time(0);
  printf("========== Totol run time %10.0f sec = %5.1f min| speed:#%5.2f/sec ===========\n",
         difftime(Tend,Tstart),difftime(Tend,Tstart)/60,count/difftime(Tend,Tstart)); 
  printf("  condition %2d%s(p,2p)%2d%s   Ti:%7.2f MeV \n", MA, symbolZ(Z) ,MA-1,symbolZ(Z-1),  Ti );
  printf("  JA = %3.1f,  JB = %3.1f\n", JA, JB);
  printf("  assume Binding energy of orbital proton is %7.2f \n", BE );
  printf("  output: %s \n", filename);  

  
  return 0;

}
