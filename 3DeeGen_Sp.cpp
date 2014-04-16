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

//TODO
// 1st knockout(Sp | k, angk, angNN)
// 2nd detector acceptance
// 3nd fram transform  
// add all dsc to get count

using namespace std;

int main(int argc, char *argv[]){
  time_t Tstart=time(0);

  if(argc < 9) {
    printf("===============Generating Flourine (p,pn) knockout data======================\n");
    printf("          Only for F(p,2p)O knockout [A(a,cd)b]\n");
    printf("Usage: ./3DeeGee_Sp.o MA Z JA JB theta_k theta_NN k dSp\n");
    printf("      MA = Mass number of isotop \n");
    printf("      Z  = charge number of isotop \n");
    printf("      JA = spin of Flourine isotop\n");
    printf("      JB = spin of Residual \n");
    printf(" theta_k = theta of resdiual \n");
    printf("theta_NN = 2-body CM scattering angle \n");    
    printf("      BE = binding energy of orbital proton \n");
    printf("      dk = step of k \n\n");
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
  float theta_k = atof(argv[5]);
  float theta_NN = atof(argv[6]);
  float BE = atof(argv[7]);
  int  kStep = atoi(argv[8]);
  
  int N = 0;
  int L = 0;
  float J = 0.5;

  int totCount = 6*300/kStep;

  float *jaco = new float[2]; // jacobian for x-sec transform
  float *output = new float[8]; // knockout output 

  char filename[50];
  sprintf(filename, "paraOut_%2d%s_JA%3.1f_JB%3.1f_angk%d_angNN%d_Sp%4.1f.dat",  MA, symbolZ(Z), JA, JB, (int)theta_k, (int)theta_NN, BE);
  
//#############################  display input condition
  printf("===========================\n");
  printf(" %d%s(p,2p)%d%s \n",MA, symbolZ(Z), MA-1, symbolZ(Z-1));
  printf("JA = %3.1f,  JB = %3.1f\n", JA, JB);
  printf(" theta_k = %6.2f deg\n", theta_k);
  printf("theta_NN = %6.2f deg\n", theta_NN);
  printf("Be = %10.4f MeV\n", BE);
  printf("k Step = %2d deg, total loops = %10d \n", kStep, totCount);
  printf(" output: %s \n", filename);  
  printf("------------------------------------------------------\n");

  //#######################rewrite output file
  FILE * paraOut;
  paraOut = fopen (filename, "w");
  // file header
  fprintf(paraOut, "#A(a,cd)B = %2dF(p,2p)%2dO, JA=%3.1f  JB=%3.1f\n", MA, MA-1, JA, JB);
  fprintf(paraOut, "#BE=%5.1f  Ti=%9.3f theta_k=%9.3f theta_NN=%9.3f\n", BE, Ti, theta_k, theta_NN);
  fprintf(paraOut, "# DWIA \n");
  fprintf(paraOut, "#%8d   ", 1);
  for (int i = 2; i <= 22 ; i ++){
    fprintf(paraOut, "%9d   ", i);
  }
  fprintf(paraOut, "\n");
  fprintf(paraOut, "%9s   %9s   %9s   %9s   %9s   %9s   %9s   %9s   %9s   ",
           "k", "T1", "theta1", "T2", "theta2", "T_c","theta_c", "T_d", "theta_d");
  for (int ID= 1; ID <=6 ; ID++){
    orbit(ID, N, L, J);
    char NLJ[9];
    sprintf(NLJ, "%1d%1s%1d/2", N, symbolL(L), (int)(2*J));
    fprintf(paraOut,"%9s   %9s   ", NLJ , NLJ);  
  }
  fprintf(paraOut, "\n");
    

  //########################### start looping
  int count = 0;
  for (float k=0; k <300; k+=kStep){

    // use knockout2D.h to calculate Tc thetac and thetad;
    output = Knockout2D(MA, Ti, k, theta_k, theta_NN, BE);

    // make_infile
    printf("%6d[%4.1f%%]| k:%6.2f | T_c:%9.3f, theta_c:%9.3f, theta_d:%9.3f\n",
             count, count*100./totCount, k, output[0], output[1], output[3]);


    // save parameters + readout
    fprintf(paraOut,"%9.3f   %9.3f   %9.3f   %9.3f   %9.3f   %9.3f   %9.3f   %9.3f   %9.3f   ",
             k, output[4],output[5],output[6],output[7], output[0], output[1], output[2], output[3]);       

    for (int ID=1; ID <=6; ID ++){

      //orbit
      orbit(ID, N, L, J);

      count ++;
      int checkmake = make_infile(MA, Z, JA, JB,  Ti, N, L, J, BE, output[0], output[1], output[3]);
      if (checkmake !=1) continue;

      // Calulate the thetaINV_c and Jacobian for X-sec
      // jaco = Jacobian(938.272, output[0], output[1], -0.6465223); // for 23F in SHARAQ04

      // run 3Dee code
      system("./threedee infile");
      
     
      // read_outfile : xsec + Ay
      read_outfile();
     
      // save parameters + readout
      //      fprintf(paraOut,"%9.6f   %9.6f   %9.6f   ",DWIA, PWIA, Ay ); 
      fprintf(paraOut,"%9.7f   %9.7f   ",DWIA ,Ay); 
      
      //Delete outfile
      remove("outfile");

    }

    fprintf(paraOut,"\n");
  }
  fprintf(paraOut,"##################################\n");

  fclose(paraOut);

  //########################## display result
  time_t Tend=time(0);
  printf("========== Totol run time %10.0f sec = %5.1f min| speed:#%5.2f/sec ===========\n",
         difftime(Tend,Tstart),difftime(Tend,Tstart)/60,count/difftime(Tend,Tstart)); 
  printf("  condition %2dF(p,2p)%2dO   nlJ=%1d%1s%3.1f  Ti:%7.2f MeV \n", MA, MA-1, N, symbolL(L), J, Ti );
  printf("  JA = %3.1f,  JB = %3.1f\n", JA, JB);
  printf("  assume Binding energy of orbital proton is %7.2f \n", BE );
  printf(" output: %s \n", filename);  
  // append infile to  output file
  char command[100];
  sprintf(command,"cat %s >> %s","infile", filename);
  system(command);
  
  return 0;

}
