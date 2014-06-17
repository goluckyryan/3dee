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
#include "constant.h"
#include "XsecTransform.h" //load Jacobian
#include "knockout2D.h"
#include "nuclei_mass.h"
#include "3DeeGenLibrary.h" 

using namespace std;

int main(int argc, char *argv[]){
  time_t Tstart=time(0);

  if(argc < 9) {
    printf("===============Generating Flourine (p,pn) knockout data======================\n");
    printf("          Only for F(p,2p)O knockout [A(a,cd)b]\n");
    printf("Usage: ./3DeeGee_k_angk.o MA Z JA JB BE dk dAngk dAngNN\n");
    printf("      MA = Mass number of isotop \n");
    printf("      Z  = charge number of isotop \n");
    printf("      JA = spin of Flourine isotop\n");
    printf("      JB = spin of Residual \n");
    printf("     dT1 = step of KE of proton 1 \n");
    printf("   dAng1 = step of ang1 of proton 1 \n");
    printf("     dT2 = step of KE of proton 2 \n");
    printf("   dAng2 = step of ang2 of proton 2 \n\n");
    exit(0);
  }

  //##################### variables
  int MA   = atoi(argv[1]);
  float Ti;
  switch (MA) {
  case 16: Ti = 200.0; break;
  case 23: Ti = 289.44; break;
  case 25: Ti = 276.779; break;  
  default: Ti = 200.0; break;
  }

  int Z = atoi(argv[2]);
  float JA = atof(argv[3]);
  float JB = atof(argv[4]);
  int   T1Step = atoi(argv[5]);
  int ang1Step = atoi(argv[6]);
  int   T2Step = atoi(argv[7]);
  int ang2Step = atoi(argv[8]);
  
  int N = 0;
  int L = 0;
  float J = 0.5;

  int totCount = 6*(300/T1Step+1)*(60/ang1Step+1)*(300/T2Step+1)*(60/ang2Step+1);

  float *jaco = new float[2]; // jacobian for x-sec transform
  float *output = new float[9]; // knockout output 

  char filename[50];
  sprintf(filename, "../result/paraOut_%2d%s_JA%3.1f_JB%3.1f.dat",  MA, symbolZ(Z), JA, JB);
  
//#############################  display input condition
  printf("===========================\n");
  printf(" %d%s(p,2p)%d%s \n",MA, symbolZ(Z), MA-1, symbolZ(Z-1));
  printf("Ti = %10.3f \n", Ti);
  printf("JA = %3.1f,  JB = %3.1f\n", JA, JB);
  printf("T1 step = %2d MeV, ang1 step = %2d, T2 step = %2d MeV, ang2 step = %2d, total loops = %10d \n", T1Step, ang1Step, T2Step, ang2Step, totCount);
  printf(" output: %s \n", filename);  
  printf("------------------------------------------------------\n");

  //#######################rewrite output file
  FILE * paraOut;
  paraOut = fopen (filename, "w");
  // file header
  fprintf(paraOut, "#A(a,cd)B = %2dF(p,2p)%2dO, JA=%3.1f  JB=%3.1f\n", MA, MA-1, JA, JB);
  fprintf(paraOut, "#Ti=%9.3f\n", Ti);
  fprintf(paraOut, "#%143s", ""); for (int ID = 1; ID<=6 ; ID++) fprintf(paraOut, "%12s%12s", "DWIA", "A00n0") ; fprintf(paraOut, "\n");
  fprintf(paraOut, "#%11d", 1); for (int i = 2; i <= 10+2*6 ; i ++) fprintf(paraOut, "%12d", i); fprintf(paraOut, "\n");
  fprintf(paraOut, "%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s", 
          "T1", "theta1", "T2", "theta2", "k", "thetak", "thetaNN", "Sp", "T_c","theta_c", "T_d", "theta_d");
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
  for (float T1=10; T1 <=310; T1+=T1Step){
    for (float ang1 = 15; ang1<=75; ang1+=ang1Step){
      for (float T2=10; T2 <=310; T2+=T2Step){
        for(float ang2 = 15; ang2<=75; ang2+=ang2Step){

          // use knockout2D.h to calculate Tc thetac and thetad;
          output = Knockout2Dinv(MA, Z,  Ti, T1, ang1, T2, ang2);

          if ( isnan(output[0])) {
            printf(" --------- skipped due to impossible kinematic. \n");
            count += 6;
            continue;
          }

          // Sp Filter
          if ( output[3] <=0 || output[3] >=100) {
            printf(" skipped due to out range of  Sp = %9.3f | T1:%9.3f ang1:%9.3f T2:%9.3f ang2:%8.3f \n", output[3], T1, ang1, T2, ang2);
            count+= 6;
            continue;
          }


          // make_infile
          printf("%6d[%4.1f%%]| T1:%9.3f ang1:%9.3f T2:%9.3f ang2:%8.3f| k:%7.3f, theta_k:%7.3f, theta_NN:%7.3f, Sp:%7.3f\n",
                 count, count*100./totCount, T1, ang1, T2, ang2, output[0], output[1], output[2], output[3]);

          printf("                                 Tc:%9.3f, theta_c:%9.3f,Td:%9.3f, theta_d:%9.3f",
                 output[4], output[5], output[6], output[7]);
          

          // Calulate the thetaINV_c and Jacobian for X-sec
          // jaco = Jacobian(938.272, output[0], output[1], -0.6465223); // for 23F in SHARAQ04
        
          jaco = Jacobian2(output[4], output[5], output[6], output[7], output[8]);
          printf("        beta:%4.2f jaco:%10.7f , ang1:%10.7f \n", output[8], jaco[0], 180-jaco[1]);
      
          // save parameters + readout
          fprintf(paraOut,"%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f",
                  T1, ang1, T2, ang2, output[0],output[1],output[2],output[3], output[4], output[5], output[6], output[7]);       

          for (int ID=1; ID <=6; ID ++){
            count ++;

            //orbit
            orbit(ID, N, L, J);

            // make infile
            make_infile(MA, Z, JA, JB,  Ti, N, L, J, output[3], output[4], output[5], output[7]);

            // run 3Dee code
            system("./threedee infile");
     
            // read_outfile : xsec + Ay
            if (read_outfile(57) == 10) continue;
          
            // save parameters + readout
            fprintf(paraOut,"%12.6f%12.6f",jaco[0]*DWIA ,A00n0); 
      
            //Delete outfile
            //        remove("outfile");
       
          }
          fprintf(paraOut, "\n");
        }
        //fprintf(paraOut, "\n");
      }
      fprintf(paraOut,"\n");
    }
    //fprintf(paraOut,"\n");
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
  printf("  T1 step = %2d MeV, ang1 step = %2d, T2 step = %2d MeV, ang2 step = %2d, total loops = %10d \n", T1Step, ang1Step, T2Step, ang2Step, totCount);
  printf("  output: %s \n", filename);  
  printf("------------------------------------------------------\n");

  
  return 0;

}
