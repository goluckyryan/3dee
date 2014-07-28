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
//#include "XsecTransform.h" //load Jacobian
#include "knockout2D.h"
#include "nuclei_mass.h"
#include "3DeeGenLibrary.h" 

using namespace std;

int main(int argc, char *argv[]){
  time_t Tstart=time(0);

  if(argc < 5) {
    printf("===============Generating Flourine (p,pn) knockout data======================\n");
    printf("          Only for F(p,2p)O knockout [A(a,cd)b]\n");
    printf("Usage: ./3DeeGee_Tc_angc_Td.o MA Z JA JB BE dTc dAngc dAngd\n");
    printf("      MA = Mass number of isotop \n");
    printf("      Z  = charge number of isotop \n");
    printf("      BE = seperation energy \n");
    printf("     dTc = step of KE of proton 1 \n");
    printf("    dAng = step of ang_c/d of proton 1/2 \n");
    exit(0);
  }

  //##################### variables
  int MA   = atoi(argv[1]);
  float Ti;
  switch (MA) {
  case 14: Ti = 245.46; break;
  case 16: Ti = 200.0; break;
  case 23: Ti = 289.44; break;
  case 25: Ti = 276.779; break;  
  default: Ti = 200.0; break;
  }

  int Z = atoi(argv[2]);
  float Sp = atof(argv[3]);
  int  TcStep = atoi(argv[4]);
  int angStep = atoi(argv[5]);
  
  int N = 0;
  int L = 0;
  float J = 0.5;

  float JA = 1.0;
  float JB = 1.5;

  const int orbRange = 6;
  const int TcStart = 10;
  const int TcEnd   = 300;
  const int angcRange = 180;
  const int angdRange = 180;


  int totCount = 0;

  for (float Tc=TcStart; Tc <=TcEnd+TcStart; Tc+=TcStep){
    for (float angc = 0; angc<=angcRange; angc+=angStep){
        for(float angd = 0; angd<=angdRange; angd+=angStep){
          totCount += orbRange;
        }
    }
  }

  float *output = new float[9]; // knockout output 

  char filename[50];
  sprintf(filename, "../result/paraOut_%2d%s_Sp%04.1f_Tc%03d_ang%03d.dat",  MA, symbolZ(Z), Sp,TcStep, angStep);
  
//#############################  display input condition
  printf("===========================\n");
  printf(" %d%s(p,2p)%d%s \n",MA, symbolZ(Z), MA-1, symbolZ(Z-1));
  printf("Sp = %6.3f, Ti = %10.3f \n", Sp, Ti);
  printf("JA = %3.1f,  JB = %3.1f\n", JA, JB);
  printf("Tc step = %2d MeV, angc step = %2d,  angd step = %2d, total loops = %10d \n", TcStep, angStep, angStep, totCount);
  printf(" output: %s \n", filename);  
  printf("------------------------------------------------------\n");

  //#######################rewrite output file
  FILE * paraOut;
  paraOut = fopen (filename, "w");
  // file header
  fprintf(paraOut, "##A(a,cd)B = %2dF(p,2p)%2dO, JA=%3.1f  JB=%3.1f\n", MA, MA-1, JA, JB);
  fprintf(paraOut, "##Sp=%6.3f, Ti=%9.3f\n", Sp, Ti);
  fprintf(paraOut, "#X%15s%2d MeV, Range (%6.1f, %6.1f)\n", "Tc step =", TcStep, TcStart, TcEnd+TcStart); 
  fprintf(paraOut, "#Y%15s%2d deg, Range (%6.1f, %6.1f)\n", "angc step =", angStep, 0, angcRange); 
  fprintf(paraOut, "#Z%15s%2d deg, Range (%6.1f, %6.1f)\n", "angd step =", angStep, 0, angdRange); 
  fprintf(paraOut, "#%*s", 12*11,""); for (int ID = 1; ID<=orbRange ; ID++) fprintf(paraOut, "%12s%12s", "DWIA", "A00n0") ; fprintf(paraOut, "\n");
  fprintf(paraOut, "#%12d", 1); for (int i = 2; i <= 10+2*orbRange ; i ++) fprintf(paraOut, "%12d", i); fprintf(paraOut, "\n");
  fprintf(paraOut, "#%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s", 
          "Tc", "theta_c", "Td", "theta_d", "k", "thetak", "thetaNN", "T_1","theta_1", "T_2", "theta_2");
  for (int ID= 1; ID <=orbRange ; ID++){
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
  int effCount = 0;
  for (float Tc=TcStart; Tc <=TcEnd+TcStart; Tc+=TcStep){
    for (float angc = 0; angc<=angcRange; angc+=angStep){
        for(float angd = 0; angd<=angdRange; angd+=angStep){

          // use knockout2D.h to calculate Tc thetac and thetad;
          output = Knockout2Dinv3(MA, Z,  Ti, Tc, angc, -angd, Sp);

          if ( isnan(output[0])) {
            //printf(" --------- skipped due to impossible kinematic. \n");
            count += orbRange;
            continue;
          }

          // Accpetance Filter
          if (  bool jjj = !AccpetanceFilter2D(output[4], output[5], output[6], output[7])) {
            //printf(" Accpatance Filter, %2d|(%3.1f,%3.1f,%3.1f), T1:%9.3f ang1:%9.3f T2:%9.3f ang2:%8.3f \n",jjj, Tc, angc, -angd, output[4], output[5], output[6], output[7]);
            count+= orbRange;
            continue;
          }

          // print condition
          printf("\e[32m==== %6d[\e[31m%4.1f%%\e[32m]| Tc:%9.3f angc:%9.3f Td:%9.3f angd:%8.3f| k:%7.3f, theta_k:%7.3f, theta_NN:%7.3f \e[m\n",
                 count, count*100./totCount, Tc, angc, output[3], -angd, output[0], output[1], output[2]);

          printf("                                 T1:%9.3f, theta_1:%9.3f,T2:%9.3f, theta_2:%9.3f\n",
                 output[4], output[5], output[6], output[7]);
         
          // save parameters + readout
          fprintf(paraOut," %12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f",
                  Tc, angc, output[3], -angd, output[0],output[1],output[2], output[4], output[5], output[6], output[7]);       

          for (int ID=1; ID <=orbRange; ID ++){
            count ++;
            effCount ++;
            //orbit
            orbit(ID, N, L, J);

            // make infile
            make_infile(MA, Z, JA, JB,  Ti, N, L, J, Sp, Tc, angc, -angd, 0);

            // run 3Dee code
            system("./threedee infile");
     
            // read_outfile : xsec + Ay
            if (read_outfile(57) == 10) continue;
          
            // save parameters + readout
            fprintf(paraOut,"%12.6f%12.6f",DWIA ,A00n0); 

            if ( DWIA > 1) printf("\e[31m");
            printf("                                                                DWIA(%1d%1s%1d/2):%14.10f \n", N, symbolL(L), (int)(2*J), DWIA);
            if ( DWIA > 1) printf("\e[m");
            //Delete outfile
            //        remove("outfile");
       
          }
          fprintf(paraOut, "\n");
      }
        //fprintf(paraOut,"\n");
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
  printf("\e[32m[%5.1f%%]========== Totol run time %10.0f sec = %5.1f min = %5.1f hr| speed:#%5.2f(%5.2f)/sec ===========\e[m\n",
         count*100./totCount,difftime(Tend,Tstart),difftime(Tend,Tstart)/60,difftime(Tend,Tstart)/3600,count/difftime(Tend,Tstart),effCount/difftime(Tend,Tstart)); 
  printf("  condition %2d%s(p,2p)%2d%s   Ti:%7.2f MeV \n", MA, symbolZ(Z) ,MA-1,symbolZ(Z-1),  Ti );
  printf("  JA = %3.1f, JB = %3.1f\n, Sp", JA, JB, Sp);
  printf("  Tc step = %2d MeV, ang step = %2d, total loops = %10d \n", TcStep, angStep, totCount);
  printf("  output: %s \n", filename);  
  printf("------------------------------------------------------\n");

  
  return 0;

}
