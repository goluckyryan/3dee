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
#include "knockout3D.h"
#include "nuclei_mass.h"
#include "3DeeGenLibrary.h" 

using namespace std;

int main(int argc, char *argv[]){
  time_t Tstart=time(0);

  if(argc < 7) {
    printf("===============Generating Nucleus (p,pn) knockout data======================\n");
    printf("          Only for A(p,2p)O knockout [A(a,cd)b]\n");
    printf("Usage: ./3DeeGee_Tc_angc_Td.o MA Z BE dTc dAng dPhi\n");
    printf("      MA = Mass number of isotop \n");
    printf("      Z  = charge number of isotop \n");
    printf("      BE = seperation energy \n");
    printf("     dTc = step of KE of proton 1 \n");
    printf("    dAng = step of ang of protons \n");
    printf("    dPhi = step of phi of protons \n");
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
  int   TcStep = atoi(argv[4]);
  int  angStep = atoi(argv[5]);
  int  phiStep = atoi(argv[6]);
  
  int N = 0;
  int L = 0;
  float J = 0.5;

  float JA = 1.0;
  float JB = 1.5;

  bool runTHREEDEE = 1;

  const int orbRange = 6;
  const int TcStart = 10;
  const int TcRange = 300;
  const int angcRange = 180;
  const int angdRange = 180;
  const int phicRange = 25; //+-
  const int phidRange = 25; //+-

  int totCount = 0;

  for (float Tc=TcStart; Tc <=TcRange+TcStart; Tc+=TcStep){
    for (float angc = 15; angc<=angcRange; angc+=angStep){
      for(float angd = 15; angd<=angdRange; angd+=angStep){
        for(float phic = -phicRange; phic<=phicRange; phic+=phiStep){
          for(float phid = -phidRange-180; phid<=phidRange-180; phid+=phiStep){
            totCount += orbRange;

            //printf("Tc:%5.1f, theta_c:%5.1f, phi_c:%5.1f, theta_d:%5.1f, phi_d:%5.1f |", Tc, angc, phic, angd, phid);

            //printf(" totcount = %d\n", totCount);
          }
        }
      }
    }
  }

  float *output = new float[13]; // knockout output 

  char filename[50];
  sprintf(filename, "../result/paraOut_%2d%s_Sp%04.1f_Tc%03d_ang%03d_phi%03d.dat",  MA, symbolZ(Z), Sp,TcStep, angStep, phiStep);
  
  //#############################  display input condition
  printf("===========================\n");
  printf(" %d%s(p,2p)%d%s \n",MA, symbolZ(Z), MA-1, symbolZ(Z-1));
  printf("Sp = %6.3f, Ti = %10.3f \n", Sp, Ti);
  printf("JA = %3.1f,  JB = %3.1f\n", JA, JB);
  printf("Tc step = %2d MeV, angc step = %2d,  phic step = %2d\n", TcStep, angStep, phiStep);
  printf("angd step = %2d, phid Step = %2d, total loops = %10d \n", angStep, phiStep, totCount);
  printf(" output: %s \n", filename);  
  printf("------------------------------------------------------\n");

  //#######################rewrite output file
  FILE * paraOut;
  paraOut = fopen (filename, "w");
  // file header
  fprintf(paraOut, "##A(a,cd)B = %2dF(p,2p)%2dO, JA=%3.1f  JB=%3.1f\n", MA, MA-1, JA, JB);
  fprintf(paraOut, "##Sp=%6.3f, Ti=%9.3f\n", Sp, Ti);
  fprintf(paraOut, "#X%15s%2d MeV, Range (%6d, %6d)\n", "Tc step =", TcStep, TcStart, TcRange); 
  fprintf(paraOut, "#Y%15s%2d deg, Range (%6d, %6d)\n", "angc step =", angStep, 0, angcRange); 
  fprintf(paraOut, "#Z%15s%2d deg, Range (%6d, %6d)\n", "angd step =", angStep, 0, angdRange);
  fprintf(paraOut, "#B%15s%2d deg, Range (%6d, %6d)\n", "phic step =", phiStep, -phicRange, phicRange);
  fprintf(paraOut, "#A%15s%2d deg, Range (%6d, %6d)\n", "phid step =", phiStep, -phidRange, phidRange);
  int paraNum = 18;
  fprintf(paraOut, "#%*s", 12*paraNum," cross section unit = ub"); for (int ID = 1; ID<=orbRange ; ID++) fprintf(paraOut, "%12s%12s", "DWIA", "A00n0") ; fprintf(paraOut, "\n");
  fprintf(paraOut, "#"); for (int i = 1; i <= paraNum+2*orbRange ; i ++) fprintf(paraOut, "%12d", i); fprintf(paraOut, "\n");
  fprintf(paraOut, "#%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s", 
          "Tc", "theta_c", "phi_c","Td", "theta_d", "phi_d", "beta_d", "k", "theta_k", "phi_k", "theta_NN", "phi_NN", "T_1","theta_1", "phi_1", "T_2", "theta_2", "phi_2");
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
  for (float Tc=TcStart; Tc <=TcRange+TcStart; Tc+=TcStep){
    for (float angc = 15; angc<=angcRange; angc+=angStep){
      for(float angd = 15; angd<=angdRange; angd+=angStep){
        for(float phic = -phicRange; phic<=phicRange; phic+=phiStep){
          for(float phid = -phidRange-180; phid<=phidRange-180; phid+=phiStep){

            //refine phid
            float phid2 = 0;
            if(phid > 180) {
              phid2 = 360-phid; 
            }else if(phid < -180){
              phid2 = 360+phid;
            }else{
              phid2 = phid;
            }

            // use knockout2D.h to calculate Tc thetac and thetad;
            output = Knockout3Dinv3(MA, Z,  Ti, Tc, angc, phic, angd, phid2, Sp);

            float betad = 0;
            if(phid2 >= 0){
              betad = 180 + phic - phid2;
            }else{
              betad = -180 + phic - phid2;
            }
            
            // Td is nan
            if ( isnan(output[7])) {
              //printf(" --------- skipped due to impossible kinematic. \n");
              count += orbRange;
              continue;
            }

            // Accpetance Filter
            if (  bool jjj = !AccpetanceFilter3D(output[0], output[1], output[2], output[3], output[4], output[5])) {
              //printf(" Accpatance Filter, %2d|(%3.1f,%3.1f,%4.1f,%3.1f%4.1f), T1:%9.3f ang1:%9.3f phi1:%9.3f, T2:%9.3f ang2:%8.3f phi2:%8.3f \n",jjj, Tc, angc, phic,angd, phid2, output[0], output[1], output[2], output[3], output[4], output[5]);
              count+= orbRange;
              continue;
            }

            // print condition
            printf("\e[32m==== %6d(%6d)[\e[31m%4.1f%%\e[32m]| Tc:%5.1f angc:%5.1f phic:%5.1f Td:%5.1f angd:%5.1f phid:%6.1f| betad:%5.1f, k:%7.3f, theta_k:%7.3f, phi_k:%7.3f, theta_NN:%7.3f, phi_NN:%7.3f\e[m\n",
                   count, effCount, count*100./totCount, Tc, angc, phic, output[7], angd, phid2,  betad,  output[8], output[9], output[10], output[11], output[12]);

            printf("        T1:%9.3f, theta_1:%9.3f, phi_1:%9.3f, T2:%9.3f, theta_2:%9.3f, phi_2:%9.3f\n",
                   output[0], output[1], output[2], output[3], output[4], output[5]);
         
            // save parameters + readout
            fprintf(paraOut," %12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f",
                    Tc, angc, phic, output[7], angd, phid2,  betad, output[8], output[9], output[10], output[11], output[12],output[0], output[1], output[2], output[3], output[4], output[5]);       

            for (int ID=1; ID <=orbRange; ID ++){
              count ++;
              effCount ++;

              if ( runTHREEDEE){              
                //orbit
                orbit(ID, N, L, J);

                // make infile
                make_infile(MA, Z, JA, JB,  Ti, N, L, J, Sp, Tc, angc, -angd, betad);

                // run 3Dee code
                system("./threedee infile");
     
                // read_outfile : xsec + Ay
                if (read_outfile(57) == 10) continue;
          
                // save parameters + readout
                fprintf(paraOut,"%12.6f%12.6f",DWIA*1000.,A00n0); 

                if ( DWIA > 1) printf("\e[31m");
                printf("                                                                DWIA(%1d%1s%1d/2):%14.10f ub\n", N, symbolL(L), (int)(2*J), DWIA*1000.);
                if ( DWIA > 1) printf("\e[m");
                //Delete outfile
                //        remove("outfile");
              }
       
            }
            fprintf(paraOut, "\n");
          }
        }
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
  printf("\e[32m[%5.1f%%](%5.1f%%)========== Totol run time %10.0f sec = %5.1f min = %5.1f hr| speed:#%5.2f(%5.2f)/sec ===========\e[m\n",
         count*100./totCount,effCount*100./totCount,difftime(Tend,Tstart),difftime(Tend,Tstart)/60,difftime(Tend,Tstart)/3600,count/difftime(Tend,Tstart),effCount/difftime(Tend,Tstart)); 
  printf("  condition %2d%s(p,2p)%2d%s   Ti:%7.2f MeV \n", MA, symbolZ(Z) ,MA-1,symbolZ(Z-1),  Ti );
  printf("  JA = %3.1f, JB = %3.1f,  Sp =%7.2f", JA, JB, Sp);
  printf("Tc step = %2d MeV, angc step = %2d,  phic step = %2d\n", TcStep, angStep, phiStep);
  printf("angd step = %2d, phid Step = %2d, total loops = %10d, effLoop = %10d \n", angStep, phiStep, totCount, effCount);
  printf("  output: %s \n", filename);  
  printf("------------------------------------------------------\n");

  
  return 0;

}
