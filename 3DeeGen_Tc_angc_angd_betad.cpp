#include <istream> 
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
 
  if(argc != 9) { 
    printf("===============Generating Nucleus (p,2p) knockout data======================\n"); 
    printf("\e[33m  (p,2p) for A(a,cd)B, A = B + b knockout in inverse kinematics \e[m\n");
    printf("\e[32m  using infile.2p.temp \e[m\n"); 
    printf("Usage: ./3DeeGee_Tc_angc_Td.o MA Z Ti BE dTc dAng dPhi comment\n"); 
    printf("      MA : Mass number of isotop \n"); 
    printf("      Z  : charge number of isotop \n");
    printf("      Ti : incident KE of a \n"); 
    printf("      BE : seperation energy \n"); 
    printf("     dTc : step of KE of a \n"); 
    printf("    dAng : step of ang of protons \n"); 
    printf("    dPhi : step of phi of protons \n"); 
    printf(" comment : comment \n"); 
    exit(0); 
  } 
 
  //##################### variables
  string temp_file = "infile.2p.EDAD1_NL_nLS.temp";
   
  int MA   = atoi(argv[1]); 
  int Z = atoi(argv[2]); 
  float Ti = atof(argv[3]);
  
  switch (MA) { 
  case 23: Ti = 289.44; break; 
  case 25: Ti = 276.779; break;   
  }

  float Sp = atof(argv[4]); 
  int   TcStep = atoi(argv[5]); 
  int  angStep = atoi(argv[6]); 
  int  phiStep = atoi(argv[7]); 
  char* comment = argv[8];
   
  int N = 1; 
  int L = 0; 
  float J = 0.5; 
 
  float JA = 1.0; 
  float JB = 1.5; 
 
  bool runTHREEDEE = 1;

  const int xIA = 2; // 1 = PWIA, 2 = DWIA

  const int orbStart = 4;
  const int orbEnd   = 4; 
  float TcStart   = 10; 
  float TcEnd     = Ti-Sp;
  float angcStart = 40; 
  float angcEnd   = 40;
  float angdStart = 40; 
  float angdEnd   = 40; 
  float phicStart = -0; // local angle
  float phicEnd   =  0;  
  float phidStart = -0; // local angle in 3D, it will change to global in Kinematics cal.
  float phidEnd   =  0;  
 
  int totCount = 0; 
 
  for (float Tc=TcStart; Tc <=TcEnd; Tc+=TcStep){ 
    for (float angc = angcStart; angc<=angcEnd; angc+=angStep){
      for(float angd = angdStart; angd<=angdEnd; angd+=angStep){ 
        for(float phic = phicStart; phic<=phicEnd; phic+=phiStep){ 
          for(float phid = -phidStart-180; phid<=phidEnd-180; phid+=phiStep){ 
            totCount += orbEnd - orbStart + 1; 
 
            //printf("Tc:%5.1f, theta_c:%5.1f, phi_c:%5.1f, theta_d:%5.1f, phi_d:%5.1f |", Tc, angc, phic, angd, phid); 
 
            //printf(" totcount = %d\n", totCount); 
          } 
        } 
      } 
    } 
  } 
 
  float *output = new float[13]; // knockout output  
 
  char filename[200]; 
  sprintf(filename, "../result/3d_%2d%s_Ti%04.0f_Sp%04.1f_Tc%03d_ang%03d_phi%03d_%s.dat",  MA, symbolZ(Z, MA), Ti, Sp,TcStep, angStep, phiStep, comment); 
   
  //#############################  display input condition 
  printf("===========================\n"); 
  printf(" %d%s(p,2p)%d%s \n",MA, symbolZ(Z, MA), MA-1, symbolZ(Z-1, MA)); 
  printf(" Sp = %6.3f, Ti = %10.3f \n", Sp, Ti); 
  printf(" JA = %3.1f,  JB = %3.1f\n", JA, JB);
  printf("Tc step  :%4d MeV, Range (%6.1f, %6.1f)\n", TcStep, TcStart, TcEnd); 
  printf("angc step:%4d deg, Range (%6.1f, %6.1f)\n", angStep, angcStart, angcEnd); 
  printf("angd step:%4d deg, Range (%6.1f, %6.1f)\n", angStep, angdStart, angdEnd);  
  printf("phic step:%4d deg, Range (%6.1f, %6.1f)\n", phiStep, phicStart, phicEnd); 
  printf("phid step:%4d deg, Range (%6.1f, %6.1f)\n", phiStep, phidStart, phidEnd);  
  printf(" orb start = %2d, orb End = %2d\n", orbStart, orbEnd); 
  printf(" total loops = %10d \n", totCount); 
  printf(" using template : %s \n", temp_file.c_str());
  printf(" output: %s \n", filename);   
  printf("------------------------------------------------------\n"); 
 
  //####################### rewrite output file 
  FILE * paraOut; 
  paraOut = fopen (filename, "w+"); 

  // file header 
  fprintf(paraOut, "##A(a,cd)B = %2dF(p,2p)%2dO, JA=%3.1f  JB=%3.1f\n", MA, MA-1, JA, JB); 
  fprintf(paraOut, "##Sp=%6.3f, Ti=%9.3f\n", Sp, Ti); 
  fprintf(paraOut, "#X%15s%2d MeV, Range (%6.1f, %6.1f)\n", "Tc step =", TcStep, TcStart, TcEnd);  
  fprintf(paraOut, "#Y%15s%2d deg, Range (%6.1f, %6.1f)\n", "angc step =", angStep, angcStart, angcEnd);  
  fprintf(paraOut, "#Z%15s%2d deg, Range (%6.1f, %6.1f)\n", "angd step =", angStep, angdStart, angdEnd); 
  fprintf(paraOut, "#B%15s%2d deg, Range (%6.1f, %6.1f)\n", "phic step =", phiStep, phicStart, phicEnd); 
  fprintf(paraOut, "#A%15s%2d deg, Range (%6.1f, %6.1f)\n", "phid step =", phiStep, phidStart, phidEnd); 
  int paraNum = 18; 
  fprintf(paraOut, "#%*s", 12*paraNum," cross section unit = ub"); 
  for (int ID = orbStart; ID<=orbEnd ; ID++) {
    if( xIA == 1) fprintf(paraOut, "%12s%12s", "PWIA", "A00n0") ;
    if( xIA == 2) fprintf(paraOut, "%12s%12s", "DWIA", "A00n0") ; 
  }
  fprintf(paraOut, "\n"); 
  fprintf(paraOut, "#"); 
  for (int i = 1; i <= paraNum+2*(orbEnd - orbStart + 1) ; i ++) fprintf(paraOut, "%12d", i); fprintf(paraOut, "\n"); 
  fprintf(paraOut, "#%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s",  
          "Tc", "theta_c", "phi_c","Td", "theta_d", "phi_d", "beta_d", "k", "theta_k", "phi_k", "theta_NN", "phi_NN", "T_1","theta_1", "phi_1", "T_2", "theta_2", "phi_2"); 
  for (int ID= orbStart; ID <=orbEnd ; ID++){ 
    orbit(ID, N, L, J); 
    char NLJ[9]; 
    sprintf(NLJ, "%1d%1s%1d/2", N, symbolL(L), (int)(2*J)); 
    for (int i = 1; i <= 2; i++){ 
      fprintf(paraOut,"%12s", NLJ);   
    } 
  } 
  fprintf(paraOut, "\n"); 
  fflush(paraOut); 
  //fclose(paraOut); 
 
  //########################### start looping 
  int count = 0; 
  int effCount = 0; 
  for (float Tc=TcStart; Tc <=TcEnd; Tc+=TcStep){ 
    for (float angc = angcStart; angc<=angcEnd; angc+=angStep){ 
      for(float angd = angdStart; angd<=angdEnd; angd+=angStep){ 
        for(float phic = phicStart; phic<=phicEnd; phic+=phiStep){ 
          for(float phid = phidStart-180; phid<=phidEnd-180; phid+=phiStep){ 
 
            //if( count % 600000 == 0 ) printf("count : %15d | Tc:%7.2f, angc:%7.2f, angd:%7.2f, phic:%7.2f, phid:%7.2f \n", count, Tc, angc, angd, phic, phid); 
 
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
            int kineticCal  = Knockout3Dinv3(MA, Z,  Ti, Tc, angc, phic, angd, phid2, Sp, output); 
 
            float betad = 0; 
            if(phid2 >= 0){ 
              betad = 180 + phic - phid2; 
            }else{ 
              betad = -180 + phic - phid2; 
            } 
             
            // Td is nan 
             
            if ( kineticCal == 0 ||  isnan(output[7])) { 
              //printf(" --------- skipped due to impossible kinematic. \n"); 
              count += orbEnd - orbStart + 1; 
              continue; 
            } 
 
            // Accpetance Filter
            /*
            bool jjj = !AccpetanceFilter3D(output[0], output[1], output[2], output[3], output[4], output[5]); 
            if (jjj) { 
              //printf(" Accpatance Filter, %2d|(%3.1f,%3.1f,%4.1f,%3.1f%4.1f), T1:%9.3f ang1:%9.3f phi1:%9.3f, T2:%9.3f ang2:%8.3f phi2:%8.3f \n",jjj, Tc, angc, phic,angd, phid2, output[0], output[1], output[2], output[3], output[4], output[5]); 
              count+= orbEnd; 
              continue; 
            }/**/ 
 
            // print condition 
            printf("\e[32m==== %d(%d)[\e[31m%4.1f%%\e[32m]| Tc:%5.1f angc:%5.1f phic:%5.1f Td:%5.1f angd:%5.1f phid:%6.1f| betad:%5.1f, k:%7.3f, theta_k:%7.3f, phi_k:%7.3f, theta_NN:%7.3f, phi_NN:%7.3f\e[m\n", 
                   count, effCount, count*100./totCount, Tc, angc, phic, output[7], angd, phid2,  betad,  output[8], output[9], output[10], output[11], output[12]); 
 
            printf("        T1:%9.3f, theta_1:%9.3f, phi_1:%9.3f, T2:%9.3f, theta_2:%9.3f, phi_2:%9.3f\n", 
                   output[0], output[1], output[2], output[3], output[4], output[5]); 
          
            // save parameters + readout 
            fprintf(paraOut," %12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f", 
                    Tc, angc, phic, output[7], angd, phid2,  betad, output[8], output[9], output[10], output[11], output[12],output[0], output[1], output[2], output[3], output[4], output[5]);        
 
            for (int ID=orbStart; ID <=orbEnd; ID ++){ 
              count ++; 
              effCount ++; 
 
              if ( runTHREEDEE){               
                //set N, L, J according to orbit ID 
                orbit(ID, N, L, J); 
                // make infile 
                make_infile(temp_file,MA, Z, JA, JB,  Ti, N, L, J, Sp, Tc, angc, -angd, betad); 
 
                // run 3Dee code 
                system("./threedee infile"); 
      
                // read_outfile : xsec + Ay 
                // 57 for Dirac 
                // 81 for Ogata 
                //if (read_outfile(81) == 10){ 
                if (read_outfile(57) == 10) { 
                  printf("READ outfile ERROR ==========> Please adjust the read Line, or revisit infile\n"); 
                  exit(1); 
                  continue; 
                } 
           
                // save parameters + readout 
                if( xIA == 1) fprintf(paraOut,"%12.6f%12.6f",PWIA*1000.,A00n0);  
                if( xIA == 2) fprintf(paraOut,"%12.6f%12.6f",DWIA*1000.,A00n0);  
               
                //printf("\e[31m%30sDWIA(%1d%1s%1d/2):%14.10f ub\e[m\n", "", N, symbolL(L), (int)(2*J), DWIA*1000.); 
                
              } 
        
            } 
            fprintf(paraOut, "\n"); 
          } 
          fflush(paraOut); 
        } 
      } 
    } 
  } 
  fprintf(paraOut,"##################################\n"); 
  fclose(paraOut); 
 
  // append infile to  output file 
  char command[100]; 
  sprintf(command,"cat %s >> %s","infile", filename); 
  system(command); 
 
  //########################## display result 
  time_t Tend=time(0); 
  printf("\e[34m==== %d%s(p,2p)%d%s @ %.3f MeV, Sp = %.3f MeV =====\e[m\n",MA, symbolZ(Z, MA), MA-1, symbolZ(Z-1, MA), Ti, Sp); 
  printf("tot count %d, count %d[%.1f%%], eff count %d[%.1f%%])\n",totCount,  count, count*100./totCount,effCount, effCount*100./totCount);
  printf("Totol run time %.0f sec = %.1f min = %.1f hr| speed:#%.2f(%.2f)/sec\n", difftime(Tend,Tstart),difftime(Tend,Tstart)/60,difftime(Tend,Tstart)/3600,count/difftime(Tend,Tstart),effCount/difftime(Tend,Tstart));    
  printf(" JA = %3.1f,  JB = %3.1f\n", JA, JB);
  printf(" orb start = %2d, orb End = %2d\n", orbStart, orbEnd); 
  printf("  Tc step:%4d MeV, Range (%6.1f, %6.1f)\n", TcStep, TcStart, TcEnd); 
  printf("angc step:%4d deg, Range (%6.1f, %6.1f)\n", angStep, angcStart, angcEnd); 
  printf("angd step:%4d deg, Range (%6.1f, %6.1f)\n", angStep, angdStart, angdEnd);  
  printf("phic step:%4d deg, Range (%6.1f, %6.1f)\n", phiStep, phicStart, phicEnd); 
  printf("phid step:%4d deg, Range (%6.1f, %6.1f)\n", phiStep, phidStart, phidEnd);  
  printf(" output: %s \n", filename);   
  printf(" %s \n", comment);
  printf(" using template : %s \n", temp_file.c_str());
  printf("------------------------------------------------------\n"); 
   
  return 0; 
 
} 
