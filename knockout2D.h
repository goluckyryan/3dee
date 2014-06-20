/*  Knockout of T(p,2p)R

The Lorentz Tranform on aribtary direction is well discrpted in Goldstien P .281

k(vector): momentum of a proton in Target nucleus
Theta_k:   direction of the momentum of a proton in O nucleus
Theta_NN:  scattering angle in 2-nucleon center of mass system
S:         separation energys

***** suffix notation ****
suffix "L" is Lab frame
suffix "" is nucleus frame
suffix "c" is CM frame
suffix "X" is experiment

***** notation ****
k is for momentum
P is for 4-momentum
E is Total relativistic energy
T is for Kinetic energy

***************** Procedure
we can transform the mulit-bpdy to 2-body scattering

in nucleus frame
A(a, 1 + 2) Resi = > Neo(a, 1) 2, such that Neo = A - Resi

by Lorentz transform of center of momentum 4-vector of New and a

then set the Theta_NN and Phi_NN in CM frame by rotation operator.

transform back to nucleus frame and lab frame

done

***************** Reverse problem

given Pi, P1, P2, find k, S, Theta_NN, Phi

Pi={Ei,0,0,ki}
P1 = {E1, k1*sin(Theta_P1), 0, k1*cos(Theta_P1)}
P2 = {E2, k2*sin(Theta_P2)*cos(Phi_P2), k2*sin(Theta_P2)*sin(Phi_P2), k2*cos(Theta_P2)}

S = sqrt((Ei+Mo-E1-E2)^2-(ki-k1-k2)^2)-Mo+mn

k = ki-k1-k2
***********************************/

#include <sstream> // enable string control
#include <cmath> // math , enable simple math function
#include <iomanip>
#include <stdio.h>  // standard I/O
#include <stdlib.h> // standard General Utlilities Library
#include "constant.h"
#include "nuclei_mass.h"

using namespace std; // declare a namespace "std", every variable in this code is inside "std"

//function___________________________________________________________________
float* Lorentz2D(float *V, float beta, float theta){
  float *U = new float[3];
    
  float gamma = 1/sqrt(1-beta*beta);
    
  U[0] = V[0]*gamma + V[1]*beta*gamma*sin(theta) + V[2]*beta*gamma*cos(theta);
  U[1] = V[0]*beta*gamma*sin(theta) + V[1]*(cos(theta)*cos(theta) + gamma*sin(theta)*sin(theta)) 
    + V[2]*(cos(theta)*sin(theta)*gamma-cos(theta)*sin(theta));
  U[2] = V[0]*beta*gamma*cos(theta) + V[1]*(cos(theta)*sin(theta)*gamma-cos(theta)*sin(theta))
    + V[2]*(sin(theta)*sin(theta) + gamma*cos(theta)*cos(theta));
    
  return U;
}


void PrintV2D(float *V, string Vname){
  printf("%10s : %10.3f, %10.3f, %10.3f\n", Vname.c_str(), V[0], V[1], V[2]); 
}

float Momentum2D(float *V){
  return sqrt(V[1]*V[1] + V[2]*V[2]);
}

float Mass2D(float *V){
  return sqrt(V[0]*V[0] - V[1]*V[1] - V[2]*V[2]);
}

float Angle(float *V){
  return atan2(V[1], V[2]);
}

float* Knockout2D(int MA, int Z, float TKEA, float k, float theta_k, float theta_NN, float Sp){
    
  float mass = Nucleus_Mass(Z, MA);

  if ( mass == -404) return(0);
  
  theta_k = theta_k/rad2deg;
  theta_NN = theta_NN/rad2deg;

  //  ########################### Nuclei's frame
  float *Pi = new float[3];
  float *Pt = new float[3];
  float *P1 = new float[3];
  float *P2 = new float[3];
  float *Pr = new float[3];
  float *Pk = new float[3];
    
  Pi[0] = mp+TKEA;
  Pi[1] = 0;
  Pi[2] = sqrt(2*mp*TKEA + TKEA*TKEA);
    
  Pt[0] = mass;
  Pt[1] = 0;
  Pt[2] = 0;
    
  Pr[0] = sqrt(pow(Sp+mass-mp,2)+k*k);
  Pr[1] = k*sin(theta_k);
  Pr[2] = k*cos(theta_k);
    
  Pk[0] = Pt[0] - Pr[0];
  Pk[1] = Pt[1] - Pr[1];
  Pk[2] = Pt[2] - Pr[2];

  float *Pc = new float[3];
  Pc[0] = (Pi[0] + Pk[0])/2;
  Pc[1] = (Pi[1] + Pk[1])/2;
  Pc[2] = (Pi[2] + Pk[2])/2;

  // ########################### 2-body CM frame
  float *Pic = new float[3];
  float *P1c = new float[3];
  float *P2c = new float[3];
  float *Pkc = new float[3];

  float betaPc = Momentum2D(Pc)/Pc[0];
  float thetaPc = Angle(Pc);
    
  Pic = Lorentz2D(Pi, -betaPc, thetaPc);
  Pkc = Lorentz2D(Pk, -betaPc, thetaPc);
    
  float thetaPic = Angle(Pic);
    
  float Etotal = sqrt(1-pow(betaPc,2))*(Pi[0]+Pk[0]); // CM frame total energy
  float p1c = sqrt((Etotal-mp-mp)*(Etotal+mp+mp))/2;
    
  P1c[0] = sqrt(mp*mp + p1c*p1c);
  P1c[1] = p1c*sin(thetaPic+theta_NN);
  P1c[2] = p1c*cos(thetaPic+theta_NN);
    
  P2c[0] = P1c[0];
  P2c[1] = -P1c[1];
  P2c[2] = -P1c[2];
    
  P1 = Lorentz2D(P1c, betaPc, thetaPc);
  P2 = Lorentz2D(P2c, betaPc, thetaPc);
  
  // ########################### Lab Frame Inversed Kinematics
  float *PiL = new float[3];
  float *P1L = new float[3];
  float *P2L = new float[3];
  float *PtL = new float[3];
  float *PrL = new float[3];
    
  float beta = Momentum2D(Pi)/Pi[0];
    
  PiL = Lorentz2D(Pi, -beta, 0);
  PtL = Lorentz2D(Pt, -beta, 0);
  P1L = Lorentz2D(P1, -beta, 0);
  P2L = Lorentz2D(P2, -beta, 0);
  PrL = Lorentz2D(Pr, -beta, 0);
    
  float *output = new float[9];
    
  output[0] = P1[0] - mp; // T_1
  output[1] = Angle(P1)*rad2deg; // theta_1
  output[2] = P2[0] - mp; // T_2
  output[3] = Angle(P2)*rad2deg; // theta_2

  output[4] = P1L[0] - mp;
  output[5] = 180 - Angle(P1L)*rad2deg ; // invert the z-axis 
  output[6] = P2L[0] - mp;
  output[7] = 180 + Angle(P2L)*rad2deg;

  output[8] = -beta;

  /*// display 
  printf("----------------- Lab frame\n");
  PrintV2D(PiL, " PiL");
  PrintV2D(PtL, " PtL");
  PrintV2D(PrL, " PrL");
  PrintV2D(P1L, " P1L");
  PrintV2D(P2L, " P2L");
  printf("T1L:%10.4f, theta_1L:%10.4f, T2L:%10.4f, theta_2L:%10.4f \n", output[4], output[5], output[6], output[7]);
  printf("----------------- Nuclues frame\n");
  PrintV2D(Pi, " Pi");
  PrintV2D(Pk, " Pk");
  PrintV2D(Pt, " Pt");
  PrintV2D(Pr, " Pr");
  PrintV2D(P1, " P1");
  PrintV2D(P2, " P2");
  printf("Tc:%10.4f, theta_c:%10.4f, Td:%10.4f, theta_d:%10.4f \n", output[0], output[1], output[2], output[3]);
  printf("k : %10.3f, theta_k:%10.3f, theta_NN:%10.3f, Sp:%10.5f\n", k, theta_k*rad2deg, theta_NN*rad2deg, Sp);
  printf("Etotal : %10.3f, 2*mp:%10.3f,  p1c:%10.3f, beta:%10.8f\n", Etotal,2*mp, p1c, beta);
  PrintV2D(Pic, " Pic");
  PrintV2D(Pkc, " Pkc");
  PrintV2D(P1c, " P1c");
  PrintV2D(P2c, " P2c");
  printf("-----------------------------------\n");
  */
  //delete pointer
  delete Pi;
  delete Pt;
  delete P1;
  delete P2;
  delete Pr;
  delete Pk;
  delete Pc;
    
  delete Pic;
  delete Pkc;
  delete P1c;
  delete P2c;

  delete PiL;
  delete PtL;
  delete P1L;
  delete P2L;
  delete PrL;
    
  return output;
}

//################################################################################################

float* Knockout2Dinv(int MA, int Z, float TKEA, float T1, float theta_1, float T2, float theta_2){
    
  float mass = Nucleus_Mass(Z, MA);

  if ( mass == -404) return(0);
  
  theta_1 = (180.-theta_1)/rad2deg;
  theta_2 = (theta_2-180.)/rad2deg;
  
  // ########################### Lab Frame Inversed Kinematics
  float *PiL = new float[3];
  float *P1L = new float[3];
  float *P2L = new float[3];
  float *PtL = new float[3];
  float *PrL = new float[3];
      
  PiL[0] = mp;
  PiL[1] = 0;
  PiL[2] = 0;

  float beta = sqrt(2*mp*TKEA + TKEA*TKEA)/(mp + TKEA);  

  P1L[0] = mp + T1;
  P1L[1] = sqrt(2*mp*T1 + T1*T1) * sin(theta_1);
  P1L[2] = sqrt(2*mp*T1 + T1*T1) * cos(theta_1);

  P2L[0] = mp + T2;
  P2L[1] = sqrt(2*mp*T2 + T2*T2) * sin(theta_2);
  P2L[2] = sqrt(2*mp*T2 + T2*T2) * cos(theta_2);


  //  ########################### Nucleus's frame
  float *Pi = new float[3];
  float *Pt = new float[3];
  float *P1 = new float[3];
  float *P2 = new float[3];
  float *Pr = new float[3];
  float *Pk = new float[3];
    
  Pi[0] = mp+TKEA;
  Pi[1] = 0;
  Pi[2] = sqrt(2*mp*TKEA + TKEA*TKEA);
    
  Pt[0] = mass;
  Pt[1] = 0;
  Pt[2] = 0;

  PtL = Lorentz2D(Pt, -beta, 0);
  P1  = Lorentz2D(P1L, beta, 0);
  P2  = Lorentz2D(P2L, beta, 0);

  PrL[0] = PiL[0]+PtL[0]-P1L[0]-P2L[0];
  PrL[1] = PiL[1]+PtL[1]-P1L[1]-P2L[1];
  PrL[2] = PiL[2]+PtL[2]-P1L[2]-P2L[2];

  Pr = Lorentz2D(PrL, beta, 0);
    
  Pk[0] = Pt[0] - Pr[0];
  Pk[1] = Pt[1] - Pr[1];
  Pk[2] = Pt[2] - Pr[2];

  float *Pc = new float[3];

  Pc[0] = (P1[0] + P2[0])/2;
  Pc[1] = (P1[1] + P2[1])/2;
  Pc[2] = (P1[2] + P2[2])/2;

  // ########################### 2-body CM frame
  float *Pic = new float[3];
  float *P1c = new float[3];
  float *P2c = new float[3];
  float *Pkc = new float[3];

  float betaPc = Momentum2D(Pc)/Pc[0];
  float thetaPc = Angle(Pc);
    
  Pic = Lorentz2D(Pi, -betaPc, thetaPc);
  Pkc = Lorentz2D(Pk, -betaPc, thetaPc);    
  P1c = Lorentz2D(P1, -betaPc, thetaPc);
  P2c = Lorentz2D(P2, -betaPc, thetaPc);
    
  //  ########################### Output
  float *output = new float[9];
    
  output[0] = Momentum2D(Pk); // k
  output[1] = Angle(Pk)*rad2deg; // theta_k
  output[2] = acos((Pic[1]*P1c[1]+Pic[2]*P1c[2])/Momentum2D(Pic)/Momentum2D(P1c))*rad2deg; // theta_NN
  output[3] = Mass2D(Pr) - mass + mp; // Sp

  output[4] = P1[0] - mp;   // Tc
  output[5] = Angle(P1)*rad2deg ; // theta_c 
  output[6] = P2[0] - mp;   // Td
  output[7] = Angle(P2)*rad2deg;  // theta_d

  output[8] = -beta;

  float gamma = 1/sqrt(1-beta*beta);
  float SpOnline = mp*(1 - gamma) - gamma*(T1+T2) + gamma*beta*(abs(P1L[2]+P2L[2])); 

  /*  // display 
  printf("----------------- Lab frame\n");
  PrintV2D(PiL, " PiL");
  PrintV2D(PtL, " PtL");
  PrintV2D(PrL, " PrL");
  PrintV2D(P1L, " P1L");
  PrintV2D(P2L, " P2L");
  printf("T1L:%10.4f, theta_1L:%10.4f, T2L:%10.4f, theta_2L:%10.4f \n", T1, 180-theta_1*rad2deg, T2, theta_2*rad2deg+180);
  printf("----------------- Nuclues frame\n");
  PrintV2D(Pi, " Pi");
  PrintV2D(Pk, " Pk");
  PrintV2D(Pt, " Pt");
  PrintV2D(Pr, " Pr");
  PrintV2D(P1, " P1");
  PrintV2D(P2, " P2");
  printf("Tc:%10.4f, theta_c:%10.4f, Td:%10.4f, theta_d:%10.4f \n", output[4], output[5], output[6], output[7]);
  printf("k : %10.3f, theta_k:%10.3f, theta_NN:%10.3f, Sp:%10.5f\n", output[0],output[1], output[2], output[3]);
  printf(" beta:%10.8f, Sp_Online:%10.5f\n", beta, SpOnline);
  PrintV2D(Pic, " Pic");
  PrintV2D(Pkc, " Pkc");
  PrintV2D(P1c, " P1c");
  PrintV2D(P2c, " P2c");
  printf("-----------------------------------\n");
  */

  //delete pointer
  delete Pi;
  delete Pt;
  delete P1;
  delete P2;
  delete Pr;
  delete Pk;
  delete Pc;
    
  delete Pic;
  delete Pkc;
  delete P1c;
  delete P2c;

  delete PiL;
  delete PtL;
  delete P1L;
  delete P2L;
  delete PrL;
    
  return output;
}


//################################################################################################

float* Knockout2Dinv2(int MA, int Z, float TKEA, float T1, float theta_1, float theta_2, float Sp){
    
  printf("-------------------------------------\n");
  printf(" T1: %10.5f, theta_1: %10.5f, theta_2: %10.5f, SP:%10.5f\n", T1, theta_1, theta_2, Sp);
  float mass = Nucleus_Mass(Z, MA);

  if ( mass == -404) return(0);
  
  theta_1 = (180.-theta_1)/rad2deg;
  theta_2 = (theta_2-180.)/rad2deg;
  
  // ########################### Lab Frame Inversed Kinematics
  float *PiL = new float[3];
  float *P1L = new float[3];
  float *P2L = new float[3];
  float *PtL = new float[3];
  float *PrL = new float[3];
      
  PiL[0] = mp;
  PiL[1] = 0;
  PiL[2] = 0;

  float beta = sqrt(2*mp*TKEA + TKEA*TKEA)/(mp + TKEA);  
  float gamma = 1/sqrt(1-beta*beta);

  P1L[0] = mp + T1;
  float p1 = sqrt(2*mp*T1 + T1*T1);
  P1L[1] = p1 * sin(theta_1);
  P1L[2] = p1 * cos(theta_1);

  float a = PiL[0] + PtL[0] - P1L[0];
  float mr = mass + Sp - mp;
  float b = mr*mr+ PtL[2]*PtL[2] + p1*p1 - 2*PtL[2]*P1L[2];
  float d = 2*p1*cos(theta_1-theta_2) - 2*PtL[2]*cos(theta_2);

  float A = 4*a*a - d*d;
  float B = 2*(a*a + mp*mp - b)*d;
  float C = 4*a*a*mp*mp - (a*a + mp*mp- b)*(a*a + mp*mp- b);

  float p2_a = (-B+sqrt(B*B-4*A*C))/(2*A);
  float p2_b = (-B-sqrt(B*B-4*A*C))/(2*A);

  float T2_a = sqrt(mp*mp + p2_a*p2_a) - mp; 
  float T2_b = sqrt(mp*mp + p2_b*p2_b) - mp;

  printf("p2 : %10.4f , %10.4f \n", p2_a, p2_b);
  printf("T2 : %10.4f , %10.4f \n", T2_a, T2_b);

  float T2 = T2_b;
  float p2 = p2_b;

  printf("T2: %10.5f, p2: %10.5f, theta_2: %10.5f\n",T2, p2,  theta_2*rad2deg);

  P2L[0] = mp + T2;
  P2L[1] = p2 * sin(theta_2);
  P2L[2] = p2 * cos(theta_2);

  //  ########################### Nucleus's frame
  float *Pi = new float[3];
  float *Pt = new float[3];
  float *P1 = new float[3];
  float *P2 = new float[3];
  float *Pr = new float[3];
  float *Pk = new float[3];
    
  Pi[0] = mp+TKEA;
  Pi[1] = 0;
  Pi[2] = sqrt(2*mp*TKEA + TKEA*TKEA);
    
  Pt[0] = mass;
  Pt[1] = 0;
  Pt[2] = 0;

  PtL = Lorentz2D(Pt, -beta, 0);
  P1  = Lorentz2D(P1L, beta, 0);
  P2  = Lorentz2D(P2L, beta, 0);

  PrL[0] = PiL[0]+PtL[0]-P1L[0]-P2L[0];
  PrL[1] = PiL[1]+PtL[1]-P1L[1]-P2L[1];
  PrL[2] = PiL[2]+PtL[2]-P1L[2]-P2L[2];

  Pr = Lorentz2D(PrL, beta, 0);
    
  Pk[0] = Pt[0] - Pr[0];
  Pk[1] = Pt[1] - Pr[1];
  Pk[2] = Pt[2] - Pr[2];

  float *Pc = new float[3];

  Pc[0] = (P1[0] + P2[0])/2;
  Pc[1] = (P1[1] + P2[1])/2;
  Pc[2] = (P1[2] + P2[2])/2;

  // ########################### 2-body CM frame
  float *Pic = new float[3];
  float *P1c = new float[3];
  float *P2c = new float[3];
  float *Pkc = new float[3];

  float betaPc = Momentum2D(Pc)/Pc[0];
  float thetaPc = Angle(Pc);
    
  Pic = Lorentz2D(Pi, -betaPc, thetaPc);
  Pkc = Lorentz2D(Pk, -betaPc, thetaPc);    
  P1c = Lorentz2D(P1, -betaPc, thetaPc);
  P2c = Lorentz2D(P2, -betaPc, thetaPc);
    
  //  ########################### Output
  float *output = new float[9];
    
  output[0] = Momentum2D(Pk); // k
  output[1] = Angle(Pk)*rad2deg; // theta_k
  output[2] = acos((Pic[1]*P1c[1]+Pic[2]*P1c[2])/Momentum2D(Pic)/Momentum2D(P1c))*rad2deg; // theta_NN
  output[3] = Mass2D(Pr) - mass + mp; // Sp

  output[4] = P1[0] - mp;   // Tc
  output[5] = Angle(P1)*rad2deg ; // theta_c 
  output[6] = P2[0] - mp;   // Td
  output[7] = Angle(P2)*rad2deg;  // theta_d

  output[8] = -beta;

  float SpOnline = mp*(1 - gamma) - gamma*(T1+T2) + gamma*beta*(abs(P1L[2]+P2L[2])); 

  /* // display 
  printf("----------------- Lab frame\n");
  PrintV2D(PiL, " PiL");
  PrintV2D(PtL, " PtL");
  PrintV2D(PrL, " PrL");
  PrintV2D(P1L, " P1L");
  PrintV2D(P2L, " P2L");
  printf("T1L:%10.4f, theta_1L:%10.4f, T2L:%10.4f, theta_2L:%10.4f \n", T1, 180-theta_1*rad2deg, T2, theta_2*rad2deg+180);
  printf("----------------- Nuclues frame\n");
  PrintV2D(Pi, " Pi");
  PrintV2D(Pk, " Pk");
  PrintV2D(Pt, " Pt");
  PrintV2D(Pr, " Pr");
  PrintV2D(P1, " P1");
  PrintV2D(P2, " P2");
  printf("Tc:%10.4f, theta_c:%10.4f, Td:%10.4f, theta_d:%10.4f \n", output[4], output[5], output[6], output[7]);
  printf("k : %10.3f, theta_k:%10.3f, theta_NN:%10.3f, Sp:%10.5f\n", output[0],output[1], output[2], output[3]);
  printf(" beta:%10.8f, Sp_Online:%10.5f\n", beta, SpOnline);
  PrintV2D(Pic, " Pic");
  PrintV2D(Pkc, " Pkc");
  PrintV2D(P1c, " P1c");
  PrintV2D(P2c, " P2c");
  printf("-----------------------------------\n");
  */

  //delete pointer
  delete Pi;
  delete Pt;
  delete P1;
  delete P2;
  delete Pr;
  delete Pk;
  delete Pc;
    
  delete Pic;
  delete Pkc;
  delete P1c;
  delete P2c;

  delete PiL;
  delete PtL;
  delete P1L;
  delete P2L;
  delete PrL;
    
  return output;
}
//################################################################################################

float* Knockout2Dinv3(int MA, int Z, float TKEA, float Tc, float theta_c, float theta_d, float Sp){
    
  //printf("---------------- knockout inv 3---------------------\n");
  //printf(" Tc: %10.5f, theta_c: %10.5f, theta_d: %10.5f, SP:%10.5f\n", Tc, theta_c, theta_d, Sp);
  float mass = Nucleus_Mass(Z, MA);

  if ( mass == -404) return(0);
  
  theta_c = theta_c/rad2deg;
  theta_d = theta_d/rad2deg;

  //  ########################### Nucleus's frame
  float *Pi = new float[3];
  float *Pt = new float[3];
  float *P1 = new float[3];
  float *P2 = new float[3];
  float *Pr = new float[3];
  float *Pk = new float[3];
    
  Pi[0] = mp+TKEA;
  Pi[1] = 0;
  Pi[2] = sqrt(2*mp*TKEA + TKEA*TKEA);
    
  Pt[0] = mass;
  Pt[1] = 0;
  Pt[2] = 0;

  P1[0] = mp + Tc;
  float p1 = sqrt(2*mp*Tc + Tc*Tc);
  P1[1] = p1*sin(theta_c);
  P1[2] = p1*cos(theta_c);

  float mr = mass + Sp - mp;
  //calculate Td
  float a = Pi[0] + Pt[0] - P1[0];
  float b = mr*mr + Pi[2]*Pi[2] +  p1*p1 - 2*Pi[2]*P1[2];
  float d = 2*p1*cos(theta_c-theta_d) - 2*Pi[2]*cos(theta_d);

  float A = 4*a*a-d*d;
  float B = 2*(a*a+mp*mp-b)*d;
  float C = 4*a*a*mp*mp-(a*a+mp*mp-b)*(a*a+mp*mp-b);

  float p2a = (-B-sqrt(B*B-4*A*C))/2/A;
  float p2b = (-B+sqrt(B*B-4*A*C))/2/A;

  float p2 = p2b;

  float T2a = 0;
  if ( p2a >= 0 ) T2a = sqrt(mp*mp + p2a*p2a)-mp;
  float T2 = sqrt(mp*mp + p2*p2)-mp;

  P2[0] = mp + T2;
  P2[1] = p2 * sin(theta_d);
  P2[2] = p2 * cos(theta_d);

  Pr[0] = Pi[0]+Pt[0]-P1[0]-P2[0];
  Pr[1] = Pi[1]+Pt[1]-P1[1]-P2[1];
  Pr[2] = Pi[2]+Pt[2]-P1[2]-P2[2];
    
  Pk[0] = Pt[0] - Pr[0];
  Pk[1] = Pt[1] - Pr[1];
  Pk[2] = Pt[2] - Pr[2];

  float *Pc = new float[3];

  Pc[0] = (P1[0] + P2[0])/2;
  Pc[1] = (P1[1] + P2[1])/2;
  Pc[2] = (P1[2] + P2[2])/2;
  
  // ########################### Lab Frame Inversed Kinematics
  float *PiL = new float[3];
  float *P1L = new float[3];
  float *P2L = new float[3];
  float *PtL = new float[3];
  float *PrL = new float[3];

  float beta = Pi[2]/Pi[0];  

  PiL = Lorentz2D(Pi, -beta, 0);
  PtL = Lorentz2D(Pt, -beta, 0);
  PrL = Lorentz2D(Pr, -beta, 0);
  P1L = Lorentz2D(P1, -beta, 0);
  P2L = Lorentz2D(P2, -beta, 0);

  // ########################### 2-body CM frame
  float *Pic = new float[3];
  float *P1c = new float[3];
  float *P2c = new float[3];
  float *Pkc = new float[3];

  float betaPc = Momentum2D(Pc)/Pc[0];
  float thetaPc = Angle(Pc);
    
  Pic = Lorentz2D(Pi, -betaPc, thetaPc);
  Pkc = Lorentz2D(Pk, -betaPc, thetaPc);    
  P1c = Lorentz2D(P1, -betaPc, thetaPc);
  P2c = Lorentz2D(P2, -betaPc, thetaPc);
    
  //  ########################### Output
  float *output = new float[9];
    
  output[0] = Momentum2D(Pk); // k
  output[1] = Angle(Pk)*rad2deg; // theta_k
  output[2] = acos((Pic[1]*P1c[1]+Pic[2]*P1c[2])/Momentum2D(Pic)/Momentum2D(P1c))*rad2deg; // theta_NN
  output[3] = T2; // Td

  output[4] = P1L[0] - mp;   // T1
  output[5] = 180 - Angle(P1L)*rad2deg; // theta_1 
  output[6] = P2L[0] - mp;   // T2
  output[7] = 180 + Angle(P2L)*rad2deg;  // theta_2

  output[8] = T2a;

  //float SpOnline = 0; 

  /* // display 
  printf("----------------- Lab frame\n");
  PrintV2D(PiL, " PiL");
  PrintV2D(PtL, " PtL");
  PrintV2D(PrL, " PrL");
  PrintV2D(P1L, " P1L");
  PrintV2D(P2L, " P2L");
  printf("T1:%10.4f, theta_1:%10.4f, T2:%10.4f, theta_2:%10.4f \n", output[4], output[5], output[6], output[7]);
  printf("----------------- Nuclues frame\n");
  PrintV2D(Pi, " Pi");
  PrintV2D(Pk, " Pk");
  PrintV2D(Pt, " Pt");
  PrintV2D(Pr, " Pr");
  PrintV2D(P1, " P1");
  PrintV2D(P2, " P2");
  printf("Tc:%10.4f, theta_c:%10.4f, theta_d:%10.4f, Sp:%7.3f \n", Tc, theta_c*rad2deg, theta_d*rad2deg, Sp);
  printf("k : %10.3f, theta_k:%10.3f, theta_NN:%10.3f, Td:%10.5f\n", output[0],output[1], output[2], output[3]);
  //printf(" beta:%10.8f, Sp_Online:%10.5f\n", beta, SpOnline);
  PrintV2D(Pic, " Pic");
  PrintV2D(Pkc, " Pkc");
  PrintV2D(P1c, " P1c");
  PrintV2D(P2c, " P2c");
  printf("-----------------------------------\n");
  */

  //delete pointer
  delete Pi;
  delete Pt;
  delete P1;
  delete P2;
  delete Pr;
  delete Pk;
  delete Pc;
    
  delete Pic;
  delete Pkc;
  delete P1c;
  delete P2c;

  delete PiL;
  delete PtL;
  delete P1L;
  delete P2L;
  delete PrL;
    
  return output;
}
