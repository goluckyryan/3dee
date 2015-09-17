/*  Knockout of T(p,2p)R

The Lorentz Tranform on aribtary direction is well discrpted in Goldstien P .281

k(vector): momentum of a proton in Target nucleus
Theta_k:   direction of the momentum of a proton in O nucleus
Theta_NN:  scattering angle in 2-nucleon center of mass system
S:         separation energys

***** suffix notation ****
suffix "L" is Lab frame
suffix "" is Target frame
suffix "c" is CM frame
suffix "X" is experiment

***** notation ****
k is for momentum
P is for 4-momentum
E is Total relativistic energy
T is for Kinetic energy

***********************************/

#include <sstream> // enable string control
#include <cmath> // math , enable simple math function
#include <iomanip>
#include <stdio.h>  // standard I/O
#include <stdlib.h> // standard General Utlilities Library
#include "constant.h"
#include "nuclei_mass.h"

using namespace std; // declare a namespace "std", every variable in this code is inside "std"

// 4 vector = (E, px, py, pz)

  float *Pi = new float[4];
  float *Pt = new float[4];
  float *P1 = new float[4];
  float *P2 = new float[4];
  float *Pr = new float[4];
  float *Pk = new float[4];

  float *Pc = new float[4];

  float *Pic = new float[4];
  float *P1c = new float[4];
  float *P2c = new float[4];
  float *Pkc = new float[4];
  float *Pcc = new float[4];

  float *PiL = new float[4];
  float *P1L = new float[4];
  float *P2L = new float[4];
  float *PtL = new float[4];
  float *PrL = new float[4];

  float * Prot = new float[4];

//function___________________________________________________________________
int Lorentz(float *V, float beta, float theta, float phi, float *U){

  float gamma = 1/sqrt(1-beta*beta);
  float alpha = beta*gamma;
  float cTh   = cos(theta); float ccTh = cTh*cTh;
  float sTh   = sin(theta); float ssTh = sTh*sTh;
  float cPh   = cos(phi); float ccPh = cPh*cPh;
  float sPh   = sin(phi); float ssPh = sPh*sPh;

  float L11 = gamma;
  float L12 = alpha*sTh*cPh;
  float L13 = alpha*sTh*sPh;
  float L14 = alpha*cTh;

  float L22 = ccPh*(ccTh + gamma*ssTh)+ssPh;
  float L23 = cPh*sPh*(ccTh+gamma*ssTh-1);
  float L24 = cTh*sTh*cPh*(gamma-1);

  float L33 = ccPh + ssPh*(ccTh + gamma*ssTh);
  float L34 = cTh*sTh*sPh*(gamma-1);

  float L44 = gamma*ccTh + ssTh;

  U[0] = V[0]*L11 + V[1]*L12 + V[2]*L13 + V[3]*L14;
  U[1] = V[0]*L12 + V[1]*L22 + V[2]*L23 + V[3]*L24;
  U[2] = V[0]*L13 + V[1]*L23 + V[2]*L33 + V[3]*L34;
  U[3] = V[0]*L14 + V[1]*L24 + V[2]*L34 + V[3]*L44;
    
  return 1;
}

void PrintV(float *V, string Vname){
  printf("%10s : %10.3f, %10.3f, %10.3f, %10.3f\n", Vname.c_str(), V[0], V[1], V[2], V[3]); 
}

float Momentum(float *V){
  return sqrt(V[1]*V[1] + V[2]*V[2] + V[3]*V[3]);
}

float Theta(float *V){
  return atan2(sqrt(V[1]*V[1] + V[2]*V[2]), V[3]);
}

float Phi(float *V){
  return atan2(V[2], V[1]);
}

float DAng(float *V, float *U){
  float dotproduct = V[1]*U[1]+V[2]*U[2]+V[3]*U[3];
  dotproduct = dotproduct/Momentum(V)/Momentum(U);
  return acos(dotproduct);
}

/*
float* Knockout3D(int MA, int Z, float TKEA, float k, float theta_k, float phi_k, float theta_NN, float phi_NN, float Sp){
  /*  
  printf("-----KNOCKOUT 3D-----------\n");
  printf("     k: %10.6f,     Sp: %10.6f \n", k, Sp);
  printf(" ang_k: %10.6f, ang_NN: %10.6f \n", theta_k, theta_NN);
  printf(" phi_k: %10.6f, phi_NN: %10.6f \n", phi_k, phi_NN);
  
  float mass = Nucleus_Mass(Z, MA);

  if ( mass == -404) return(0);
  
  theta_k = theta_k/rad2deg;
  theta_NN = theta_NN/rad2deg;
  phi_k = phi_k/rad2deg;
  phi_NN = phi_NN/rad2deg;
 
  //  ########################### Nuclei's frame
  float *Pi = new float[4];
  float *Pt = new float[4];
  float *P1 = new float[4];
  float *P2 = new float[4];
  float *Pr = new float[4];
  float *Pk = new float[4];
    
  Pi[0] = mp+TKEA;
  Pi[1] = 0;
  Pi[2] = 0;
  Pi[3] = sqrt(2*mp*TKEA + TKEA*TKEA);
    
  Pt[0] = mass;
  Pt[1] = 0;
  Pt[2] = 0;
  Pt[3] = 0;
    
  Pr[0] = sqrt(pow(Sp+mass-mp,2)+k*k);
  Pr[1] = k*sin(theta_k)*cos(phi_k);
  Pr[2] = k*sin(theta_k)*sin(phi_k);
  Pr[3] = k*cos(theta_k);
    
  Pk[0] = Pt[0] - Pr[0];
  Pk[1] = Pt[1] - Pr[1];
  Pk[2] = Pt[2] - Pr[2];
  Pk[3] = Pt[3] - Pr[3];

  float *Pc = new float[4];
  Pc[0] = (Pi[0] + Pk[0])/2;
  Pc[1] = (Pi[1] + Pk[1])/2;
  Pc[2] = (Pi[2] + Pk[2])/2;
  Pc[3] = (Pi[3] + Pk[3])/2;

  // ########################### 2-body CM frame
  float *Pic = new float[4];
  float *P1c = new float[4];
  float *P2c = new float[4];
  float *Pkc = new float[4];
  float *Pcc = new float[4];

  float betaPc = Momentum(Pc)/Pc[0];
  float thetaPc = Theta(Pc);
  float phiPc = Phi(Pc);
  
  Pcc = Lorentz(Pc, -betaPc, thetaPc, phiPc);
    
  Pic = Lorentz(Pi, -betaPc, thetaPc, phiPc);
  Pkc = Lorentz(Pk, -betaPc, thetaPc, phiPc);
    
  float thetaPic = Theta(Pic);
  float phiPic = Phi(Pic);
    
  float Etotal = sqrt(1-pow(betaPc,2))*(Pi[0]+Pk[0]); // CM frame total energy
  float p1c = sqrt((Etotal-mp-mp)*(Etotal+mp+mp))/2;
    
  P1c[0] = sqrt(mp*mp + p1c*p1c);
  P1c[1] = p1c*(cos(theta_NN)*cos(phiPic)*sin(thetaPic)+sin(theta_NN)*(cos(thetaPic)*cos(phiPic)*cos(phi_NN) - sin(phiPic)*sin(phi_NN)));
  P1c[2] = p1c*(cos(theta_NN)*sin(phiPic)*sin(thetaPic)+sin(theta_NN)*(cos(thetaPic)*sin(phiPic)*cos(phi_NN) + cos(phiPic)*sin(phi_NN)));
  P1c[3] = p1c*(cos(theta_NN)*cos(thetaPic) - cos(phi_NN) * sin(thetaPic) * sin(theta_NN));
    
  P2c[0] = P1c[0];
  P2c[1] = -P1c[1];
  P2c[2] = -P1c[2];
  P2c[3] = -P1c[3];
    
  P1 = Lorentz(P1c, betaPc, thetaPc, phiPc);
  P2 = Lorentz(P2c, betaPc, thetaPc, phiPc);
  
  // ########################### Lab Frame Inversed Kinematics
  float *PiL = new float[4];
  float *P1L = new float[4];
  float *P2L = new float[4];
  float *PtL = new float[4];
  float *PrL = new float[4];
  float *PkL = new float[4];
    
  float beta = Momentum(Pi)/Pi[0];
    
  PiL = Lorentz(Pi, -beta, 0, 0);
  PtL = Lorentz(Pt, -beta, 0, 0);
  P1L = Lorentz(P1, -beta, 0, 0);
  P2L = Lorentz(P2, -beta, 0, 0);
  PrL = Lorentz(Pr, -beta, 0, 0);
  
  PkL[0] = PtL[0] - PrL[0];
  PkL[1] = PtL[1] - PrL[1];
  PkL[2] = PtL[2] - PrL[2];
  PkL[3] = PtL[3] - PrL[3];
    
  float *output = new float[13];
    
  output[0] = P1[0] - mp;        // T_c
  output[1] = Theta(P1)*rad2deg; // theta_c
  output[2] = Phi(P1)*rad2deg;
  output[3] = P2[0] - mp;          // T_d
  output[4] = - Theta(P2)*rad2deg; // theta_d
  output[5] = Phi(P2)*rad2deg;
  output[6] = output[2]-output[5];// offplane angle, rotation around Kd x(Ka x Kc)

  output[7] = P1L[0] - mp;
  output[8] = Theta(P1L)*rad2deg;
  output[9] = Phi(P1L)*rad2deg;
  output[10] = P2L[0] - mp;
  output[11] = Theta(P2L)*rad2deg;
  output[12] = Phi(P2L)*rad2deg;

  // display 
  printf("----------------- Lab frame\n");
  PrintV(PiL, " PiL");
  PrintV(PtL, " PtL");
  PrintV(PrL, " PrL");
  PrintV(PkL, " PkL");
  PrintV(P1L, " P1L");
  PrintV(P2L, " P2L");
  printf("T1L:%10.4f, theta_1L:%10.4f, phi_1L:%10.4f \nT2L:%10.4f, theta_2L:%10.4f, phi_2L:%10.4f\n", output[7], output[8], output[9], output[10], output[11], output[12]);
  printf("----------------- Nuclues frame\n");
  PrintV(Pi, " Pi");
  PrintV(Pk, " Pk");
  PrintV(Pt, " Pt");
  PrintV(Pr, " Pr");
  PrintV(P1, " P1");
  PrintV(P2, " P2");
  PrintV(Pc, " Pc");
  printf("betaPc: %10.4f, thetaPc: %10.4f, phiPc: %10.4f \n", betaPc, thetaPc*rad2deg, phiPc*rad2deg);
  printf("Tc:%10.4f, theta_c:%10.4f, phi_c:%10.4f \nTd:%10.4f, theta_d:%10.4f, phi_d:%10.4f\noffplane ang:%10.4f\n", output[0], output[1], output[2], output[3], output[4], output[5], output[6]);
  printf("k : %10.3f, theta_k:%10.3f, theta_NN:%10.3f, Sp:%10.5f\n", k, theta_k*rad2deg, theta_NN*rad2deg, Sp);
  printf("Etotal : %10.3f, 2*mp:%10.3f,  p1c:%10.3f, beta:%10.8f\n", Etotal,2*mp, p1c, beta);
  PrintV(Pcc, " Pcc");
  PrintV(Pic, " Pic");
  PrintV(Pkc, " Pkc");
  PrintV(P1c, " P1c");
  PrintV(P2c, " P2c");
  printf("Dang(Pic, P1c) = %6.2f deg \n", DAng(Pic,P1c)*rad2deg );
  printf("=======================================\n");
  
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
  delete PkL;
    
  return output;
}*/

//#################################################################################

int Knockout3Dinv3(int MA, int Z, float TKEA, float Tc, float theta_c, float phi_c, float theta_d, float phi_d, float Sp, float *output){
  /*
  printf("------KNOCKOUT3D inv3---------\n");
  printf("   T_c: %10.6f,     Sp: %10.6f \n", Tc, Sp);
  printf(" ang_c: %10.6f,  ang_d: %10.6f \n", theta_c, theta_d);
  printf(" phi_c: %10.6f,  phi_d: %10.6f \n", phi_c, phi_d);
  */
  //float mass = Nucleus_Mass(Z, MA);

  float mass;
  switch (MA) {
  case 14: mass = 13044.836422; break;
  case 16: mass = 14895.080096; break;
  case 22: mass = 20498.064784; break;
  case 23: mass = 21423.077578; break;
  case 25: mass = 23294.116162; break;
  case 24: mass = 22370.269976; break; 
  default: mass = -404; break;
  }

  if ( mass == -404) return(0);
  
  theta_c = theta_c/rad2deg;
  theta_d = theta_d/rad2deg;
  phi_c = phi_c/rad2deg;
  phi_d = phi_d/rad2deg;
 
  //  ########################### Nuclei's frame

    
  Pi[0] = mp+TKEA;
  Pi[1] = 0;
  Pi[2] = 0;
  Pi[3] = sqrt(2*mp*TKEA + TKEA*TKEA);
    
  Pt[0] = mass;
  Pt[1] = 0;
  Pt[2] = 0;
  Pt[3] = 0;

  P1[0] = mp+Tc;
  float pc = sqrt(2*mp*Tc+Tc*Tc);
  P1[1] = pc*sin(theta_c)*cos(phi_c);
  P1[2] = pc*sin(theta_c)*sin(phi_c);
  P1[3] = pc*cos(theta_c);

  //Find Td
  float mr = mass + Sp - mp;
  float a = Pi[0] + Pt[0] - P1[0];
  float b = mr*mr + Pi[3]*Pi[3] +  pc*pc - 2*Pi[3]*P1[3];
  float DD = cos(theta_c)*cos(theta_d)+sin(theta_c)*sin(theta_d)*cos(phi_c-phi_d);
  float d = 2*pc*DD - 2*Pi[3]*cos(theta_d);

  float A = 4*a*a-d*d;
  float B = 2*(a*a+mp*mp-b)*d;
  float C = 4*a*a*mp*mp-(a*a+mp*mp-b)*(a*a+mp*mp-b);

  float Delta = B*B - 4*A*C;

  if( Delta < 0 ){
    printf(" Delta < 0, Tc:%10.4f, theta_c:%10.4f, theta_d:%10.4f, phi_c:%10.4f, phi_d:%10.4f \n", Tc, theta_c*rad2deg, theta_d*rad2deg, phi_c*rad2deg, phi_d*rad2deg);
    return 0;
  }

  float p2a = (-B-sqrt(Delta))/2/A;
  float p2b = (-B+sqrt(Delta))/2/A;

  float p2 = p2b;

  float Tda = sqrt(mp*mp + p2a*p2a)-mp;
  if( p2a > 0) printf(" | p2min = %f, Tdmin = %f \n", p2a, Tda); 
  float Td = sqrt(mp*mp + p2*p2)-mp;
  
  P2[0] = mp+Td;
  P2[1] = p2*sin(theta_d)*cos(phi_d);
  P2[2] = p2*sin(theta_d)*sin(phi_d);
  P2[3] = p2*cos(theta_d);

  Pr[0] = Pi[0] + Pt[0] - P1[0] - P2[0];
  Pr[1] = Pi[1] + Pt[1] - P1[1] - P2[1];
  Pr[2] = Pi[2] + Pt[2] - P1[2] - P2[2];  Pr[3] = Pi[3] + Pt[3] - P1[3] - P2[3];
    
  Pk[0] = Pt[0] - Pr[0];
  Pk[1] = Pt[1] - Pr[1];
  Pk[2] = Pt[2] - Pr[2];
  Pk[3] = Pt[3] - Pr[3];

  Pc[0] = (Pi[0] + Pk[0])/2;
  Pc[1] = (Pi[1] + Pk[1])/2;
  Pc[2] = (Pi[2] + Pk[2])/2;
  Pc[3] = (Pi[3] + Pk[3])/2;

  // ########################### 2-body CM frame


  float betaPc = Momentum(Pc)/Pc[0];
  float thetaPc = Theta(Pc);
  float phiPc = Phi(Pc);
  
  Lorentz(Pc, -betaPc, thetaPc, phiPc, Pcc);
  Lorentz(Pi, -betaPc, thetaPc, phiPc, Pic);
  Lorentz(Pk, -betaPc, thetaPc, phiPc, Pkc);
  Lorentz(P1, -betaPc, thetaPc, phiPc, P1c);
  Lorentz(P2, -betaPc, thetaPc, phiPc, P2c);
    
  // ########################### Lab Frame Inversed Kinematics
  float beta = Momentum(Pi)/Pi[0];
    
  Lorentz(Pi, -beta, 0, 0, PiL);
  Lorentz(Pt, -beta, 0, 0, PtL);
  Lorentz(P1, -beta, 0, 0, P1L);
  Lorentz(P2, -beta, 0, 0, P2L);
  Lorentz(Pr, -beta, 0, 0, PrL);

  //Find theta_NNm phi_NN
  float theta_ic = Theta(Pic);
  float phi_ic = Phi(Pic);

  Prot[0] = 0;
  Prot[1] =  (P1c[1]*cos(-phi_ic) - P1c[2]*sin(-phi_ic))*cos(-theta_ic)+P1c[3]*sin(-theta_ic);
  Prot[2] =   P1c[1]*sin(-phi_ic) + P1c[2]*cos(-phi_ic);
  Prot[3] = -(P1c[1]*cos(-phi_ic) - P1c[2]*sin(-phi_ic))*sin(-theta_ic)+P1c[3]*cos(-theta_ic);
  
  //===================================  
    
  output[0] = P1L[0] - mp; // T_1
  output[1] = 180 - Theta(P1L)*rad2deg; // theta_1
  output[2] = Phi(P1L)*rad2deg;
  output[3] = P2L[0] - mp; // T_2
  output[4] = 180 - Theta(P2L)*rad2deg; // theta_2
  output[5] = Phi(P2L)*rad2deg;

  output[6] = phi_c - phi_d;// offplane angle, rotation around Kd x(Ka x Kc)
  output[7] = Td;
  
  output[8] = Momentum(Pk);//k
  output[9] = Theta(Pr)*rad2deg;//theta_k
  output[10] = Phi(Pr)*rad2deg;//phi_k
  output[11] = Theta(Prot) *rad2deg;//theta_NN
  output[12] = Phi(Prot)*rad2deg;//phi_NN

  /*// display 
  printf("----------------- Lab frame\n");
  PrintV(PiL, " PiL");
  PrintV(PtL, " PtL");
  PrintV(PrL, " PrL");
  PrintV(P1L, " P1L");
  PrintV(P2L, " P2L");
  printf("T1L:%10.4f, theta_1L:%10.4f, phi_1L:%10.4f\nT2L:%10.4f, theta_2L:%10.4f, phi_2L:%10.4f\n", output[0], output[1], output[2], output[3], output[4], output[5]);
  printf("----------------- Nuclues frame\n");
  PrintV(Pi, " Pi");
  PrintV(Pk, " Pk");
  PrintV(Pt, " Pt");
  PrintV(Pr, " Pr");
  PrintV(P1, " P1");
  PrintV(P2, " P2");
  PrintV(Pc, " Pc");
  printf("betaPc: %10.4f, thetaPc: %10.4f, phiPc: %10.4f \n", betaPc, thetaPc*rad2deg, phiPc*rad2deg);
  printf("Tc:%10.4f, theta_c:%10.4f, phi_c:%10.4f\nTd:%10.4f, theta_d:%10.4f, phi_d:%10.4f\n", Tc, theta_c*rad2deg, phi_c*rad2deg, output[7], theta_d*rad2deg, phi_d*rad2deg);
  printf("k : %10.3f, theta_k:%10.3f, phi_k%10.3f, theta_NN:%10.3f, phi_NN:%10.3f,Sp:%10.5f\n", output[8], output[9], output[10], output[11], output[12], Sp);
  printf(" off plane angle: %10.4f deg\n", output[6]); 
  PrintV(Pcc, " Pcc");
  PrintV(Pic, " Pic");
  PrintV(Pkc, " Pkc");
  PrintV(P1c, " P1c");
  PrintV(P2c, " P2c");
  printf("=====================================\n");
  */
  //delete pointer

  /*
  delete [] Pi;
  delete [] Pt;
  delete [] P1;
  delete [] P2;
  delete [] Pr;
  delete [] Pk;
  delete [] Pc;
         
  delete [] Pic;
  delete [] Pkc;
  delete [] P1c;
  delete [] P2c;
         
  delete [] PiL;
  delete [] PtL;
  delete [] P1L;
  delete [] P2L;
  delete [] PrL;
         
  delete [] Prot;

  */
    
  return 1;
}


