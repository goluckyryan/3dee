#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <stdlib.h>
#include <cmath>
#include <cstring>
#include "constant.h"
#include "3DeeGenLibrary.h"
#include "XsecTransform.h"
#include "knockout2D.h"
#include "knockout3D.h"

using namespace std;

int main(int argc, char *argv[]){
  
  float *output = new float[13]; // knockout output 
  float *outputINV = new float[13]; // knockout inver output
  
  float TKA = 289.44;
  int Ma = 23;
  int Z  = 8;
  float k = 260.391;
  float angk = 72.753;
  float phik = 14.528;
  float angNN = 88.976;
  float phiNN = -18.994;
  float BE    = 13.3;
  output = Knockout3D(Ma, Z,  TKA, k, angk, phik, angNN, phiNN, BE);

  printf("T1L:%10.4f, theta_1L:%10.4f, phi_1L:%10.4f\n", output[7], output[8], output[9]);
  printf("T2L:%10.4f, theta_2L:%10.4f, phi_2L:%10.4f\n", output[10], output[11], output[12]);

  printf("T_c:%9.3f, theta_c:%9.3f, phi_c:%9.3f\n", output[0], output[1], output[2]);
  printf("T_d:%9.3f, theta_d:%9.3f, phi_d:%9.3f\n", output[3], output[4], output[5]);    

  printf("----- input of Knockout3Dinv3 \n");
  printf("output[4]:%5.1f, output[5]:%5.1f\n", -output[4], output[5]);
  
  outputINV = Knockout3Dinv3(Ma,Z, TKA, output[0], output[1], output[2], -output[4], output[5],BE);

  printf("---------- knockout 3D inv \n");
  printf("T1:%10.4f, theta_1:%10.4f, phi_1:%10.4f\n", outputINV[0], outputINV[1], outputINV[2]);
  printf("T2:%10.4f, theta_2:%10.4f, phi_2:%10.4f\n", outputINV[3], outputINV[4], outputINV[5]);
  printf(" k:%10.4f, theta_k:%10.4f,  phi_k:%10.4f\n", outputINV[8], outputINV[9], outputINV[10]);
  printf("Sp:%10.4f,  ang_NN:%10.4f, phi_NN:%10.4f\n", BE, outputINV[11], outputINV[12]);

  /*float theta1 = atof(argv[1]);
  float phi1   = atof(argv[2]);
  float theta2 = atof(argv[3]);
  float phi2   = atof(argv[4]);
  */

  printf("T1:%7.1f, theta1: %7.1f, phi1: %7.1f\n", outputINV[0], outputINV[1], outputINV[2]);
  printf("T2:%7.1f, theta2: %7.1f, phi2: %7.1f\n", outputINV[3], outputINV[4], outputINV[5]);
  
  printf("--- Accpetance Filter --- %d \n",AccpetanceFilter3D(outputINV[0], outputINV[1], outputINV[2], outputINV[3], outputINV[4], outputINV[5]));
  
  return 0;

}


