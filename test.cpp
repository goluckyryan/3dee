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
  /*
  float *output = new float[9]; // knockout output 
  float *outputINV = new float[9]; // knockout inver output
  
  float k = 270;
  float angk = 50;
  float phik = 0;
  float angNN = 90;
  float phiNN = 00;
  float BE    = 4.63;
  output = Knockout2D(16, 8,  200, k, angk, angNN, BE);
 
  printf("k:%6.2f angk:%6.2f angNN:%6.2f  Sp: %6.2f\n ---> T_c:%9.3f, theta_c:%9.3f, theta_d:%9.3f\n",
         k, angk, angNN, BE,  output[0], output[1], output[3]);    
  printf("T1L:%10.4f, theta_1L:%10.4f, T2L:%10.4f, theta_2L:%10.4f \n", output[4], output[5], output[6], output[7]);

  outputINV = Knockout2Dinv3(16,8, 200, output[0], output[1], output[3],BE);
  printf("Tc:%10.4f, theta_c:%10.4f, Td:%10.4f, theta_d:%10.4f \n", outputINV[4], outputINV[5], outputINV[6], outputINV[7]);

  */

  float theta1 = atof(argv[1]);
  float phi1   = atof(argv[2]);
  float theta2 = atof(argv[3]);
  float phi2   = atof(argv[4]);

  printf("theta1: %5.1f, phi1: %5.1f, theta2: %5.1f, phi2: %5.1f \n", theta1, phi1, theta2, phi2);

  AccpetanceFilter3D(100, theta1, phi1, 100, theta2, phi2);

  return 0;

}


