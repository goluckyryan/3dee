#include <string>
#include <sstream>
#include <cstdlib>
#include <stdlib.h>
#include <cmath>
#include <cstring>
#include "XsecTransform.h"
#include "knockout2D.h"

using namespace std;

int main(int argc, char *argv[]){

  /*
  float Tinc = atof(argv[1]);
  float theta = atof(argv[2]);
  float beta = atof(argv[3]);

  const float mp = 938.272;
  
  float *jaco = new float[2];

  jaco = Jacobian(mp, Tinc, theta, beta);

  printf("%10.3f, %10.3f\n", jaco[0], jaco[1]);
  */
  float *output = new float[8]; // knockout output 
  
  float k = 220;
  float angk = 170;
  float angNN = 60;
  float BE    = 0;
  output = Knockout2D(16, 8,  200, k, angk, angNN, BE);

  printf("k:%6.2f angk:%6.2f angNN:%6.2f  Sp: %6.2f\n ---> T_c:%9.3f, theta_c:%9.3f, theta_d:%9.3f\n",
         k, angk, angNN, BE,  output[0], output[1], output[3]);    

  return 0;

}


