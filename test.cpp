#include <istream> 
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
#include "3DeeGenLibrary.h"
//#include "XsecTransform.h"
#include "knockout2D.h"
#include "nuclei_mass.h"
//#include "knockout3D.h"

using namespace std;

int main(int argc, char *argv[]){
  const int nline = 100;
  string line[nline];
  ifstream file_in;
  file_in.open("outfile");
  
  if (!file_in) {
    printf("Unable to open outfile\n");
    exit(1);
  }

  for(int i=0; i<nline; i++){
    getline(file_in,line[i]);

    size_t pos = line[i].find("PWIA x-sec mb"); //get the last    
    if( pos != std::string::npos){
      double PWIA = atof(line[i].substr(pos+13, 16).c_str());
      printf(" PWIA : %.10f \n", PWIA);
    }

    pos = line[i].find("DWIA x-sec mb"); //get the last    
    if( pos != std::string::npos){
      double DWIA = atof(line[i].substr(pos+13, 16).c_str());
      printf(" DWIA : %.10f \n", DWIA);
    }

    pos = line[i].find("Beam a A00n0"); //get the last    
    if( pos != std::string::npos){
      double A00n0 = atof(line[i].substr(pos+28, 16).c_str());
      printf(" A00n0: %.10f \n", A00n0);
    }


    if( file_in.eof() ) break;
  }

  

  
  file_in.close();

  /*
  float *output = new float[13]; // knockout output 
  float *outputINV = new float[13]; // knockout inver output
  
  float TKA = 289.44;
  int Ma = 23;
  int Z  = 9;
  float k = 100;
  float angk = 60;
  float phik = 14.528;
  float angNN = 70;
  float phiNN = -18.994;
  float BE    = 13.26;

  cout << "---"<< Nucleus_Name(9,23) << "---"<<endl;

  printf("%d---%3s----\n",23, Nucleus_Name(9,23).c_str());
  printf("---%s----\n", symbolZ(9,23));
  printf("---%s----\n", symbolZ(20,40));
  printf("---%s----\n", symbolZ(50,120));

  //output = Knockout2D(Ma, Z, TKA, k, angk, angNN, BE);
  
  //outputINV = Knockout2Dinv(Ma, Z, TKA, output[4], output[5], output[6], output[7]);
  //  outputINV = Knockout2Dinv3(Ma, Z, TKA, output[0], output[1], output[3], BE);
  //outputINV = Knockout3Dinv3(Ma, Z, TKA, 90, 80, 70, 2,  BE);

/*
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

  float theta1 = atof(argv[1]);
  float phi1   = atof(argv[2]);
  float theta2 = atof(argv[3]);
  float phi2   = atof(argv[4]);
  

  printf("T1:%7.1f, theta1: %7.1f, phi1: %7.1f\n", outputINV[0], outputINV[1], outputINV[2]);
  printf("T2:%7.1f, theta2: %7.1f, phi2: %7.1f\n", outputINV[3], outputINV[4], outputINV[5]);
  
  printf("--- Accpetance Filter --- %d \n",AccpetanceFilter3D(outputINV[0], outputINV[1], outputINV[2], outputINV[3], outputINV[4], outputINV[5]));
  */
  return 0;

}


