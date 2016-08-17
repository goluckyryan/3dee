#ifndef THREEDEELibrary
#define THREEDEELibrary

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <cstring>
#include "nuclei_mass.h"

using namespace std;

double DWIA;
double PWIA;
double Ay;
double A00n0, Pn000, P0n00;

//***********************************************************************************
const char* symbolL(int L);
const char* symbolZ(int Z, int A);
void orbit(int ID, int &N, int &L, float &J);
int make_infile(string filename, int MA, int Z, float JA, float JB, float Ta, int N, int L, float J, float BE, float Tc, float theta_c, float theta_d, float ang_d);
int read_outfile(int linePWIA);
int AccpetanceFilter2D(float T1, float theta1, float T2, float theta2);
int AccpetanceFilter3D(float T1, float theta1, float phi1, float T2, float theta2, float phi2);
float mwdcY(float x);
float* RotY( float * V, float ang );

//***********************************************************************************
const char* symbolL(int L){
  switch (L){
  case 0: return "s";break;
  case 1: return "p";break;
  case 2: return "d";break;
  case 3: return "f";break;
  case 4: return "g";break;
  case 5: return "h";break;
  case 6: return "i";break;
  }
}

const char* symbolZ(int Z, int A){
/*
  switch (Z){
  case 1: return "H";break;
  case 2: return "He";break;
  case 3: return "Li";break;
  case 4: return "Be";break;
  case 5: return "B";break;
  case 6: return "C";break;
  case 7: return "N";break;
  case 8: return "O";break;
  case 9: return "F";break;
  }
*/

  string sym = Nucleus_Name(Z,A);
  return sym.c_str();
  
}

void orbit(int ID, int &N, int &L, float &J){
  switch(ID){
  case  1: N = 1; L = 0; J = 0.5; break; // 1s1/2
  // 2
  case  2: N = 1; L = 1; J = 1.5; break; // 1p3/2
  case  3: N = 1; L = 1; J = 0.5; break; // 1p1/2
  // 8
  case  4: N = 1; L = 2; J = 2.5; break; // 1d5/2
  case  5: N = 2; L = 0; J = 0.5; break; // 2s1/2
  case  6: N = 1; L = 2; J = 1.5; break; // 1d3/2
  // 20
  case  7: N = 1; L = 3; J = 3.5; break; // 1f7/2
  // 28
  case  8: N = 2; L = 1; J = 1.5; break; // 2p3/2
  case  9: N = 1; L = 3; J = 2.5; break; // 1f5/2
  case 10: N = 2; L = 1; J = 0.5; break; // 2p1/2
  // 40
  case 11: N = 1; L = 4; J = 4.5; break; // 1g9/2
  // 50
  case 12: N = 2; L = 2; J = 2.5; break; // 2d5/2
  case 13: N = 1; L = 4; J = 3.5; break; // 1g7/2
  case 14: N = 3; L = 1; J = 0.5; break; // 3s1/2
  case 15: N = 2; L = 2; J = 1.5; break; // 2d3/2
  case 16: N = 1; L = 5; J = 5.5; break; // 1h11/2
  //82
  case 17: N = 2; L = 3; J = 3.5; break; // 2f7/2
  case 18: N = 1; L = 5; J = 4.5; break; // 1h9/2
  case 19: N = 1; L = 6; J = 6.5; break; // 1i13/2
  case 20: N = 3; L = 1; J = 1.5; break; // 3p3/2
  case 21: N = 2; L = 3; J = 2.5; break; // 2f5/2
  case 22: N = 3; L = 1; J = 0.5; break; // 3p1/2
  //126
  }
}

int make_infile(string filename, int MA, int Z, float JA, float JB, float Ta, int N, int L, float J, float BE, float Tc, float theta_c, float theta_d, float ang_d){

  stringstream MAstr, ZAstr, Tastr, JAstr, JBstr, Nstr, Lstr, Jstr, BEstr, Tcstr, theta_cstr,theta_dstr, ang_dstr;
  MAstr << MA;
  ZAstr<<Z;
  Tastr<<Ta;
  JAstr<<JA;
  JBstr<<JB;
  Nstr << N;
  Lstr <<L;
  Jstr <<J;
  BEstr <<BE;
  Tcstr<<Tc;
  theta_cstr<<theta_c;
  theta_dstr<<theta_d;
  ang_dstr<<ang_d;

  if(theta_c <0){
    printf(" ### theta_c should be positive. \n");
    return 1;
  }

  //read infile.temp___________________
  string line[19];
  ifstream file_in;
  //file_in.open("infile.2p.temp");
  file_in.open(filename.c_str());

  if (!file_in) {
    cerr << "Unable to open file infile";
    return 0;   // call system to stop
  }

  for(int i=1; i<=18; i++){
    getline(file_in,line[i]);
  }
  file_in.close();

  //modify___________________________

  std::ostringstream line1;
  line1 << MA << symbolZ(Z, MA) << "(p,2p)" << MA-1 << symbolZ(Z-1, MA-1) <<"   " 
        << N << symbolL(L) << (int)2*J << "/2   Relativitic, Ta =" << Ta << " MeV";
  line[1]= line1.str();

  // string.replace( start pos, number of pos, string )
  line[4]=line[4].replace(10,1,Nstr.str());
  line[4]=line[4].replace(20,1,Lstr.str());
  line[4]=line[4].replace(30,3,Jstr.str());

  line[5]=line[5].replace(0,2,MAstr.str());
  line[5]=line[5].replace(10,1,ZAstr.str());
  int lengthTa=(Tastr.str()).length();
  if (lengthTa >10){
    printf( "Ta is larger then 8 digit. \n");
    return 0;
  }
  line[5]=line[5].replace(40,10,"          "); //reset;
  line[5]=line[5].replace(40,lengthTa,Tastr.str());

  std::stringstream MBstr;
  MBstr << MA-1;
  line[8]=line[8].replace(0,2,MBstr.str()) ;
  std::stringstream ZBstr;
  ZBstr << Z-1;
  line[8]=line[8].replace(10,1,ZBstr.str()) ;
  

  //JA, JB
  int lengthJAstr = (JAstr.str()).length();
  int lengthJBstr = (JBstr.str()).length();
  line[3]=line[3].replace(0,lengthJAstr,JAstr.str());
  line[3]=line[3].replace(10,lengthJBstr,JBstr.str());


  //line14 BE
  int lengthBE=(BEstr.str()).length();
  line[14]=line[14].replace(0,10,"          "); //reset
  line[14]= line[14].replace(40,lengthBE,BEstr.str());


  //line17 Tc, theta_c, theta_d, ang_d
  int lengthTc=(Tcstr.str()).length();
  int lengththeta_c=(theta_cstr.str()).length();
  int lengththeta_d=(theta_dstr.str()).length();
  int lengthang_d=(ang_dstr.str()).length();
  line[17]=line[17].replace(0,10,"          "); //reset
  line[17]=line[17].replace(10,10,"          "); //reset
  line[17]=line[17].replace(20,10,"          "); //reset
  line[17]=line[17].replace(30,10,"          "); //reset
  line[17]=line[17].replace(0,lengthTc,Tcstr.str());
  line[17]=line[17].replace(10,lengththeta_c,theta_cstr.str());
  line[17]=line[17].replace(20,lengththeta_d,theta_dstr.str());
  line[17]=line[17].replace(30,lengthang_d, ang_dstr.str());

  // Save to infile
  ofstream file_out;
  file_out.open("infile");

  for(int i=1; i<=18; i++){
    file_out << line[i] << endl;
  }
  file_out.close();

  return 1;

}

int read_outfile(int linePWIA){
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
    if( file_in.eof() ) break;
  }
  
  file_in.close();

//  for( int i = 0; i< 10; i++){
//    printf("%s\n", line[linePWIA+i].c_str());
//  }

  int lineAdj = 0;
  if ( line[linePWIA].length() == 0){
    lineAdj = 1;
  }else if (line[linePWIA+3].length() == 0){
    lineAdj = -1;
  }

  linePWIA += lineAdj;

//  printf("----------1\n");

  if ( line[linePWIA+6].length()==0 || line[linePWIA+9].length()==0 || line[linePWIA+10].length()==0)  return 10;

  /* printf("%2d | %2d |%s \n", linePWIA   , line[linePWIA   ].length(), line[linePWIA].c_str());
     printf("%2d | %2d |%s \n", linePWIA+3 , line[linePWIA+3 ].length(), line[linePWIA+3].c_str());
     printf("%2d | %2d |%s \n", linePWIA+6 , line[linePWIA+6 ].length(), line[linePWIA+6].c_str());
     printf("%2d | %2d |%s \n", linePWIA+9 , line[linePWIA+9 ].length(), line[linePWIA+9].c_str());
     printf("%2d | %2d |%s \n", linePWIA+10, line[linePWIA+10].length(), line[linePWIA+10].c_str());
  */
  PWIA  = atof((line[linePWIA]).substr(23,15).c_str());
  DWIA  = atof((line[linePWIA+3]).substr(23,15).c_str());
  A00n0 = atof((line[linePWIA+6]).substr(40,15).c_str());
  Pn000 = atof((line[linePWIA+9]).substr(40,15).c_str());
  P0n00 = atof((line[linePWIA+10]).substr(40,15).c_str());

  Ay = Pn000;
   
  printf("PWIA = %f, DWIA = %f, Ay = %f\n", PWIA, DWIA, Ay);

  return 0;
}

int AccpetanceFilter2D(float T1, float theta1, float T2, float theta2){

  // return 1 for accepted, 0 for rejected

  //printf("------------------  T1:%9.3f ang1:%9.3f T2:%9.3f ang2:%8.3f \n", T1, theta1, T2, theta2);

  if(T1 < 14 || T2<14 || T1>350 || T2 > 350) return 0;

  if(theta1 < 17 || theta1 > 70) return 0;

  if(theta2 < 17 || theta2 > 70) return 0;

  return 1;

}

int AccpetanceFilter3D(float T1, float theta1, float phi1, float T2, float theta2, float phi2){

  // return 1 for accepted, 0 for rejected
  bool debug = 0;

  if(T1 < 14 || T2<14 || T1>400 || T2 > 400) return 0;
  theta1 = theta1/rad2deg;
  phi1   = phi1/rad2deg;
  theta2 = theta2/rad2deg;
  phi2   = phi2/rad2deg;

  const float mwdcDis = 1022.37;
  const float mwdcAng = 60/rad2deg;

  // the vector to MWDC
  float len1;
  len1 = mwdcDis/(sin(theta1)*cos(phi1)*sin(mwdcAng)+cos(theta1)*cos(mwdcAng)); 
  float v1[3];
  v1[0] = len1*cos(phi1)*sin(theta1);
  v1[1] = len1*sin(phi1)*sin(theta1);
  v1[2] = len1*cos(theta1);

  if(debug)printf("len1:%10.3f, v1: {%10.3f, %10.3f, %10.3f}\n", len1, v1[0], v1[1], v1[2]);

  float len2;
  len2 = mwdcDis/(sin(theta2)*cos(phi2)*sin(-mwdcAng)+cos(theta2)*cos(mwdcAng)); 
  float v2[3];
  v2[0] = len2*cos(phi2)*sin(theta2);
  v2[1] = len2*sin(phi2)*sin(theta2);
  v2[2] = len2*cos(theta2);

  if(debug)printf("len2:%10.3f, v2: {%10.3f, %10.3f, %10.3f}\n", len2, v2[0], v2[1], v2[2]);
  
  float* vL = new float [3];
  vL = RotY(v1, 60);
  if(debug)printf("vL: {%10.3f, %10.3f, %10.3f}\n", vL[0], vL[1], vL[2]);
  if(debug)printf("  mwdcY(%10.3f)=%10.3f , Y %10.3f\n",vL[0], mwdcY(-vL[0]), vL[1]);

  float* vR = new float[3];
  vR = RotY(v2, -60);
  if(debug)printf("vR: {%10.3f, %10.3f, %10.3f}\n", vR[0], vR[1], vR[2]);
  if(debug)printf("  mwdcY(%10.3f)=%10.3f , Y %10.3f\n", vR[0], mwdcY(vR[0]), vR[1]);
  
  if ( mwdcY(-vL[0]) < abs(vL[1]) ) return 0; 
  if ( mwdcY(vR[0]) < abs(vR[1]) ) return 0; 

  //printf("...... accpected ..... \n");

  return 1;

}

float mwdcY(float x){
  // for MWDC-L , x -> -x
  // for MWDC-R , x is normal

  if (x >= -170 && x < 200.*2./3.-170.){ // (-170, -36.67), D= 133.333
    return 3.*(x+170.)/2.;
  }else if (x >= 200.*2./3.-170. && x < 200.*(-2./3.)+950.){ //(-36.67,816.67), D = 853.34
    return 200;
  }else if (x >= 200.*(-2./3.)+950. && x <= 950){ //(816.67, 950), D = 133.333
    return (x-950)*(-3./2.);
  }else{
    return -1;
  }

}

float* RotY( float * V, float ang ){
  // Left Hand rotate in 3D
  float *U = new float[3];
  ang = ang/rad2deg;  

  U[0] = V[0]*cos(ang) - V[2]*sin(ang);
  U[1] = V[1];
  U[2] = V[0]*sin(ang) + V[2]*cos(ang);

  return U;
}


#endif
