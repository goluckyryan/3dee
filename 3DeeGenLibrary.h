#include <iostream>
#include <string>
#include <stdio.h>
#include <cstring>

using namespace std;

double DWIA;
double PWIA;
double Ay;
double A00n0, Pn000, P0n00;

//***********************************************************************************
//***********************************************************************************
char* symbolL(int L){
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

char* symbolZ(int Z){
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
}

void orbit(int ID, int &N, int &L, float &J){
  switch(ID){
  case 1: N = 1; L = 0; J = 0.5; break;
  case 2: N = 1; L = 1; J = 1.5; break;
  case 3: N = 1; L = 1; J = 0.5; break;
  case 4: N = 1; L = 2; J = 2.5; break;
  case 5: N = 2; L = 0; J = 0.5; break;
  case 6: N = 1; L = 2; J = 1.5; break;
  case 7: N = 1; L = 3; J = 3.5; break;
  case 8: N = 2; L = 1; J = 1.5; break;
  }
}

int make_infile(int MA, int Z, float JA, float JB, float Ta, int N, int L, float J, float BE, float Tc, float theta_c, float theta_d){

  stringstream MAstr, ZAstr, Tastr, JAstr, JBstr, Nstr, Lstr, Jstr, BEstr, Tcstr, theta_cstr,theta_dstr;
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

  if(theta_c <0){
    printf(" ### theta_c should be positive. \n");
    return 1;
  }

  //read infile.temp___________________
  string line[19];
  ifstream file_in;
  file_in.open("infile.2p.temp");

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
  line1 << MA << symbolZ(Z) << "(p,2p)" << MA-1 << symbolZ(Z-1) <<"   " 
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


  //line17 Tc, theta_c, theta_d
  int lengthTc=(Tcstr.str()).length();
  int lengththeta_c=(theta_cstr.str()).length();
  int lengththeta_d=(theta_dstr.str()).length();
  line[17]=line[17].replace(0,10,"          "); //reset
  line[17]=line[17].replace(10,10,"          "); //reset
  line[17]=line[17].replace(20,10,"          "); //reset
  line[17]=line[17].replace(0,lengthTc,Tcstr.str());
  line[17]=line[17].replace(10,lengththeta_c,theta_cstr.str());
  line[17]=line[17].replace(20,lengththeta_d,theta_dstr.str());

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

  string line[84];
  ifstream file_in;
  file_in.open("outfile");

  if (!file_in) {
    printf("Unable to open outfile\n");
    exit(1);
  }
  
  for(int i=1; i<84; i++){
    getline(file_in,line[i]);
  }
  
  file_in.close();
  
  int lineAdj = 0;
  if ( line[linePWIA].length() == 0){
    lineAdj = 1;
  }else if (line[linePWIA+3].length() == 0){
    lineAdj = -1;
  }

  linePWIA += lineAdj;

  if ( line[linePWIA+6].length()==0 || line[linePWIA+9].length()==0 || line[linePWIA+10].length()==0) return 10;

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
   
  //printf("PWIA = %f, DWIA = %f, Ay = %f\n", PWIA, DWIA, Ay);

  return 0;
}

int AccpetanceFilter2D(float T1, float theta1, float T2, float theta2){

  // return 1 for accepted, 0 for rejected

  if(T1 < 30 || T2<30) return 0;

  if(theta1 < 20 || theta1 > 70) return 0;
  if(theta2 < 20 || theta2 > 70) return 0;

  return 1;

}

int AccpetanceFilter3D(float T1, float theta1, float phi1, float T2, float theta2, float phi2){

  // return 1 for accepted, 0 for rejected

  if(T1 < 30 || T2<30) return 0;

  if(theta1 < 20 || theta1 > 70) return 0;
  if(theta2 < 20 || theta2 > 70) return 0;

  return 1;

}

float mwdcY(float x){
  // for MWDC-L , x -> -x
  // for MWDC-R , x is normal

  if (x >= -170 && x < 200.*2./3.-170.){
    return 3.*(x+170.)/2.;
  }else if (x >= 200.*2./3.-170. && x < 200.*(-2./3.)+950.){
    return 200;
  }else if (x >= 200.*(-2./3.)+950. && x <= 950){
    return (x-950)*(-3./2.);
  }else{
    return 0;
  }
    

}


