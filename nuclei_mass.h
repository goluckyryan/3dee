#ifndef NUCLEIMASS_H
#define NUCLEIMASS_H

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>  //atoi

using namespace std;

string dataPath="mass03.txt";

double Nucleus_Mass(int Z, int A)
{
  string line;
  int    lineNum=0;
  int    list_A, list_Z;
  double list_mass, list_massDe;
  double mass;
  ifstream myfile;
  int    flag=0;

  const double mp = 938.272046;
  const double mn = 939.565378;

  myfile.open(dataPath.c_str());

  if (myfile.is_open())
    {
      while (/*! myfile.eof() &&*/ flag == 0 && lineNum <3218)
        {
          lineNum ++ ;
          getline (myfile,line);

          if (lineNum >= 40 ){
            list_Z     = atoi((line.substr(10,5)).c_str());
            list_A     = atoi((line.substr(15,5)).c_str());
            list_massDe = atof((line.substr(55,11)).c_str());
            list_mass  = list_Z*mp + (list_A-list_Z)*mn - list_massDe/1000*list_A;
            if ( A == list_A && Z == list_Z) {
              mass = list_mass;
              flag = 1;
            }

          }
        }
      myfile.close();
    }
  else cout << "Unable to open mass03.txt \n";

  if (flag == 1){
    return mass;
  }else{
  	return -404;
  }

  return 0;
}

string Nucleus_Name(int Z, int A){

  string line, name;
  int    lineNum=0;
  int    list_A, list_Z;
  ifstream myfile;
  int    flag=0;

  myfile.open(dataPath.c_str());

  if (myfile.is_open())
    {
      while ( flag == 0 && lineNum <3218)
        {
          lineNum ++ ;
          getline (myfile,line);

          if (lineNum >= 40 ){
            list_Z     = atoi((line.substr(10,5)).c_str());
            list_A     = atoi((line.substr(15,5)).c_str());
            name       = line.substr(19,3);
            if ( A == list_A && Z == list_Z) {
              flag = 1;
            }

          }
        }
      myfile.close();
    }
  else cout << "Unable to open mass03.txt \n";

  if (flag == 1){
    return name;
  }else{
    return "---";
  }

}

#endif
