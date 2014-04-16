#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <stdlib.h>
#include <cmath>
#include <cstring>

using namespace std;

string reaction="F(p,2p)", residual="O", reactance="Fluorine";

int main(int argc, char *argv[]){

     if(argc <=9) {
 	cout << "=======================================================" << endl
 	     << "          Only for "<< reaction << residual<<" knockout [A(a,cd)b]" << endl << endl
	     << "Usage: ./make_infile.o\e[32m MA T_inc N L J BE Tc theta_c theta_d \e[m" << endl << endl
//	     << "     ZA = Charge of isotop (X.)" << endl
	     << "     MA = Mass number of isotop (X.)" << endl
	     << "  T_inc = incident energy proton (X.X)" << endl
	     << "  N L J = shell number for orbital proton (X X X.X)" << endl
	     << "     BE = binding energy of orbital proton (X.X)" << endl
	     << "     Tc = K.E. of knockouted neutron (X.X)" << endl
	     << "theta_c = detection angle of knockout neutron (positive angle) (X.X)" << endl
	     << "theta_d = detection angle of scattered proton (X.X)"<<endl
	     << "  **** for (p,pn) or detail, change code line62"<< endl << endl;
	exit(0);
    }


    int MA = atoi(argv[1]);
    float T_inc = atof(argv[2]);
    int N = atoi(argv[3]);
    int L = atoi(argv[4]);
    float J = atof(argv[5]);
    float BE = atof(argv[6]);
    float Tc = atof(argv[7]);
    float theta_c =atof(argv[8]);
    float theta_d = atof(argv[9]);

    if(theta_c <0){
        cout << " ### theta_c should be positive. " << endl;
        exit(0);
    }
    if(MA <14 || MA>24 ){
        cout << " only "<< reactance <<" 14 to 24 !! " << endl;
        exit(0);
    }
    if((10*MA)%10 !=0 ) {
        cout << " Mass number (int) is ok." << endl;
        exit(0);
    }
    if((10*N)%10 !=0 || (10*L)%10 != 0 || fmod(J*10,5) !=0 || abs(L-J) != 0.5 ){
        cout << "N L J is not well defined" << endl;
        exit(0);
    }


    //read infile.XX.temp___________________
    string line[19];
    ifstream file_in;
    file_in.open("infile.2p.temp");

    if (!file_in) {
        cerr << "Unable to open file infile.2p.temp";
        exit(1);   // call system to stop
    }

    for(int i=1; i<=18; i++){
        getline(file_in,line[i]);
        //cout << "line " << i << "\t: " << line[i] << endl;
    }
    file_in.close();

    //modify___________________________
    string orb;

    switch (L){
    	case 0: orb="s";break;
    	case 1: orb="p";break;
    	case 2: orb="d";break;
    	case 3: orb="f";break;
    	case 4: orb="g";break;
    	case 5: orb="h";break;
    	case 6: orb="i";break;

    }

    std::ostringstream line1;
    line1 << MA << reaction << MA-1 << residual <<"    " << N << orb << J << "   Relativitic, T_inc =" << argv[2] ;
    line[1]= line1.str();

    // string.replace( start pos, number of pos, string )
    line[4]=line[4].replace(10,1,argv[3]); //N
    line[4]=line[4].replace(20,1,argv[4]); //L
    line[4]=line[4].replace(30,3,argv[5]); //J

    line[5]=line[5].replace(0,2,argv[1]);
    int lengthTa=strlen(argv[2]);
    if (lengthTa >10){
        cout << "T_inc is larger then 8 digit." << endl;
        exit(4);
    }
    line[5]=line[5].replace(40,10,"          "); //reset;
    line[5]=line[5].replace(40,lengthTa,argv[2]);

    std::stringstream MBstr;
    MBstr << MA-1;
    line[8]=line[8].replace(0,2,MBstr.str()) ;


    //MA , JA, JB
    string JA,JB; // JX = spin of nucleus X

/* Oxygen spin 
    switch(MA){
    	case 14: JA="0. "; JB="0.5";break;
    	case 15: JA="1.5"; JB="0. ";break;
    	case 16: JA="0. "; JB="1.5";break;
        case 17: JA="2.5"; JB="0. ";break;
        case 18: JA="0. "; JB="2.5";break;
        case 19: JA="2.5"; JB="0. ";break;
        case 20: JA="0. "; JB="2.5";break;
        case 21: JA="2.5"; JB="0. ";break;
        case 22: JA="0. "; JB="0.5";break;
        case 23: JA="1.5"; JB="0. ";break;
        case 24: JA="0. "; JB="1.5";break;
    }
*/

//Fluorine
    switch(MA){
    	case 14: JA="2. "; JB="1.5";break;
    	case 15: JA="1.5"; JB="0. ";break;
    	case 16: JA="0. "; JB="1.5";break;
        case 17: JA="2.5"; JB="0. ";break;
        case 18: JA="1. "; JB="2.5";break;
        case 19: JA="1.5"; JB="0. ";break;
        case 20: JA="2. "; JB="2.5";break;
        case 21: JA="2.5"; JB="0. ";break;
        case 22: JA="4. "; JB="0.5";break;
        case 23: JA="2.5"; JB="0. ";break;
        case 24: JA="1. "; JB="1.5";break;
    }
    line[3]=line[3].replace(0,3,JA);
    line[3]=line[3].replace(10,3,JB);


    //line14 BE
    int lengthBE=strlen(argv[6]);
    line[14]=line[14].replace(0,10,"          "); //reset
    line[14]= line[14].replace(40,lengthBE,argv[6]);


    //line17 Tc, theta_c, theta_d
    int lengthTc=strlen(argv[7]);
    int lengththeta_c=strlen(argv[8]);
    int lengththeta_d=strlen(argv[9]);
    line[17]=line[17].replace(0,10,"          "); //reset
    line[17]=line[17].replace(10,10,"          "); //reset
    line[17]=line[17].replace(20,10,"          "); //reset
    line[17]=line[17].replace(0,lengthTc,argv[7]);
    line[17]=line[17].replace(10,lengththeta_c,argv[8]);
    line[17]=line[17].replace(20,lengththeta_d,argv[9]);

    cout << "========write -> infile ==========" << endl;
    for (int i=1; i <=18 ; i++){
    	cout << "line " << i << "\t: " << line[i] << endl;
    }


    ofstream file_out;
    file_out.open("infile");

    for(int i=1; i<=18; i++){
        file_out << line[i] << endl;
        //cout << "line " << i << " : " << line[i] << endl;
    }
    file_out.close();


    return 0;

}
