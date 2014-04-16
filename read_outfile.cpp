#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

using namespace std;


int main(){

    string x,str;
    string xsecstr,PWIAstr,DWIAstr,A00n0str;
    double xsec,PWIA,DWIA,A00n0;

    ifstream file_in;
    file_in.open("outfile");

    if (!file_in) {
        cerr << "Unable to open file outfile";
        exit(1);   // call system to stop
    }

    cout << "outfile opened." << endl;

    while(file_in >> x){
        str = str+x;
    }

    //cout << str <<endl;

    size_t found;
    found=str.find("(p,2p)");
    string react;
    if (found!=string::npos){
        react=str.substr(found-3,12);
//        cout << "         reaction = " << react << endl; 
    }
    string norb;
    if (found!=string::npos){
        norb=str.substr(found+9,5);
//        cout << " orbit of knockout proton = " << norb << endl;
        
    }

    found=str.find("INCIDENT");
    double InKE;
    if (found!=string::npos){
        InKE=atof((str.substr(found+15,7)).c_str());
//        cout << "  Incident energy = " << InKE << " [MeV]" <<endl;
    }

    found=str.find("BINDING");
    double BE;
    if (found!=string::npos){
        BE=atof((str.substr(found+14,7)).c_str());
//        cout << "   Binding Energy = " << BE<< " [MeV]" <<endl;
    }
    
    found=str.find("ENERGY");
    found=str.find("ENERGY", found+10);
    found=str.find("ENERGY", found+10);
    double engdect;
    if (found!=string::npos){
        engdect=atof((str.substr(found+7,7)).c_str());
//        cout << "detected energy = " << engdect<< " [MeV]" << endl;
    }

    found=str.find("ANGLE1");
    double ang1;
    if (found!=string::npos){
        ang1=atof((str.substr(found+7,7)).c_str());
//        cout << "angle of knockout = " << ang1<< " [deg]" << endl;
    }

    found=str.find("ANGLE2");
    double ang2;
    if (found!=string::npos){
        ang2=atof((str.substr(found+7,7)).c_str());
//        cout << " angle of scatter = " << ang2<< " [deg]" << endl;
    }

    found=str.find("Cross");
    if (found!=string::npos){
        xsecstr=str.substr(found+10,10);
        xsec= atof(xsecstr.c_str()) ;
//        cout << "    Cross Section = " << xsec  << " [mb]" << endl;
    }

    found=str.find("PWIA");
    if (found!=string::npos){
        PWIAstr=str.substr(found+11,10);
        PWIA= atof(PWIAstr.c_str()) ;
//        cout << "   PWIA x-Section = " << PWIA  << " [mb]" << endl;
    }

    found=str.find("DWIA");
    if (found!=string::npos){
        DWIAstr=str.substr(found+11,10);
        DWIA= atof(DWIAstr.c_str()) ;
//        cout << "   DWIA x-Section = " << DWIA  << " [mb]" << endl;
    }

    found=str.find("A00n0");
    if (found!=string::npos){
        A00n0str=str.substr(found+5,10);
        A00n0= atof(A00n0str.c_str()) ;
//        cout << "       PWIA A00n0 = " << A00n0  << endl << endl;
    }

	//cout << endl << endl << str ;

    file_in.close();
    ofstream fout("test_readout.txt", ios::app);
    fout << react << "\t" << norb <<"\t" << BE << "\t" << engdect << "\t" << ang1 << "\t" << ang2 << "\t";
    fout << xsec  << "\t" << PWIA <<"\t" <<DWIA << "\t" <<A00n0 << endl;
    fout.close();
    return 0;

}
