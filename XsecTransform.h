#include <cmath> // math , enable simple math function
#include <stdlib.h> // standard General Utlilities Library

using namespace std;

const float deg2rad = 3.14159/180;

float* LorentzTransform(float mass, float Tinc, float theta, float beta){

    float Energy = mass + Tinc;
    float Momt = sqrt(2*mass*Tinc+Tinc*Tinc);
    float Momtx = Momt*sin(theta*deg2rad);
    float Momtz = Momt*cos(theta*deg2rad);
    float gamma = 1/sqrt(1-beta*beta);
    
    // Lorentz transform 
    float EnergyA = gamma*Energy + gamma*beta*Momtz;
    float MomtzA  = gamma*beta*Energy + gamma*Momtz;
    float MomtA   = sqrt(MomtzA*MomtzA + Momtx*Momtx);
    float thetaA  = acos(MomtzA/MomtA)/deg2rad;

    float *A = new float[3];
    A[0]=EnergyA;
    A[1]=MomtA*sin(thetaA*deg2rad);
    A[2]=MomtzA;

    return  A;

}

float* Jacobian(float mass, float Tinc, float theta, float beta){
    

    float Energy = mass + Tinc;
    float Momt = sqrt(2*mass*Tinc+Tinc*Tinc);
    float Momtx = Momt*sin(theta*deg2rad);
    float Momtz = Momt*cos(theta*deg2rad);
    float gamma = 1/sqrt(1-beta*beta);
    
    // Lorentz transform 
    float EnergyA = gamma*Energy + gamma*beta*Momtz;
    float MomtzA  = gamma*beta*Energy + gamma*Momtz;
    float MomtA   = sqrt(MomtzA*MomtzA + Momtx*Momtx);
    float thetaA  = acos(MomtzA/MomtA)/deg2rad;
    
    // Jacobian frome C fram to A frame;
    float *jaco = new float[2];
    jaco[0] =  MomtA*MomtA*MomtA/Momt/Momt/gamma/(Momt + Energy*beta*cos(theta*deg2rad));   // d(sigma)/d(Omega)

    //The Journal of Chemical Physics 69, 1737 (1978)
    //jaco[0] = MomtA*MomtA*Energy/Momt/Momt/EnergyA;    // d^3(sigma)/dv/d(Omega)
    jaco[1] = thetaA;

    return jaco;

}
