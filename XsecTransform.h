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

    float *A = new float[5];
    A[0] = EnergyA;                   // E
    A[1] = MomtA*sin(thetaA*deg2rad); //px
    A[2] = MomtzA;                    //pz
    A[3] = MomtA;                     //|p|
    A[4] = thetaA;                    // angle

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
    //jaco[0] =  MomtA*MomtA*MomtA/Momt/Momt/gamma/(Momt + Energy*beta*cos(theta*deg2rad));   // d(sigma)/d(Omega)

    //The Journal of Chemical Physics 69, 1737 (1978)
    //jaco[0] = MomtA*MomtA*Energy/Momt/Momt/EnergyA;    // d^3(sigma)/dp/d(Omega)
    //jaco[1] = thetaA;

    //little modified // d^3(sigma)/dE/d(Omega)
    jaco[0] = Momt*Energy/MomtA/EnergyA;

    //
    jaco[1] = thetaA; //in deg

    return jaco;

}

float* Jacobian2(float T1, float theta1, float T2, float theta2, float beta){
    
    const float mass = 938.272; // proton

    theta2 = abs(theta2);

    float Energy1 = mass + T1;
    float Momt1 = sqrt(2*mass*T1+T1*T1);
    float Momt1x = Momt1*sin(theta1*deg2rad);
    float Momt1z = Momt1*cos(theta1*deg2rad);

    float Energy2 = mass + T2;
    float Momt2 = sqrt(2*mass*T2+T2*T2);
    float Momt2x = Momt2*sin(theta2*deg2rad);
    float Momt2z = Momt2*cos(theta2*deg2rad);

    float gamma = 1/sqrt(1-beta*beta);
    
    // Lorentz transform of A
    float *trans = new float[5];

    trans = LorentzTransform(mass, T1, theta1, beta);
    
    float EnergyA = trans[0];
    float MomtA   = trans[3];
    float thetaA  = trans[4];

    trans = LorentzTransform(mass, T2, theta2, beta);
    
    float EnergyB = trans[0];
    float MomtB   = trans[3];
    float thetaB  = trans[4];
    
    // Jacobian frome C fram to A frame;
    float *jaco = new float[2];

    float jaco1 = MomtA*EnergyA/Momt1/Energy1; // J(EA,OmegaA, E1, Omega1)

    float jaco2 = MomtB*MomtB*MomtB/Momt2/Momt2/gamma/(Momt2 + Energy2*beta*cos(theta2*deg2rad)); // J(OmegaB. Omega2)

    //jaco[0] =  MomtA*MomtA*MomtA/Momt/Momt/gamma/(Momt + Energy*beta*cos(theta*deg2rad));   // d(sigma)/d(Omega)

    //The Journal of Chemical Physics 69, 1737 (1978)
    //jaco[0] = MomtA*MomtA*Energy/Momt/Momt/EnergyA;    // d^3(sigma)/dp/d(Omega)
    //jaco[1] = thetaA;

    //little modified // d^3(sigma)/dE/d(Omega)
    jaco[0] = jaco1*jaco2;

    //
    jaco[1] = thetaA; //in deg

    return jaco;

}
