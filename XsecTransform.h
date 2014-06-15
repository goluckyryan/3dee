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

float* Jacobian2(float T_c, float theta_c, float T_d, float theta_d, float beta){
    // frome CM to Lab  
    const float mass = 938.272; // proton

    //T_c, theta_c  = Lab energy and angle
    //theta_d = Lab angle

    //A, B are CM frame

    theta_d = abs(theta_d);

    float Energy_c = mass + T_c;
    float Momt_c = sqrt(2*mass*T_c+T_c*T_c);
    float Momt_cx = Momt_c*sin(theta_c*deg2rad);
    float Momt_cz = Momt_c*cos(theta_c*deg2rad);

    float Energy_d = mass + T_d;
    float Momt_d = sqrt(2*mass*T_d+T_d*T_d);
    float Momt_dx = Momt_d*sin(theta_d*deg2rad);
    float Momt_dz = Momt_d*cos(theta_d*deg2rad);

    float gamma = 1/sqrt(1-beta*beta);
    
    // Lorentz transform to Lab
    float *trans = new float[5];

    trans = LorentzTransform(mass, T_c, theta_c, beta);
    
    float Energy_1 = trans[0];
    float Momt_1   = trans[3];
    float theta_1  = trans[4];

    trans = LorentzTransform(mass, T_d, theta_d, beta);
    
    float Energy_2 = trans[0];
    float Momt_2   = trans[3];
    float theta_2  = trans[4];
    
    // Jacobian frome CM fram to Lab frame;
    float *jaco = new float[2];

    float jaco1 = Momt_1*Energy_1/Momt_c/Energy_c; // J(E_c,Omega_c, E_1, Omega_1)

    float jaco2 = Momt_2*Momt_2*Momt_2/Momt_d/Momt_d/gamma/abs(Momt_d + Energy_d*beta*cos(theta_d*deg2rad)); // J(Omega_d. Omega_2)

    jaco[0] = jaco1*jaco2;

    //axuillary
    jaco[1] = theta_1; //in deg 

    return jaco;

}

float* Jacobian3(const float *output){
    // frome CM to Lab  
    const float mass = 938.272; // proton

    // 0 = T_c, 1 = theta_c, 2 = T_d, 3 = theta_d //CM frame
    // 4 , 5, 6, 7, // Lab frame
    // 8 =  beta , from CM to Lab

    float Energy_c = mass + output[0];
    float Momt_c = sqrt(2*mass*output[0]+output[0]*output[0]);

    float Energy_d = mass + output[2];
    float Momt_d = sqrt(2*mass*output[2]+output[2]*output[2]);
    float theta_d = abs(output[3]);

    float beta = output[8];
    float gamma = 1/sqrt(1-output[8]*output[8]);
    
    float Energy_1 = mass + output[4];
    float Momt_1   = sqrt(2*mass*output[4]+output[4]*output[4]);
    
    float Momt_2   = sqrt(2*mass*output[5]+output[5]*output[5]);
    
    // Jacobian frome CM fram to Lab frame;
    float *jaco = new float[2];

    float jaco1 = Momt_1*Energy_1/Momt_c/Energy_c; // J(E_c,Omega_c, E_1, Omega_1)

    float jaco2 = Momt_2*Momt_2*Momt_2/Momt_d/Momt_d/gamma/abs(Momt_d + Energy_d*beta*cos(theta_d*deg2rad)); // J(Omega_d. Omega_2)

    jaco[0] = jaco1*jaco2;

    //axuillary, not use
    jaco[1] = jaco1; 

    return jaco;

}
