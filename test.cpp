#include <string>
#include <sstream>
#include <cstdlib>
#include <stdlib.h>
#include <cmath>
#include <cstring>
#include "XsecTransform.h"

using namespace std;

int main(int argc, char *argv[]){

    float Tinc = atof(argv[1]);
    float theta = atof(argv[2]);
    float beta = atof(argv[3]);

    const float mp = 938.272;

    float *jaco = new float[2];

    jaco = Jacobian(mp, Tinc, theta, beta);

    printf("%10.3f, %10.3f\n", jaco[0], jaco[1]);

    return 0;

}


