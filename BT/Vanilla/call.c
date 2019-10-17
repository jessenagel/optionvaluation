#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


int main(int argc, char **argv) {
    if(argc != 7){
        printf("Usage: ./call <K> <S_0> <T> <r> <q> <sigma>");
        exit(1);
    }
    srand(time(0));//Set seed for random draws
    clock_t start = clock();
    double K = atof(argv[1]);
    double S_0 = atof(argv[2]);
    double T = atof(argv[3]);
    double r = atof(argv[4]);
    double q = atof(argv[5]);
    double sigma = atof(argv[6]);
    int N = 20000;  //Amount of Points in discrete time
    double* p = malloc(sizeof(double)*N);
    double delta_t = T / (double) N;
    double up = exp(sigma* sqrt(delta_t));
    double p0 = (up*exp(-q * delta_t) - exp(-r * delta_t)) / (up*up - 1);
    double p1 = exp(-r*delta_t) - p0;
    for(int i = 0;i<N;i++){
        p[i] = S_0 * pow(up,(2*i - N)) - K;
        if(p[i]< 0){
            p[i]= 0;
        }
    }
    for(int j = N-2;j>=0;j--){
        for(int i =0;i<=j; i++){
            p[i]=p0 * p[i+1] + p1 * p[i];
            double exercise =  S_0 * pow(up,2*i - j) -K ;
            if(p[i] < exercise){
                p[i] = exercise;
            }
        }
    }

    clock_t end = clock();
    float seconds = (float)(end - start) / CLOCKS_PER_SEC;
//    printf("Price = %lf, calculated in %f seconds", p[0], seconds);
    printf("%f, %f, %f, %f, %f, %f, %f, %f\n", K,S_0,T,r,q,sigma, p[0], seconds);

    return 0;
}
