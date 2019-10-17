#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


int main(int argc, char **argv) { //Down and out
    if(argc != 8){
        printf("Usage: ./call <K> <S_0> <T> <r> <q> <sigma> <H>");
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
    double H = atof(argv[7]);

    int N = 20000;  //Amount of Points in discrete time
    double* p = malloc(sizeof(double)*N);
    double delta_t = T / (double) N;
    double u = exp(sigma* sqrt(delta_t));
    double d = 1/u;
    double rho = exp((r-q)*delta_t);
    double pi= (rho-d)/(u-d);
    double p0 = (u*exp(-q * delta_t) - exp(-r * delta_t)) / (u*u - 1);
    double p1 = exp(-r*delta_t) - p0;
    int m = -round(-log(H/S_0)/log(u));
    int n_hat;
    int* lambda = malloc(sizeof(int)*(N+1));
    for(int i =0;i<=N;i++){
        lambda[i] = -round(-(m+i)/2);
    }
    int m_abs = m;
    if(m_abs<0) {
        m_abs=-m_abs;
    }
    n_hat = N - m_abs;
    double** Q = malloc(sizeof(double)*(n_hat+1));
    for(int i = 0;i<=n_hat;i++){
        Q[i]= malloc(sizeof(double)*i);
    }
    double** P = malloc(sizeof(double)*(N+1));
    for(int i = 0;i<=N;i++){
        P[i]= malloc(sizeof(double)*i);
    }
    for(int j = 0;j<=n_hat;j++){
        Q[n_hat][j] = K - S_0*pow(u,m) * pow(u,(2*j - n_hat));
        if(Q[n_hat][j]< 0){
            Q[n_hat][j]= 0;
        }
    }

    for(int i = n_hat-1;i>=0;i--){
        for(int j =0;j<=i;j++){
            Q[i][j]=(p0 * Q[i+1][j+1] + p1 * Q[i+1][j]);
            double exercise = K - S_0*pow(u,m) * pow(u,2*j - i);
            if(Q[i][j] < exercise){
                Q[i][j] = exercise;
            }
        }
    }
    for(int j = 0;j<=lambda[N]-1;j++){
        P[N][j]=0;
    }
    P[N][lambda[N]]=Q[N-m_abs][(int) -round(-(N-m_abs)/2)];

    for(int i = N-1;i>=m_abs;i--){
        for(int j =0;j<=lambda[N]-1; j++){
            P[i][j]=(p0 * P[i+1][j+1] + (p1) * P[i+1][j]);
        }
        P[i][lambda[i]]=Q[i-m_abs][(int) -round(-(i-m_abs)/2)];
    }
    for(int i =m_abs-1;i>=0;i--){
        for(int j=0;j<=i;j++){
            P[i][j]=(p0 * P[i+1][j+1] + (p1) * P[i+1][j]);
        }
    }

    clock_t end = clock();
    float seconds = (float)(end - start) / CLOCKS_PER_SEC;
//   printf("Price = %lf, calculated in %f seconds", p[0], seconds);
//    printf("Price = %lf, calculated in %f seconds", p[0], seconds);
    printf("%f, %f, %f, %f, %f, %f,%f, %f, %f\n", K,S_0,T,r,q,sigma,H, P[0][0], seconds);

    return 0;
}
