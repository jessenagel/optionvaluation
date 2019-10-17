#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
double chebyshev(int n, double x,double* c){

    double* w=malloc(sizeof(double)*(n+2));
    w[n+1]=0;
    w[n]=0;
    for(int j=n-1 ;j>=0;j--){
        w[j] = c[j]+2*x*w[j+1]-w[j+2];
    }
    return w[0]-x*w[1];
}
double findminimum(double* V,int n){
    double min = V[0];
    for(int i =1;i<n;i++){
        if(min > V[i]){
            min = V[i];
        }
    }
    return min;
}
double findmaximum(double* V,int n){
    double max = V[0];
    for(int i =1;i<n;i++){
        if(max < V[i]){
            max = V[i];
        }
    }
    return max;
}
double* solve(int n,double** A,double* B,int* lambda){
    double* x = malloc(sizeof(double)*n);
    double sum;
    for(int k=0;k<n-1;k++){
        for(int i = k+1;i<n;i++){
            B[lambda[i]]=B[lambda[i]]-B[lambda[k]]*A[lambda[i]][k];
        }
    }
    x[n-1]=B[lambda[n-1]]/A[lambda[n-1]][n-1];
    for(int i = n-2;i>=0;i--){
        sum = B[lambda[i]];
        for(int j= i+1;j<n;j++){
            sum = sum -A[lambda[i]][j]*x[j];
        }
        x[i]= sum/A[lambda[i]][i];
    }
    return x;

}
double* Gauss(int n,double** A,double* B,int* lambda){
    double* S= malloc(sizeof(double)*4);
    double r,rmax,smax,xmult;
    int j;
    for(int i=0;i<n;i++){
        lambda[i] =i;
        smax =0;
        for(j=0;j<n;j++){
            if(A[i][j]>smax){
                smax = A[i][j];
            }
            if(-A[i][j]>smax){
                smax = -A[i][j];
            }
        }
        S[i]=smax;
    }
    for(int k=0;k<n-1;k++) {
        rmax = 0;
        for(int i=k;i<n;i++) {
            if(A[lambda[i]][k]/S[lambda[i]] >0){
                r = A[i][k]/S[i];
            }
            else{
                r = -A[i][k]/S[i];
            }
            if (r> rmax){
                rmax = r;
                j = i;
            }
        }

        int temp = lambda[j];
        lambda[j] = lambda[k];
        lambda[k] = temp;
        for(int i =k+1;i<n;i++){
            xmult = A[lambda[i]][k]/A[lambda[k]][k];
            A[lambda[i]][k]=xmult;
            for(j =k+1;j<n;j++){
                A[lambda[i]][j] = A[lambda[i]][j]-(xmult)*A[lambda[k]][j];
            }
        }
    }
    free(S);
    double* x = solve(n,A,B,lambda);
    return x;
}

double* fit(double* X, double* Y,int m){//up-and-in

    int* lambda= malloc(sizeof(int)*64);
    double a = findminimum(X,m);
    double b = findmaximum(X,m);
    for(int i = 1;i<m;i++){
        X[i] = (2*X[i] - a - b) / (b-a);
    }
    int n = 4;
    double** T;
    T = (double**)malloc(sizeof(double*)*(m));
    for (int i =0; i< m;i++){
        T[i] = (double*)malloc(sizeof(double)*n);
    }
    for(int k =0; k<m;k++){
        T[k][0]=1;
        T[k][1]=X[k];
        for(int j =2;j<n;j++){
            T[k][j]=2*X[k]*T[k][j-1]-T[k][j-2];
        }
    }

    double** A;
    A = (double**)malloc(sizeof(double*)*(n));
    for (int i =0; i< n;i++){
        A[i] = (double*)malloc(sizeof(double)*n);
    }
    double* B = malloc(sizeof(double)*n);
    double s;
    for(int i=0; i<n;i++){
        s =0;
        for(int k=0; k<m;k++){
            s = s +Y[k]*T[k][i];
        }
        B[i]=s;
        for(int j=i; j<n;j++){
            s = 0;
            for(int k=0; k<m;k++){
                s = s +T[k][j]*T[k][i];
            }
            A[i][j]=s;
            A[j][i]=s;
        }
    }
    double* x= Gauss(n,A,B,lambda);
    return x;
}
double normal2(){
    double y_1 = (double) rand()/ ((double) RAND_MAX );
    double y_2 = (double) rand()/ ((double) RAND_MAX );
    return sqrt(-2.0*log(y_1))*cos(M_PI*y_2);
}

int main(int argc, char **argv) {
    if(argc != 8){
        printf("Usage: ./call <K> <S_0> <T> <r> <q> <sigma> <H>");
        exit(1);
    }
    srand(time(0));
    clock_t start = clock();
    double K = atof(argv[1]);
    double S_0 = atof(argv[2]);
    double T = atof(argv[3]);
    double r = atof(argv[4]);
    double q = atof(argv[5]);
    double sigma = atof(argv[6]);
    double H = atof(argv[7]);
    int M = 100000; //I=M
    int N = 300;  //J=N
    double phi_alt[M];
    double** S;
    S = (double**)malloc(sizeof(double*)*M);
    for (int i =0; i< M;i++){
        S[i] = (double*)malloc(sizeof(double)*N);
    }
    bool** K_I = (bool**)malloc(sizeof(bool*)*M);
    for (int i =0; i< M;i++){
        K_I[i] = (bool*)malloc(sizeof(bool)*N);
    }
    double delta_t = T / (double) N;

    for(int i= 0; i < M; i++ ) {
        for (int j = 0; j < N; j++) {
            double Z= normal2();
            if(j == 0){
                S[i][j]= S_0;
            }
            else{
                S[i][j]=S[i][j-1] * exp((r - q - (sigma * sigma) / 2) * delta_t) * exp(sigma * sqrt(delta_t)  * Z);

            }

        }
    }
    for(int i= 0; i < M; i++ ) {
        bool knockedIn = false;
        for (int j = 0; j < N; j++) {
            if(knockedIn==true){
                K_I[i][j] = true;

            }else{
                if(S[i][j] < H){
                    K_I[i][j] =false;
                }else{
                    knockedIn= true;
                    K_I[i][j] = true;
                }
            }
        }
    }
    double** Phi =malloc(sizeof(double*)*M);
    for (int i =0; i< M;i++){
        Phi[i] = (double*)malloc(sizeof(double)*N);
    }
    for(int i=0; i< M-1; i++){
        if(K- S[i][N-1]>0 && K_I[i][N-1]==true){
            Phi[i][N-1]= K- S[i][N-1];
        }else{
            Phi[i][N-1]=0;
        }
    }
    for(int j= N-2; j > 0; j-- ){
        double a,b;
        double* x= malloc(sizeof(double)*4);
        double* X = malloc(sizeof(double)* M);
        double* Y = malloc(sizeof(double)* M);
        int counter = 0;
        for (int i = 0; i < M; i++) {
            if((K- S[i][j+1])>0&& K_I[i][j+1]==true) {
                X[counter] = S[i][j];

                Y[counter] = exp(-r * delta_t) * Phi[i][j+1];
                counter++;
            }
        }

        if(counter > 2) {
            a = findminimum(X,counter);
            b = findmaximum(X,counter);
            x  = fit(X, Y, counter);
        }
        for(int i = 0;i<M-1;i++){
            double phi = 0;
            if(K- S[i][j]>0&& K_I[i][j]==true){
                phi=K - S[i][j];
            }else{
                phi=0;
            }
            double Psi = 0;
            double g = chebyshev(4, (2 * S[i][j] - a - b) / (b - a), x);
            Psi = g;
            phi_alt[i] = 0;
            if(phi > Psi && phi > 0){
                Phi[i][j] = phi;
            }else{
                Phi[i][j] = exp(-r * delta_t) * Phi[i][j+1];
            }
            phi_alt[i] = Phi[i][j];
        }
        free(X);
        free(Y);
    }
    double sum =0;
    for (int i =0; i< M;i++){
        sum = sum +phi_alt[i];
        free(S[i]);
    }
    double average = sum / (double) M;
    clock_t end = clock();
    float seconds = (float)(end - start) / CLOCKS_PER_SEC;
//    printf("Price = %lf, calculated in %f seconds", average, seconds);
    printf("%f, %f, %f, %f, %f, %f, %f, %f, %f\n", K,S_0,T,r,q,sigma,H ,average, seconds);

    free(S);
    return 0;
}
