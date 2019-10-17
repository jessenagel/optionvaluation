#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


//Function for calculating normal CDF as provided by John D. Cook, code in public domain
double phi(double x) {
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}

//Box-muller
double normal2(){
    double y_1 = (double) rand()/ ((double) RAND_MAX );
    double y_2 = (double) rand()/ ((double) RAND_MAX );
    return sqrt(-2.0*log(y_1))*cos(M_PI*y_2);
}
int main(int argc, char **argv) {
    if(argc != 8){
        printf("Usage: ./put <K> <S_0> <T> <r> <q> <sigma> <K>");
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


    int N = 100000;
    int M = 1000;
    double V=0;

    double delta_t = T / (double) M;
    double m = 2 * r / (sigma * sigma);
    double n = (2 * (r - q)) / (sigma * sigma);
    double delta_t_sqrt = sqrt(delta_t);
    for(int i=0; i < N ;i++){
        double S = S_0;
        double S_min1 = S;
        int check = 0;
        for (int j=0; j < M; j++){

            double t = ((double) j + 1)/ (double) M * T;
            double Theta = T - t;
            double Theta_sqrt = sqrt(Theta);
            if(S>= H){
                check =1;
            }
            if(check ==1) {
                if (t == T) {
                    if (K - S > 0) {
                        V = V + exp(-r * T) * (K - S);
                        break;
                    } else {
                        break;
                    }
                }
                double k = 1 - exp(-r * Theta);

                double Q_1 = (-(n - 1) - sqrt((n - 1) * (n - 1) + (4 * m) / k)) / 2;
                double d1 = (log((S / K)) + (Theta) * (r - q + (sigma * sigma) / 2)) / (sigma * Theta_sqrt);
                double d2 = d1 - sigma * Theta_sqrt;


                double e_qt = exp(-q * (Theta));
                double P = phi(-d2) * K * exp(-r * Theta) - S * e_qt * phi(-d1);

                double P_accent = e_qt * (phi(d1) - 1);
                double S_roof = (Q_1 * (K - P)) / (Q_1 - (1 + P_accent));
                if (S_roof > S) {
                    if (K - S > 0) {
                        V = V + exp(-r * t) * (K - S);
                        break;
                    } else {
                        break;
                    }
                }
            }
            double Z = normal2();
            S = S_min1 * exp((r - q - (sigma * sigma) / 2) * delta_t) * exp(sigma * delta_t_sqrt  *Z );
            S_min1 = S;
        }
    }


    clock_t end = clock();
    float seconds = (float)(end - start) / CLOCKS_PER_SEC;
    printf("%f, %f, %f, %f, %f, %f,%f, %f, %f\n", K,S_0,T,r,q,sigma,H, V/(double) N, seconds);

    return 0;
}
