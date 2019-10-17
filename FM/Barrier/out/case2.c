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

double norm(double x){
    return  1/sqrt(2*M_PI)*exp(-0.5*x*x);

}

int main(int argc, char **argv) {
    if(argc != 8){
        printf("Usage: ./put <K> <S_0> <T> <r> <q> <sigma> <H>\n");
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
    int M = 200;
    double V[N];
    double delta_t = T / (double) M;
    double m = 2 * r / (sigma * sigma);
    double n = (2 * (r - q)) / (sigma * sigma);
//    double exponentfirstpart = exp((r - q - (sigma * sigma) / 2) * delta_t);
    double delta_t_sqrt = sqrt(delta_t);
    double lambda = (r-q)/(sigma*sigma) +0.5;

    for (int i = 0; i < N; i++) {
        double S = S_0;
        double S_min1 = S_0;
        for (int j = 0; j < M; j++) {
            double t = ((double) j + 1) / (double) M * T;
            double Theta = T - t;
            double Theta_sqrt = sqrt(Theta);

            if (S > H) {
                V[i] = exp(-r*t)*(K-H);
                break;
            }

            if (t == T) {
                if (K - S > 0) {
                    V[i] = exp(-r * T) * (K - S);
                    break;
                } else {
                    V[i] = 0;
                    break;
                }
            }


            // Critical price
            double R = K-H;
            double d1 = (log((S / K)) + (Theta) * (r - q + (sigma * sigma) / 2)) / (sigma * Theta_sqrt);
            double d1_H_2 = (log((H * H / S) / K) + (Theta) * (r - q + (sigma * sigma) / 2)) / (sigma * Theta_sqrt);
            double d2 = d1 - sigma * Theta_sqrt;
            double d2_H_2 = d1_H_2 - sigma * Theta_sqrt;

            double e_qt = exp(-q * (Theta));
            double p_0_H_2 = K * exp(-r * Theta) * phi(-d2_H_2)  - (H * H / S) * e_qt * phi(-d1_H_2);
            double p_0 = K * exp(-r * Theta) * phi(-d2) - S * e_qt * phi(-d1);
            double p_1 = p_0 - pow((H / S), (2 * lambda - 2)) * p_0_H_2;

            double u = (S + 5 * H) / H;
            double m_accent = m * (1 + u * (1 / (1 - exp(-u * r * Theta)) - 1));
            double x_1= (log(S/H))/(sigma*sqrt(Theta))+lambda*sigma*sqrt(Theta);
            double x_2= x_1 - sigma * sqrt(Theta);
            double y_1= (log(H/S))/(sigma*sqrt(Theta))+lambda*sigma*sqrt(Theta);
            double y_2= y_1 - sigma * sqrt(Theta);
            double v = sqrt((lambda -1)*(lambda-1)+(2*r)/(sigma*sigma));
            double z_1= (log(H/S))/(sigma*sqrt(Theta))+v*sigma*sqrt(Theta);
            double z_2 = z_1 - 2*v*sigma*sqrt(Theta);
            double beta_minus = (-(n - 1) - sqrt((n - 1) * (n - 1) + 4 * m_accent)) / 2;
            double beta_plus = (-(n - 1) + sqrt((n - 1) * (n - 1) + 4 * m_accent)) / 2;
            double delta_1 = -((exp(-q*Theta)*phi(x_1))/(sigma*sqrt(Theta)))*(K/H-1)-exp(-q*Theta)*phi(-x_1);
            double delta_2 = pow((H/S),2*lambda)*((K*S/(H*H))*exp(-r*Theta)*(2*lambda-2)*phi(-y_2)+exp(-q*Theta)*(1-2*lambda)*phi(-y_1)-(exp(-q*Theta)*phi(y_1))/(sigma*sqrt(Theta))*(K/H-1));
            double delta_3 = (R/S)*(pow(H/S,lambda-1+v)*(phi(z_1)/(sigma*sqrt(Theta))-(lambda-1+v)*(phi(-z_1)))+(pow(H/S,lambda-1-v)*(phi(z_2)/(sigma*sqrt(Theta))-(lambda-1-v)*(phi(-z_2)))));
            double W_2 = ((-1 - ( delta_1 +delta_2 + delta_3)) );
            double p_2_1= K*exp(-r*Theta)*phi(-x_2)-S*exp(-q*Theta)*phi(-x_1);
            double p_2_2= -K*exp(-r*Theta)*pow(H/S,2*lambda-2)*phi(-y_2)+S*exp(-q*Theta)*pow(H/S,2*lambda)*phi(-y_1);
            double p_2_3= R*(pow(H/S,lambda-1+v)*phi(-z_1)+pow(H/S,lambda-1-v)*phi(-z_2));
            double p_2 = p_2_1 +p_2_2+ p_2_3;
//            printf("p1=%f,p2=%f,p3=%f\n",p_2_1,p_2_2,p_2_3);
            double S_roof_2 = K - p_2 - W_2 * (pow(S, beta_minus) - pow(H, (beta_minus - beta_plus)) * pow(S, beta_plus))/
                                              (beta_minus * pow(S, (beta_minus - 1)) - beta_plus * pow(H, (beta_minus - beta_plus)) * pow(S, (beta_plus - 1)));
//            printf("beta-=%f\n",(pow(S, beta_minus) - pow(H, (beta_minus - beta_plus)) * pow(S, beta_plus))/ (beta_minus * pow(S, (beta_minus - 1)) - beta_plus * pow(H, (beta_minus - beta_plus)) * pow(S, (beta_plus - 1))));
//            printf("deel1=%f,deel2=%f",x_1,x_2);
//            printf("S^=%f,K=%f,p_2=%f,w_2-=%f\n",S_roof_2,K,p_2,W_2);
            if (S < S_roof_2) {
                V[i] = exp(-r * t) * (K - S);
                break;

            }
            double Z = normal2();
            S = S_min1 * exp((r - q - (sigma * sigma) / 2) * delta_t + sigma * delta_t_sqrt * Z);
            S_min1 = S;
        }
    }

    double sum = 0;
    for (int a = 0; a < N; ++a) {
        sum = sum + V[a];
//        printf("%f\n",V[a]);
    }
    double average = sum / N;
    clock_t end = clock();
    float seconds = (float) (end - start) / CLOCKS_PER_SEC;
    printf("%f, %f, %f, %f, %f, %f, %f, %f, %f\n", K,S_0,T,r,q,sigma,H ,average, seconds);

    return 0;
}