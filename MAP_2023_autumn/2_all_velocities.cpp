#include "iostream"
#include <cmath>
#include "fstream"

const int step = 10;

//steps
const int T = 500, K = 100, P = 11, p0 = 5;
const int SIZE = P * K;

//inittials
const int a = 30;

// p \in [-5; 5]
int index(int k, int p){
    return (p + p0) * K + k; 
}

double initial(int k){
    return std::exp(-(k - a) * (k - a));
}

double get_gamma(int p){
    return (p * 2.0) / P;
}

void set_initial_values(double* data){
    for (int k = 0; k < K; k++){
        data[index(k, p0)] = initial(k);
    }
}

void dump_concentration(double* data, int i){
    std::ofstream out;
    char filename[100];
    sprintf(filename, "data/task2/out_%03d.dat", i);
    out.open(filename);

    for (int k = 0; k < K; k++){
        double total = 0;
        for (int p = -p0; p <= p0; p++){
            total +=data[index(k, p0)];
        }
        out << total << std::endl;
    }
    out.close();
}

void make_iteration(double* next, double* data){
    double gamma;
    for (int p = -p0; p <= p0; p++){
        gamma = get_gamma(p);
        if (gamma > 0){
            next[index(0, p)] = data[index(0, p)];
            for (int k = 1; k < K; k++){
                next[index(k, p)] = data[index(k, p)] - gamma * (data[index(k, p)] - data[(index(k - 1, p))]);    
            }
        } else {
            next[index(K - 1, p)] = data[index(K - 1, p)];
            for (int k = 0; k < K - 1; k++){
                next[index(k, p)] = data[index(k, p)] - gamma * (data[index(k + 1, p)] - data[(index(k, p))]);    
            }
        }
    }
}

int main(){
    double data[SIZE], next[SIZE];

    set_initial_values(data);

    for (int i = 0; i < T; i++){
        make_iteration(next, data);
        std::swap(data, next);

        if (i % step == 0){
            dump_concentration(data, i);
        }
    }
    
    return 0;
}
