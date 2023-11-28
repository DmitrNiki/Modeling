#include "iostream"
#include <cmath>
#include "fstream"

const int step = 10;

//steps
const int T = 5000, K = 300, P = 11, P0 = 5;
const int SIZE = P * K;

//inittials
const int width = 20;
const int center = 50;
const int p0 = 4;

// p \in [-5; 5]
int index(int k, int p){
    return (p + P0) * K + k; 
}

double initial(int k){
    return std::exp(-1. * (k - center) * (k - center) / width / width);
}

double get_gamma(int p){
    return (p * 1.0) / P;
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
        for (int p = -P0; p <= P0; p++){
            total +=data[index(k, p)];
        }
        out << total << std::endl;
    }
    out.close();
}

void make_iteration(double* next, double* data){
    double gamma;
    for (int p = -P0; p <= P0; p++){
        gamma = get_gamma(p);
        if (gamma > 0){
            next[index(0, p)] = data[index(0, p)];
            for (int k = 1; k < K; k++){
                next[index(k, p)] = data[index(k, p)] -
		       	gamma * (data[index(k, p)] - data[(index(k - 1, p))]);    
            }
        } else {
            next[index(K - 1, p)] = data[index(K - 1, p)];
            for (int k = 0; k < K - 1; k++){
                next[index(k, p)] = data[index(k, p)] -
		       	gamma * (data[index(k + 1, p)] - data[(index(k, p))]);    
            }
        }
    }
// Зеркальное отражение
/*    for (int p = -P0; p <= P0; p++){
	   gamma = get_gamma(p);
	   if (gamma > 0){
		  next[index(0, p)] = next[index(0, -p)];
	   } else {
		  next[index(K - 1, p)] = next[index(K-1, -p)];
	   }
    }
*/

//Диффузное отражение
	double nom_left = 0, nom_right, denom = 0;
	for (int p = -P0; p <= P0; p++) {
		gamma = get_gamma(p);
		if(gamma > 0){
			nom_right += next[index(K - 1, p)] * gamma;
		} else {
			nom_left += next[index(0, p)] * gamma;
		}
		if (gamma > 0){
			denom += gamma * std::exp(-gamma * gamma / 2.);
		}
	}

	for (int p = -P0; p <= P0; p++){
		gamma = get_gamma(p);
		if (gamma > 0){
			next[index(0, p)] = -nom_left / denom * std::exp(-gamma * gamma / 2.);
		} else {
			next[index(K - 1, p)] = nom_right / denom * std::exp(-gamma * gamma / 2.);
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
