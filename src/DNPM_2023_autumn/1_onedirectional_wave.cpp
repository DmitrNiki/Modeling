#include "iostream"
#include <cmath>
#include "fstream"

const int step = 10;
const int T = 100, N = 100;
const int a = 3;

// const double h = 1., tau = 1.;
// const double v = 4.;


double initial(int k){
    return std::exp(-(k - a) * (k - a));
}

int main(){
    const double gamma = 0.5;
    double data[N], next[N];

    std::ofstream out;
    char filename[100];
    for (int k =0; k < N; k++){
        data[k] = initial(k);
    }

    for (int i = 0; i < T; i++){
        next[0] = data[0];
        for(int k = 1; k < N; k++){
            next[k] = data[k] - gamma * (data[k] - data[k - 1]);
        }
        std::swap(data, next);
        if (i % step == 0){
            sprintf(filename, "data/task1/out_%03d.dat", i);
            out.open(filename);
            for (int k = 0; k < N; k++){
                out << data[k] << std::endl;
            }
            out .close();
        }
    }
    
    return 0;
}
