#include <iostream> 
#include <cmath>
#include <fstream>

const int step = 10;

// Steps
const int T = 1000, K_x = 100, K_y = 100, P_x = 7, P_y = 7, 
      P0_x = P_x / 2, P0_y = P_y / 2;

const int SIZE = P_x * P_y * K_x * K_y;


//inittials
const int width = 20;
const int center_x = 50, center_y = 50;
const int p0_x = 2, p0_y = 3;
const int wall_x = 30;

// p \in [-5; 5]
int index(int k_x, int k_y, int p_x, int p_y){
    return (p_y + P0_y) * P_x * K_x * K_y + (p_x + P0_x) * K_x * K_y + k_y * K_x + k_x;
}

double get_gamma_x(int p){
    return (p * 1.0) / P_x;
}


double get_gamma_y(int p){
    return (p * 1.0) / P_y;
}

double initial_wall(int k_x, int k_y, int p_x, int p_y, int denom){
	double gamma_x = get_gamma_x(p_x);
	double gamma_y = get_gamma_y(p_y);
			
	if (k_x < wall_x){
		return std::exp(-1. * (gamma_x * gamma_x + gamma_y * gamma_y) / 2.) / denom;
	} else {
		return 1e-9;
	}
}

double initial_peak(int k_x, int k_y, int p_x, int p_y){
	if (p_x == p0_x && p_y == p0_y){
		return std::exp(-1. * ((k_x - center_x) * (k_x - center_x) + (k_y - center_y) * (k_y - center_y)) / width / width);
	} else {
		return 1e-9;
	}
}

void set_initial_values(double* data){
	double denom = 0, gamma_x, gamma_y;
	for(int p_x = -P0_x; p_x <=	 P0_x; p_x++){
		for (int p_y = -P0_y; p_y <= P0_y; p_y++){
			gamma_x = get_gamma_x(p_x);
			gamma_y = get_gamma_y(p_y);
			denom += std::exp(-1. * (gamma_x * gamma_x + gamma_y * gamma_y) / 2.);
		}
	}

    for (int k_x = 0; k_x < K_x; k_x++){
	    for (int k_y = 0; k_y < K_y; k_y++){
		    for(int p_x = -P0_x; p_x <= P0_x; p_x++){
			    for (int p_y = -P0_y; p_y <= P0_y; p_y++){
					data[index(k_x, k_y, p_x, p_y)] = initial_peak(k_x, k_y, p_x, p_y);
			    }
		    }
	    }
    }
}

void dump_partical_count(double* data){
    std::ofstream out;
    out.open("data/task3/total.dat", std::ios_base::app);
	
	double total = 0;
	for (int k_x = 0; k_x < K_x; k_x++){
	    for (int k_y = 0; k_y < K_y; k_y++){
		    for(int p_x = -P0_x; p_x <= P0_x; p_x++){
			    for (int p_y = -P0_y; p_y <= P0_y; p_y++){
					total += data[index(k_x, k_y, p_x, p_y)];
				}
			}
		}
	}
	out << total << std::endl;
	out.close();
}

void dump_concentration(double* data, int i){
    std::ofstream out;
    char filename[100];
    sprintf(filename, "data/task3/out_%03d.dat", i);
    out.open(filename);


    for (int k_x = 0; k_x < K_x; k_x++){
	    for (int k_y = 0; k_y < K_y; k_y++){
			double total = 0;
		    for(int p_x = -P0_x; p_x <= P0_x; p_x++){
			    for (int p_y = -P0_y; p_y <= P0_y; p_y++){
					total += data[index(k_x, k_y, p_x, p_y)];
				}
			}
			out << k_x << " " << k_y << " " <<  total << std::endl;
		}
	}
    out.close();
}


void make_iteration_x(double* next, double* data){
	for(int p_x = -P0_x; p_x <= P0_x; p_x++){
		for (int p_y = -P0_y; p_y <= P0_y; p_y++){
			double gamma = get_gamma_x(p_x);
			if (gamma > 0){
				for (int k_y = 0; k_y < K_y; k_y++){
					//next[index(0, k_y, p_x, p_y)] = data[index(0, k_y, p_x, p_y)];
					for (int k_x = 1; k_x < K_x; k_x++){
						next[index(k_x, k_y, p_x, p_y)] = data[index(k_x, k_y, p_x, p_y)] -
					gamma * (data[index(k_x, k_y, p_x, p_y)] - data[(index(k_x - 1, k_y, p_x, p_y))]);	
					}
				}
			} else {
				for (int k_y = 0; k_y < K_y; k_y++){
					//next[index(K_x - 1, k_y, p_x, p_y)] = data[index(K_x - 1, k_y, p_x, p_y)];
					for (int k_x = 0; k_x < K_x - 1; k_x++){
						next[index(k_x, k_y, p_x, p_y)] = data[index(k_x, k_y, p_x, p_y)] -
						gamma * (data[index(k_x + 1, k_y, p_x, p_y)] - data[(index(k_x, k_y, p_x, p_y))]);	
					}
				}
			}
		}
	}


	// //Диффузное отражение
	double nom_left[K_y] = {0}, nom_right[K_y] = {0}, denom = 0;
	for (int p_x = -P0_x; p_x <= P0_x; p_x++) {
		for (int p_y = -P0_y; p_y <= P0_y; p_y++) {
			double gamma_x = get_gamma_x(p_x);
			
			for (int k_y = 0; k_y < K_y; k_y++){
				if(gamma_x > 0){
					nom_right[k_y] += next[index(K_x - 1, k_y, p_x, p_y)] * gamma_x;
				} else {
					nom_left[k_y] += next[index(0, k_y, p_x, p_y)] * gamma_x;
				}
			}
			if (gamma_x > 0){
				denom += gamma_x * std::exp(-gamma_x * gamma_x / 2.);
			}
		}	
	}

	for (int p_x = -P0_x; p_x <= P0_x; p_x++) {
		for (int p_y = -P0_y; p_y <= P0_y; p_y++) {
			double gamma_x = get_gamma_x(p_x);
			
			for (int k_y = 0; k_y < K_y; k_y++){
				if(gamma_x > 0){		
					next[index(0, k_y, p_x, p_y)] = -nom_left[k_y] / denom * std::exp(-gamma_x * gamma_x / 2.);
				} else {
					next[index(K_x - 1, k_y, p_x, p_y)] = -nom_right[k_y] / (-denom) * std::exp(-gamma_x * gamma_x / 2.);
				}
			}
		}	
	}

	// Зеркальное отражение
	// for (int p_x = -P0_x; p_x <= P0_x; p_x++){
	// 	for(int p_y = -P0_y; p_y <= P0_y; p_y++){
	// 		double gamma = get_gamma_x(p_x);
	// 		for (int k_y = 0; k_y < K_y; k_y++){
	// 			if (gamma > 0){
	// 				next[index(0, k_y, p_x, p_y)] = next[index(0, k_y, -p_x, p_y)];
	// 			} else {
	// 				next[index(K_x - 1, k_y, p_x, p_y)] = next[index(K_x - 1, k_y, -p_x, p_y)];
	// 			}
	// 		}
	// 	}
	// }	

}

void make_iteration_y(double* next, double* data){
	for(int p_x = -P0_x; p_x <= P0_x; p_x++){
		for (int p_y = -P0_y; p_y <= P0_y; p_y++){
			double gamma = get_gamma_y(p_y);
			if (gamma > 0){
				for (int k_x = 0; k_x < K_x; k_x++){
					next[index(k_x, 0, p_x, p_y)] = data[index(k_x, 0, p_x, p_y)];
					for (int k_y = 1; k_y < K_y; k_y++){
						next[index(k_x, k_y, p_x, p_y)] = data[index(k_x, k_y, p_x, p_y)] -
					gamma * (data[index(k_x, k_y, p_x, p_y)] - data[(index(k_x, k_y - 1, p_x, p_y))]);	
					}
				}
			} else {
				for (int k_x = 0; k_x < K_x; k_x++){
					next[index(k_x, K_y - 1, p_x, p_y)] = data[index(k_x, K_y - 1, p_x, p_y)];
					for (int k_y = 0; k_y < K_y - 1; k_y++){
						next[index(k_x, k_y, p_x, p_y)] = data[index(k_x, k_y, p_x, p_y)] -
						gamma * (data[index(k_x, k_y + 1, p_x, p_y)] - data[(index(k_x, k_y, p_x, p_y))]);	
					}
				}
			}
		}
	}

double nom_left[K_y] = {0}, nom_right[K_y] = {0}, denom = 0;
	for (int p_x = -P0_x; p_x <= P0_x; p_x++) {
		for (int p_y = -P0_y; p_y <= P0_y; p_y++) {
			double gamma = get_gamma_y(p_y);
			
			for (int k_x = 0; k_x < K_x; k_x++){
				if(gamma > 0){
					nom_right[k_x] += next[index(k_x, K_y - 1, p_x, p_y)] * gamma;
				} else {
					nom_left[k_x] += next[index(k_x, 0, p_x, p_y)] * gamma;
				}
			}
			if (gamma > 0){
				denom += gamma * std::exp(-gamma * gamma / 2.);
			}
		}	
	}

	for (int p_x = -P0_x; p_x <= P0_x; p_x++) {
		for (int p_y = -P0_y; p_y <= P0_y; p_y++) {
			double gamma = get_gamma_y(p_y);
			
			for (int k_x = 0; k_x < K_x; k_x++){
				if(gamma > 0){		
					next[index(k_x, 0, p_x, p_y)] = -nom_left[k_x] / denom * std::exp(-gamma * gamma / 2.);
				} else {
					next[index(k_x, K_y - 1, p_x, p_y)] = -nom_right[k_x] / (-denom) * std::exp(-gamma * gamma / 2.);
				}
			}
		}	
	}

	//Зеркальное отражение
	// for (int p_x = -P0_x; p_x <= P0_x; p_x++){
	// 	for(int p_y = -P0_y; p_y <= P0_y; p_y++){
	// 		double gamma = get_gamma_y(p_y);
	// 		for (int k_x = 0; k_x < K_x; k_x++){
	// 			if (gamma > 0){
	// 				next[index(k_x, 0, p_x, p_y)] = next[index(k_x, 0, p_x, -p_y)];
	// 			} else {	
	// 				next[index(k_x, K_y - 1, p_x, p_y)] = next[index(k_x, K_y - 1, p_x, -p_y)];
	// 			}
	// 		}
	// 	}
	// }
}

int main(){
    double data[SIZE], next[SIZE];

    set_initial_values(data);

    for (int i = 0; i < T; i++){
        make_iteration_x(next, data);
        std::swap(data, next);

        make_iteration_y(next, data);
        std::swap(data, next);

        if (i % step == 0){
            dump_concentration(data, i);
			dump_partical_count(data);
        }
    }
    
    return 0;
}

