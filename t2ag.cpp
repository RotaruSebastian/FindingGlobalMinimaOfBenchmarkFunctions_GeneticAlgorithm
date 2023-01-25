#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>

void set_bounds();
void init_var();
void reset_var();
double to_double(int individual, int nr);
double DeJong(int individual);
double Schwefel(int individual);
double Rastrigin(int individual);
double Michalewicz(int individual);
double func(int individual);
double fitness(int individual);
void elitism();
void clear_var();
void calc_var();
void select_trn();
void mutate();
void crossover();
void print_sol();
void Genetic();

std::vector <std::vector <bool>> population;
std::vector <std::vector <bool>> new_population;
std::vector <std::vector <bool>> population_copy;
std::vector <double> solution_copy;
std::vector <double> solution;
std::vector <double> I_fit;
std::vector <double> Q_prob;
int nr_var, nr_bits, nr_f, nr_elitism = 0, trn_K = 5, trn_J = 1;
double mutation_rate = 1.98, crossover_rate = 1;
bool b_elitism = true;

void set_bounds() {
    switch(nr_f) {
        case 0:
            nr_bits=20;
            break;
        case 1:
            nr_bits=27;
            break;
        case 2:
            nr_bits=20;
            break;
        case 3:
            nr_bits=19;
            break;
        default:
            std::cout<<"ERR1";
    }
}

void init_var() {
	std::vector <bool> line;
	for(int i = 0; i < 200; ++i) {
		for(int j = 0; j < nr_var * nr_bits; ++j) {
            line.push_back(rand() % 2);
        }
		population.push_back(line);
		line.clear();
	}
}

void reset_var() {
	population.clear();
	new_population.clear();
	population_copy.clear();
	solution_copy.clear();
	solution.clear();
	I_fit.clear();
	Q_prob.clear();
}

double to_double(int individual, int nr) {
	int n = nr_bits * (nr + 1);
	int var = 0, p = 1;
	for(int i = nr_bits * nr; i < n; ++i) {
		var += p * population[individual][i];
		p *= 2;
	}
	switch(nr_f) {
		case 0:
			return double(var) / (pow(2, nr_bits) - 1) * 10.24 - 5.12;
		case 1:
			return double(var) / (pow(2, nr_bits) - 1) * 1000 - 500;
		case 2:
			return double(var) / (pow(2, nr_bits) - 1) * 10.24 - 5.12;
		case 3:
			return double(var) / (pow(2, nr_bits) - 1) * M_PI;
		default:
			std::cout << "ERR5";
			return 0;
	}
}

double DeJong(int individual) {
	double rez = 0, xi;
	for(int i = 0; i < nr_var; ++i) {
		xi = to_double(individual, i);
		rez += pow(xi,2);
	}
	return rez;
}

double Schwefel(int individual) {
	double rez = 0, xi;
	for(int i = 0; i < nr_var; ++i) {
		xi = to_double(individual, i);
		rez += xi * sin(sqrt(std::abs(xi)));
	}
	return nr_var * 418.9829 - rez;
}

double Rastrigin(int individual) {
	double rez = 0, xi;
	for(int i = 0; i < nr_var; ++i) {
		xi = to_double(individual, i);
		rez += pow(xi,2) - 10 * cos(2 * xi * M_PI);
	}
	return 10 * nr_var + rez;
}

double Michalewicz(int individual) {
	double rez = 0, xi;
	for(int i = 0; i < nr_var; ++i) {
		xi = to_double(individual, i);
		rez -= sin(xi) * pow(sin(i * pow(xi,2) / M_PI),20);
	}
	return rez;
}

double func(int individual) {
	switch(nr_f) {
		case 0:
			return DeJong(individual);
		case 1:
			return Schwefel(individual);
		case 2:
			return Rastrigin(individual);
		case 3:
			return Michalewicz(individual);
		default:
			std::cout << "ERR2";
			return 0;
	}
}

double fitness(int individual) {
	switch(nr_f) {
		case 0:
			return 1 / pow((DeJong(individual) + 1), 2);
		case 1:
			return 1 / pow(nr_var * 418.9829 - Schwefel(individual), 2);
		case 2:
			return 1 / pow(10 * nr_var + Rastrigin(individual) + 1, 2);
		case 3:
			return 1 / (Michalewicz(individual) + 30);
		default:
			std::cout << "ERR3";
			return 0;
	}
}

void elitism() {
	nr_elitism = 25;
    population_copy = population;
    solution_copy = solution;
	bool ok = true;
	while(ok) {
		ok = false;
		for(int i = 0; i < 199; ++i) {
            if(solution_copy[i] < solution_copy[i + 1]) {
                std::swap(solution_copy[i], solution_copy[i + 1]);
                std::swap(population_copy[i], population_copy[i + 1]);
                ok = true;
            }
        }
	}
	for(int i = 0; i < nr_elitism; ++i) {
        new_population.push_back(population_copy[i]);
    }
}

void clear_var() {
	solution.clear();
	I_fit.clear();
	Q_prob.clear();
}

void calc_var() {
    clear_var();
	for(int i = 0; i < 200; ++i) {
        solution.push_back(fitness(i));
    }
}

void select_trn() {
	new_population.clear();
	if(b_elitism) {
        elitism();
    }
	for(int i = nr_elitism; i < 200; i += trn_J) {
		std::vector <int> r;
		std::vector <double> r_fit;
		for(int j = 0; j < trn_K; ++j) {
			r.push_back(rand() % 200);
			r_fit.push_back(fitness(r[j]));
			int k = j;
			while(k) {
				if(r_fit[k] > r_fit[k - 1]) {
					std::swap(r_fit[k], r_fit[k - 1]);
					std::swap(r[k],r[k - 1]);
				}
				--k;
			}
		}
		for(int k = 0; k < trn_J; ++k) {
            new_population.push_back(population[r[k]]);
        }
		r.clear();
		r_fit.clear();
	}
    population = new_population;
}

void mutate() {
	for(int i = nr_elitism; i < 200; ++i) {
        for(int j = 0; j < nr_var * nr_bits; ++j) {
            int r = rand() % 200;
            if(double(r) / 100 > mutation_rate) {
                population[i][j] = !population[i][j];
            }
        }
    }
}

void crossover() {
	int neighbours_interval = 200 - nr_elitism;
	for(int i = nr_elitism; i < 199; ++i) {
		int partner = rand() % neighbours_interval;
		for(int j = 0; j < nr_var * nr_bits; ++j) {
			int r = rand() % 201;
			if(double(r) / 100 > crossover_rate) {
                std::swap(population[i][j], population[partner + nr_elitism][j]);
            }
		}
	}
}

void print_sol() {
	double min_sol = 99999999;
	for(int i = 0; i < 200; ++i) {
		double individual_sol = func(i);
		if(individual_sol < min_sol) {
            min_sol = individual_sol;
        }
	}
	switch(nr_f) {
		case 0:
			std::cout << "DEJONG|";
			break;
		case 1:
			std::cout << "SCHWEFEL|";
			break;
		case 2:
			std::cout << "RASTRIGIN|";
			break;
		case 3:
			std::cout << "MICHALEWICZ|";
			break;
		default:
			std::cout << "ERR4|";
	}
	std::cout << "NR_VAR|" << nr_var << "|MIN_SOL|" << min_sol << ", |TIME|";
}

void Genetic() {
	auto start_time = std::chrono::high_resolution_clock::now();
	int t = 0;
    set_bounds();
    init_var();
    calc_var();
	while(t < 1000) {
        select_trn();
		mutate();
		crossover();
        clear_var();
        calc_var();
		++t;
	}
    select_trn();
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;
    print_sol();
	std::cout << time / std::chrono::milliseconds(1) << ", \n";
}

int main() {
	std::cout << std::fixed;
	std::cout.precision(5);
	for(nr_f = 0; nr_f < 4; ++nr_f) {
        nr_var = 5;
		srand(time(nullptr));
		Genetic();
        clear_var();
        reset_var();
        nr_var = 10;
		srand(time(nullptr));
		Genetic();
        clear_var();
        reset_var();
        nr_var=30;
		srand(time(nullptr));
		Genetic();
        clear_var();
        reset_var();
	}
	return 0;
}
