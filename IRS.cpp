#include <bits/stdc++.h>
using namespace std;
#define ll long long


// O(elements * generation) ~ O(n * m)
//possible optimization:
// +/ Heapify to find best individual (logn)
// +/ bin pow?
const int N_ELEMENTS = 40;           
const int POPULATION_SIZE = 100;    
const int GENERATIONS = 500;         
const double MUTATION_RATE = 0.1;    
const double PI = 3.14159265358979323846;
const double Gaussian_Noise = 1e-11;
// (from the paper)
const double BETA_MIN = 0.2;
const double K = 1.6;
const double PHI = 0.43 * PI;

random_device rd;
mt19937 gen(rd());
uniform_real_distribution<> phase_gacha(-PI, PI);
uniform_real_distribution<> percen_gacha(0.0, 1.0); 

struct Individual // pick random based on  "fitness" -> pick best -> discard the rest
{
    vector<double> phases; // (theta_1 to theta_N)
    double fitness;        // Achievable rate / signal strength
};


struct Matrix 
{
    double mat[N_ELEMENTS][N_ELEMENTS];

    Matrix()
    {
        
    }
};

Matrix multSca(Matrix A, long double k) //O(n^2)
{
    Matrix res;
    for (ll i = 0; i < 2; i++)
        for (ll j = 0; j < 2; j++)
            res.mat[i][j] = A.mat[i][j] * k;
    return res;
}

Matrix addMat(Matrix A, Matrix B) //O(n^2)
{
    Matrix res;
    for (ll i = 0; i < 2; i++)
        for (ll j = 0; j < 2; j++)
            res.mat[i][j] = A.mat[i][j] + B.mat[i][j];
    return res;
}

Matrix multMat(Matrix A, Matrix B) //O(n^3)
{
    Matrix res;
    for (ll i = 0; i < 2; ++i)
        for (ll j = 0; j < 2; ++j)
            for (ll k = 0; k < 2; ++k)
                res.mat[i][k] += A.mat[i][j] * B.mat[k][j];
    
    return res;
}



// amplitude (phase shift)
double amplitude(double theta) 
{
    // beta(theta) = (1 - beta_min) * ((sin(theta - phi) + 1) / 2)^k + beta_min
    double base = (sin(theta - PHI) + 1.0) / 2.0;
    return (1.0 - BETA_MIN) * pow(base, K) + BETA_MIN;
}

// double bin_pow(double a, int b) 
// {
//     if (!b) return 1;
//     long long x = bin_pow(a, b / 2);
//     if (b % 2 == 0) return x * x;
//     else return x * x * a;
// }

// P2 formula
void evaluate_fitness(Individual& index, const vector<double>& Psi, const vector<double>& phi_mag, const vector<double>& phi_arg) 
{
    double total_signal = 0.0;
    
    // max(theta) s.t (beta)^2 * psi 
    //                 + beta * phi_mag * cos(phi_arg - theta)
    for (int n = 0; n < N_ELEMENTS; ++n) 
    {
        double theta = index.phases[n];
        double beta = amplitude(theta);
        
        double term1 = beta * beta * Psi[n];
        double term2 = beta * phi_mag[n] * cos(phi_arg[n] - theta);
        
        total_signal += (term1 + term2);
    }
    index.fitness = log2(1.0 + total_signal);
}

// Create population
void begin_population(vector<Individual>& population) 
{
    for (int i = 0; i < POPULATION_SIZE; ++i) 
    {
        Individual index;
        //index.phases.reserve(N_ELEMENTS);
        for (int j = 0; j < N_ELEMENTS; ++j) index.phases.push_back(phase_gacha(gen));
        index.fitness = 0.0;
        population.push_back(index);
    }
}

//Battle royale
Individual select_parent(const vector<Individual>& population) 
{
    uniform_int_distribution<> index_dist(0, POPULATION_SIZE - 1);
    Individual best = population[index_dist(gen)];
    
    for (int i = 0; i < 3; ++i) // 3 chameleon royale
    { 
        Individual competitor = population[index_dist(gen)];
        if (competitor.fitness > best.fitness) best = competitor;
    }
    return best;
}

//Breed
Individual crossover(const Individual& parent1, const Individual& parent2) 
{
    Individual child;
    //child.phases.reserve(N_ELEMENTS);
    
    for (int i = 0; i < N_ELEMENTS; ++i) 
    {
        if (percen_gacha(gen) > 0.5) child.phases.push_back(parent1.phases[i]);
        else child.phases.push_back(parent2.phases[i]);
    }
    return child;
}

// Induce mutation (random noise)
void mutate(Individual& ind) 
{
    for (int i = 0; i < N_ELEMENTS; ++i) 
    {
        if (percen_gacha(gen) < MUTATION_RATE) 
        {
            ind.phases[i] += phase_gacha(gen) * 0.1; 
            
            // Keep within [-PI, PI]
            if (ind.phases[i] > PI) ind.phases[i] -= 2 * PI;
            if (ind.phases[i] < -PI) ind.phases[i] += 2 * PI;
        }
    }
}

int main() 
{
    // cai nay bo chiu, random
    vector<double> Psi(N_ELEMENTS, 0.05);
    vector<double> phi_mag(N_ELEMENTS, 0.8);
    vector<double> phi_arg(N_ELEMENTS);

    for(int i=0; i<N_ELEMENTS; ++i) phi_arg[i] = phase_gacha(gen);

    //create population
    vector<Individual> population;
    begin_population(population);
    for (auto& ind : population) evaluate_fitness(ind, Psi, phi_mag, phi_arg);
    

    //Start evolution
    Individual best_overall;
    best_overall.fitness = -1e9; 

    for (int generation = 0; generation < GENERATIONS; ++generation) 
    {
        vector<Individual> next_generation;
        
        // Battle Royale
        Individual best_in_gen;
        int current_fitness = -1e9;
        for (auto& ind : population)
        {
            if (ind.fitness > current_fitness)
            {
                current_fitness = ind.fitness;
                best_in_gen = ind;
            }
        }
        //auto best_in_gen = max_element(population.begin(), population.end(),
            //[](const Individual& a, const Individual& b) { return a.fitness < b.fitness; });
        next_generation.push_back(best_in_gen);
        
        if (best_in_gen.fitness > best_overall.fitness) best_overall = best_in_gen;


        // Breed new generation
        while (next_generation.size() < POPULATION_SIZE) 
        {
            Individual parent1 = select_parent(population);
            Individual parent2 = select_parent(population);
            
            Individual child = crossover(parent1, parent2);
            mutate(child);
            evaluate_fitness(child, Psi, phi_mag, phi_arg);
            
            next_generation.push_back(child);
        }
        
        population = next_generation;

        cout << "Generation " << generation << ": " << best_overall.fitness << '\n';
    }

    //for(int i=0; i<best_overall.phases.size(); ++i) cout << best_overall.phases[i] << " ";
    
}