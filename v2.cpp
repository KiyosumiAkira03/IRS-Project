#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <random>
#include <algorithm>

using namespace std;
typedef complex<double> cd;

// ---------------------------------------------------------
// 1. Pure Linear Algebra Engine
// ---------------------------------------------------------
struct Matrix {
    int rows, cols;
    vector<vector<cd>> data;

    Matrix(int r, int c) : rows(r), cols(c), data(r, vector<cd>(c, 0.0)) {}

    // Matrix Multiplication O(n^3)
    Matrix operator*(const Matrix& other) const {
        if (cols != other.rows) throw invalid_argument("Matrix dimension mismatch for multiplication.");
        Matrix result(rows, other.cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < other.cols; ++j) {
                for (int k = 0; k < cols; ++k) {
                    result.data[i][j] += data[i][k] * other.data[k][j];
                }
            }
        }
        return result;
    }

    // Matrix Addition
    Matrix operator+(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) throw invalid_argument("Matrix dimension mismatch for addition.");
        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result.data[i][j] = data[i][j] + other.data[i][j];
            }
        }
        return result;
    }

    // Scalar Multiplication (Matrix Scaling)
    Matrix operator*(double scalar) const {
        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result.data[i][j] = data[i][j] * scalar;
            }
        }
        return result;
    }

    // Hermitian Transpose (Conjugate Transpose)
    Matrix hermitian() const {
        Matrix result(cols, rows);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result.data[j][i] = conj(data[i][j]);
            }
        }
        return result;
    }
    
    // Frobenius / L2 Norm (With explicit std:: scope)
    double norm() const {
        double sum = 0;
        for(int i = 0; i < rows; ++i) {
            for(int j = 0; j < cols; ++j) {
                sum += std::norm(data[i][j]); // Explicit std::norm resolves scope lookup error
            }
        }
        return sqrt(sum);
    }
};

// ---------------------------------------------------------
// 2. Hardware Phase Shift Constraint Model (Equation 5)
// ---------------------------------------------------------
const double BETA_MIN = 0.2;     // Paper parameters
const double PHI = 0.43 * M_PI;
const double K_FACTOR = 1.6;

double calculateAmplitude(double phase) {
    double sin_term = (sin(phase - PHI) + 1.0) / 2.0;
    return (1.0 - BETA_MIN) * pow(sin_term, K_FACTOR) + BETA_MIN;
}

// ---------------------------------------------------------
// 3. System Environment Constants
// ---------------------------------------------------------
const int N = 40; // Number of IRS elements
const int M = 2;  // Antennas at AP (Updated to match paper's simulation setup)

const double PT = pow(10, 36.0 / 10.0);       // 36 dBm Transmit Power converted to mW (approx 3981 mW)
const double NOISE_POWER = pow(10, -94.0 / 10.0); // -94 dBm Noise floor converted to mW

// ---------------------------------------------------------
// 4. Evolutionary Optimization Structures
// ---------------------------------------------------------
struct Individual {
    vector<double> phases; // Chromosome array: theta_n
    double fitness;        // Achievable rate (R_SE in bps/Hz)

    Individual() : phases(N), fitness(0.0) {}
};

// Evaluates fitness strictly using the exact paper equations
void evaluateFitness(Individual& ind, const Matrix& Phi, const Matrix& hd_H) {
    // A. Construct the 1 x N Reflective Row Matrix Vector: v^H
    Matrix v_H(1, N);
    for (int i = 0; i < N; ++i) {
        double theta = ind.phases[i];
        double beta = calculateAmplitude(theta);
        v_H.data[0][i] = conj(polar(beta, theta)); // Conjugate transpose state element
    }

    // B. Calculate Effective Channel Matrix: H_eff_H = v^H * Phi + hd_H
    Matrix H_eff_H = (v_H * Phi) + hd_H;

    // C. Derive Optimal Transmit Beamforming Vector via Matrix Notation
    Matrix H_eff = H_eff_H.hermitian(); 
    double norm_H = H_eff.norm();
    
    // w* = sqrt(P_T) * (H_eff / ||H_eff||)
    Matrix w = H_eff * (sqrt(PT) / norm_H);

    // D. Compute Achievable Rate Metric: R_SE = log2(1 + |H_eff_H * w|^2 / sigma^2)
    Matrix received_signal = H_eff_H * w;
    double signal_power = std::norm(received_signal.data[0][0]);
    double snr = signal_power / NOISE_POWER;
    
    ind.fitness = log2(1.0 + snr);
}

// ---------------------------------------------------------
// 5. Evolutionary Operators Engine
// ---------------------------------------------------------
double randomDouble(double min, double max) {
    static mt19937 generator(1337); 
    uniform_real_distribution<double> distribution(min, max);
    return distribution(generator);
}

Individual selectParent(const vector<Individual>& pop) {
    int idx1 = randomDouble(0, pop.size() - 1);
    int idx2 = randomDouble(0, pop.size() - 1);
    return (pop[idx1].fitness > pop[idx2].fitness) ? pop[idx1] : pop[idx2];
}

Individual crossover(const Individual& p1, const Individual& p2) {
    Individual child;
    for (int i = 0; i < N; ++i) {
        child.phases[i] = (randomDouble(0, 1) > 0.5) ? p1.phases[i] : p2.phases[i];
    }
    return child;
}

void mutate(Individual& ind, double mutation_rate) {
    for (int i = 0; i < N; ++i) {
        if (randomDouble(0, 1) < mutation_rate) {
            ind.phases[i] += randomDouble(-0.4, 0.4);
            if (ind.phases[i] > M_PI) ind.phases[i] -= 2 * M_PI;
            if (ind.phases[i] < -M_PI) ind.phases[i] += 2 * M_PI;
        }
    }
}

// ---------------------------------------------------------
// 6. Main Flow Setup
// ---------------------------------------------------------
int main() {
    cout << "Simulating Channel Environments per Section V..." << endl;

    // Raw channel vectors
    Matrix G(N, M);
    Matrix hr(N, 1);
    Matrix hd_H(1, M);

    // Initializing Channel Tensors with Rayleigh fading coefficients 
    for (int i = 0; i < N; ++i) {
        hr.data[i][0] = cd(randomDouble(-0.5, 0.5), randomDouble(-0.5, 0.5));
        for (int j = 0; j < M; ++j) G.data[i][j] = cd(randomDouble(-0.5, 0.5), randomDouble(-0.5, 0.5));
    }
    for (int j = 0; j < M; ++j) hd_H.data[0][j] = cd(randomDouble(-0.1, 0.1), randomDouble(-0.1, 0.1));

    // Calculate explicit intermediate matrix: Phi = diag(h_r^H) * G
    Matrix diag_hr_H(N, N);
    for (int i = 0; i < N; ++i) {
        diag_hr_H.data[i][i] = conj(hr.data[i][0]);
    }
    Matrix Phi = diag_hr_H * G;

    // Evolutionary Strategy parameters
    const int POP_SIZE = 60;
    const int GENERATIONS = 300;
    const double MUTATION_RATE = 0.15;

    vector<Individual> population(POP_SIZE);
    for (auto& ind : population) {
        for (int i = 0; i < N; ++i) ind.phases[i] = randomDouble(-M_PI, M_PI);
        evaluateFitness(ind, Phi, hd_H);
    }

    cout << "Executing Global Genetic Array Search Operations..." << endl;
    double best_global_fitness = 0;

    for (int gen = 1; gen <= GENERATIONS; ++gen) {
        vector<Individual> next_generation;

        auto best_it = max_element(population.begin(), population.end(), 
            [](const Individual& a, const Individual& b) { return a.fitness < b.fitness; });
        next_generation.push_back(*best_it);
        best_global_fitness = best_it->fitness;

        while (next_generation.size() < POP_SIZE) {
            Individual p1 = selectParent(population);
            Individual p2 = selectParent(population);
            Individual child = crossover(p1, p2);
            mutate(child, MUTATION_RATE);
            evaluateFitness(child, Phi, hd_H);
            next_generation.push_back(child);
        }
        population = next_generation;

        if (gen % 50 == 0) {
            cout << "Generation " << gen << " | Max Achievable Spectral Efficiency: " << best_global_fitness << " bps/Hz" << endl;
        }
    }

    cout << "\nOptimization Complete. Best Achievable Rate found: " << best_global_fitness << " bps/Hz" << endl;
    return 0;
}