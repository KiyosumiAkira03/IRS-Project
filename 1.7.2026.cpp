#include <bits/stdc++.h>

using namespace std;

// --- CẤU HÌNH BÀI TOÁN ---
int N_ELEMENTS = 40;      // Số phần tử IRS (N)
const int M_ANTENNAS = 2;      // Số anten tại AP (M)
const int POPULATION_SIZE = 100;
const int GENERATIONS = 1000;
const double MUTATION_RATE = 0.15;
const double PI = 3.14159265358979323846;

// Thông số mô hình thực tế (từ paper)
const double BETA_MIN = 0.2;
const double K_CONST = 1.6;
const double PHI_CONST = 0.43 * PI;
const double PT = pow(10.0, (36 - 30.0) / 10.0);
const double NOISE_POWER = pow(10.0, (-94 - 30.0) / 10.0);

struct comp {
    double A;     // Biên độ (Amplitude)
    double theta; // Pha (Phase)
};

comp sumComplex(comp a, comp b) {
    double real = a.A * cos(a.theta) + b.A * cos(b.theta);
    double imag = a.A * sin(a.theta) + b.A * sin(b.theta);
    return {sqrt(real * real + imag * imag), atan2(imag, real)};
}

comp multComplex(comp a, comp b) {
    return {a.A * b.A, a.theta + b.theta};
}

// Hàm tính biên độ beta(theta) theo mô hình thực tế
double calculate_beta(double theta) {
    double base = (sin(theta - PHI_CONST) + 1.0) / 2.0;
    return (1.0 - BETA_MIN) * pow(base, K_CONST) + BETA_MIN;
}

// --- CẤU TRÚC THUẬT TOÁN TIẾN HÓA ---
struct Individual {
    vector<double> phases; // Bộ tham số theta
    double fitness;        // Achievable Rate

    bool operator>(const Individual& other) const {
        return fitness > other.fitness;
    }
};

// --- DỮ LIỆU KÊNH TRUYỀN (Giả lập Rayleigh) ---
vector<comp> h_d; // AP-to-user (M x 1)
vector<vector<comp>> Phi; // Ma trận kết hợp diag(h_r^H)*G (N x M)

// --- CẤU HÌNH KHOẢNG CÁCH VÀ SUY HAO ---
const double ALPHA_AP_USER = 3.8; 
const double ALPHA_IRS_USER = 2.8;
const double ALPHA_AP_IRS = 2.2;
const double C0 = 1e-4;
const double d_AP_IRS = 500.0;

random_device rd;
mt19937 gen(rd());
uniform_real_distribution<> rayleigh(0.5, 1.5);
uniform_real_distribution<> phase_dist(-PI, PI);
uniform_real_distribution<> prob_dist(0.0, 1.0);

void init_channels(double d_AP_USER) {
    // 1. Tính toán khoảng cách an toàn 
    double d_IRS_USER = abs(d_AP_IRS - d_AP_USER);
    if (d_IRS_USER < 1.0) d_IRS_USER = 1.0; 

    // 2. Tính toán Path Loss
    double path_loss_direct = sqrt(C0 * pow(d_AP_USER, -ALPHA_AP_USER));
    double path_loss_cascaded = sqrt(C0 * pow(d_AP_IRS, -ALPHA_AP_IRS)) * sqrt(C0 * pow(d_IRS_USER, -ALPHA_IRS_USER));
    
    h_d.resize(M_ANTENNAS);
    for (int m = 0; m < M_ANTENNAS; ++m) 
        h_d[m] = {path_loss_direct * rayleigh(gen), phase_dist(gen)};

    Phi.assign(N_ELEMENTS, vector<comp>(M_ANTENNAS));
    for (int n = 0; n < N_ELEMENTS; ++n) {
        for (int m = 0; m < M_ANTENNAS; ++m) {
            Phi[n][m] = {path_loss_cascaded * rayleigh(gen), phase_dist(gen)};
        }
    }
}

// Hàm đánh giá Fitness: Tính R_SE = log2(1 + SNR)
void evaluate(Individual& ind) {
    double total_norm_sq = 0;
    for (int m = 0; m < M_ANTENNAS; ++m) {
        comp sum_m = {0, 0};
        for (int n = 0; n < N_ELEMENTS; ++n) {
            double b_n = calculate_beta(ind.phases[n]);
            comp v_n_h = {b_n, -ind.phases[n]};
            sum_m = sumComplex(sum_m, multComplex(v_n_h, Phi[n][m]));
        }
        comp h_d_m_h = {h_d[m].A, -h_d[m].theta};
        sum_m = sumComplex(sum_m, h_d_m_h);
        total_norm_sq += (sum_m.A * sum_m.A);
    }
    double snr = (total_norm_sq * PT) / NOISE_POWER;
    ind.fitness = log2(1.0 + snr);
}

double run_evolution_strategy() {
    vector<Individual> population(POPULATION_SIZE);
    for (int i = 0; i < POPULATION_SIZE; ++i) {
        population[i].phases.resize(N_ELEMENTS);
        for (int n = 0; n < N_ELEMENTS; ++n) 
            population[i].phases[n] = phase_dist(gen);
        evaluate(population[i]);
    }

    Individual best_overall = population[0];

    for (int g = 1; g <= GENERATIONS; ++g) {
        sort(population.begin(), population.end(), greater<Individual>());

        if (population[0].fitness > best_overall.fitness) {
            best_overall = population[0];
        }

        vector<Individual> next_gen;
        
        // Elitism
        for (int i = 0; i < 20; ++i) {
            next_gen.push_back(population[i]);
        }

        while (next_gen.size() < POPULATION_SIZE) {
            uniform_int_distribution<> parent_dist(0, 49);
            Individual p1 = population[parent_dist(gen)];
            Individual p2 = population[parent_dist(gen)];

            Individual child;
            child.phases.resize(N_ELEMENTS);
            for (int n = 0; n < N_ELEMENTS; ++n) {
                child.phases[n] = (prob_dist(gen) < 0.5) ? p1.phases[n] : p2.phases[n];
                
                if (prob_dist(gen) < MUTATION_RATE) {
                    child.phases[n] += phase_dist(gen) * 0.2; 
                    if (child.phases[n] > PI) child.phases[n] -= 2*PI;
                    if (child.phases[n] < -PI) child.phases[n] += 2*PI;
                }
            }
            evaluate(child);
            next_gen.push_back(child);
        }
        population = next_gen;
         //if (g % 100 == 0) {

            //cout << "Generation " << g << " | Best Rate: " << best_overall.fitness << " bps/Hz" << endl;

        //}
    }
     // 4. Kết quả cuối cùng

    //cout << "\n--- OPTIMIZATION COMPLETED ---" << endl;

    //cout << "Final Best Rate: " << best_overall.fitness << " bps/Hz" << endl;

    //cout << "Best Theta sequence (all 40 elements): " << endl;

    //for (int i = 0; i < 40; ++i) cout << best_overall.phases[i] << " ";

    //cout << endl;
    return best_overall.fitness;
}

// --- CÁC BƯỚC CỦA EVOLUTION STRATEGY ---
int main() { 
    // Khảo sát khoảng cách
    cout << "--- SURVEY DISTANCE ---" << endl;
    for (double d_AP_USER = 480.0; d_AP_USER <= 500.0; d_AP_USER += 2.5) {
        init_channels(d_AP_USER);
        double best_rate = run_evolution_strategy();
        cout << "Distance d_AP_User = " << d_AP_USER << " m | Achievable Rate: " << best_rate << " bps/Hz\n";
    }

    // Khảo sát N (Số lượng phần tử IRS)
    cout << "\n--- SURVEY N_ELEMENTS ---" << endl;
    for (int current_N = 10; current_N <= 50; current_N += 5) {
        N_ELEMENTS = current_N;
        init_channels(498);
        double best_rate = run_evolution_strategy();
        cout << "N = " << N_ELEMENTS << " | Achievable Rate: " << best_rate << " bps/Hz\n";
    }
    
    return 0;
}
