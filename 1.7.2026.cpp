#include <bits/stdc++.h>

using namespace std;

// --- CẤU HÌNH BÀI TOÁN ---
const int N_ELEMENTS = 40;      // Số phần tử IRS (N)
const int M_ANTENNAS = 2;      // Số anten tại AP (M)
const int POPULATION_SIZE = 100;
const int GENERATIONS = 1000;
const double MUTATION_RATE = 0.15;
const double PI = 3.14159265358979323846;

// Thông số mô hình thực tế (từ paper)
const double BETA_MIN = 0.2;
const double K_CONST = 1.6;
const double PHI_CONST = 0.43 * PI;
const double PT = 3.98;         // Transmit power (36 dBm ~ 3.98W)
const double NOISE_POWER = 3.98e-10; // Noise power (-94 dBm)

// --- CẤU TRÚC SỐ PHỨC VÀ MA TRẬN (Tận dụng từ code gốc) ---
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
    vector<double> phases; // Bộ 40 tham số theta
    double fitness;        // Achievable Rate

    bool operator>(const Individual& other) const {
        return fitness > other.fitness;
    }
};

// --- DỮ LIỆU KÊNH TRUYỀN (Giả lập Rayleigh) ---
vector<comp> h_d; // AP-to-user (M x 1)
vector<vector<comp>> Phi; // Ma trận kết hợp diag(h_r^H)*G (N x M)

const double AP_USER = 3.8;
const double IRS_USER = 2.8;
const double AP_IRS = 2.2;
const double C0 = 1e-4;
const double d_AP_IRS = 500.0;

void init_channels(double d_AP_USER) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> rayleigh(0.5, 1.5);
    uniform_real_distribution<> phase_dist(-PI, PI);

    double path_loss_direct = sqrt(C0 * pow(d_AP_USER, -AP_USER));
    double path_loss_cascaded = sqrt(C0 * pow(d_AP_IRS, -AP_IRS)) * sqrt(C0 * pow(d_AP_IRS - d_AP_USER, -IRS_USER));
    
    h_d.resize(M_ANTENNAS);
    for (int m = 0; m < M_ANTENNAS; ++m) 
        h_d[m] = {path_loss_direct * rayleigh(gen), phase_dist(gen)};

    Phi.assign(N_ELEMENTS, vector<comp>(M_ANTENNAS));
    for (int n = 0; n < N_ELEMENTS; ++n) {
        for (int m = 0; m < M_ANTENNAS; ++m) {
            // Phi_nm = h_r_n* * G_nm
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
            // v_n^H (Conjugate transpose) = beta * e^(-j*theta)
            comp v_n_h = {b_n, -ind.phases[n]};
            sum_m = sumComplex(sum_m, multComplex(v_n_h, Phi[n][m]));
        }
        // Cộng với thành phần trực tiếp h_d^H
        comp h_d_m_h = {h_d[m].A, -h_d[m].theta};
        sum_m = sumComplex(sum_m, h_d_m_h);
        total_norm_sq += (sum_m.A * sum_m.A);
    }
    double snr = (total_norm_sq * PT) / NOISE_POWER;
    ind.fitness = log2(1.0 + snr);
}

// --- CÁC BƯỚC CỦA EVOLUTION STRATEGY ---
int main() { for (double d_AP_USER = 480.0; d_AP_USER <= 500.0; d_AP_USER += 2.5) {
    init_channels(d_AP_USER);
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> phase_dist(-PI, PI);
    uniform_real_distribution<> prob_dist(0.0, 1.0);

    // 1. Sinh ngẫu nhiên 100 bộ tham số ban đầu
    vector<Individual> population(POPULATION_SIZE);
    for (int i = 0; i < POPULATION_SIZE; ++i) {
        population[i].phases.resize(N_ELEMENTS);
        for (int n = 0; n < N_ELEMENTS; ++n) 
            population[i].phases[n] = phase_dist(gen);
        evaluate(population[i]);
    }

    Individual best_overall = population[0];

    // 2. Chạy qua 1000 thế hệ
    for (int g = 1; g <= GENERATIONS; ++g) {
        // Sắp xếp quần thể theo fitness giảm dần
        sort(population.begin(), population.end(), greater<Individual>());

        if (population[0].fitness > best_overall.fitness) {
            best_overall = population[0];
        }

        vector<Individual> next_gen;
        
        // Giữ lại Top 20 cá thể tốt nhất (Elitism)
        for (int i = 0; i < 20; ++i) {
            next_gen.push_back(population[i]);
        }

        // 3. Lai tạo và lấp đầy 80 chỗ trống để đủ 100
        while (next_gen.size() < POPULATION_SIZE) {
            // Chọn ngẫu nhiên 2 bố mẹ từ Top 50
            uniform_int_distribution<> parent_dist(0, 49);
            Individual p1 = population[parent_dist(gen)];
            Individual p2 = population[parent_dist(gen)];

            // Lai tạo (Crossover)
            Individual child;
            child.phases.resize(N_ELEMENTS);
            for (int n = 0; n < N_ELEMENTS; ++n) {
                child.phases[n] = (prob_dist(gen) < 0.5) ? p1.phases[n] : p2.phases[n];
                
                // Đột biến (Mutation) - thêm nhiễu Gaussian nhỏ
                if (prob_dist(gen) < MUTATION_RATE) {
                    child.phases[n] += phase_dist(gen) * 0.2; // Độ lệch đột biến
                    // Chuẩn hóa về [-PI, PI]
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
    cout << "Distance d_AP_User = " << d_AP_USER << " m | Achievable Rate: " << best_overall.fitness << " bps/Hz\n";
    }
    return 0;
}
