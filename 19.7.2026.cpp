#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <random>
#include <algorithm>
#include <string>
#include <iomanip>

using namespace std;

// --- CẤU HÌNH BÀI TOÁN ---
int N_ELEMENTS = 40;      // Số phần tử IRS (N)
const int M_ANTENNAS = 2; // Số anten tại AP (M)
const int POPULATION_SIZE = 100;
const int GENERATIONS = 1000;
const double MUTATION_RATE = 0.15;
const double PI = 3.14159265358979323846;

const double PT = pow(10.0, (36 - 30.0) / 10.0);
const double NOISE_POWER = pow(10.0, (-94 - 30.0) / 10.0);

// --- THÔNG SỐ VẬT LÝ THEO BÀI BÁO ---
const double W = 2.0 * PI * 2.4e9; // Tần số góc omega = 2*PI*2.4GHz
const double Z0 = 377.0;           // Trở kháng không gian tự do (Ohm)

// Ràng buộc vật lý của linh kiện
const double C_MIN = 0.466e-12; // 0.466 pF
const double C_MAX = 2.35e-12;  // 2.35 pF
const double R_MIN = 0.47;      // 0.47 Ohm 
const double R_MAX = 30.8;      // 30.8 Ohm
const double L1_MIN = 0.5e-9;   // 0.5 nH 
const double L1_MAX = 13.0e-9;  // 13.0 nH
const double L2_MIN = 0.5e-9;   // 0.5 nH 
const double L2_MAX = 13.0e-9;  // 13.0 nH

// ==============================================================
// --- CẤU HÌNH CỐ ĐỊNH THÔNG SỐ ---
// ==============================================================
const bool FIX_C = false;       // Đổi thành true nếu muốn cố định
const double FIXED_C_VAL = 1.408e-12; // 1.408 pF

const bool FIX_R = false;
const double FIXED_R_VAL = 2.0;       // 2 Ohm

const bool FIX_L1 = true;
const double FIXED_L1_VAL = 2.5e-9;   // 2.5 nH

const bool FIX_L2 = true;
const double FIXED_L2_VAL = 0.7e-9;   // 0.7 nH
// ==============================================================

// Thông số của Practical Phase Shift Model (Bài báo)
const double BETA_MIN = 0.2;
const double MODEL_K = 1.6;
const double MODEL_PHI = 0.43 * PI;

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

// --- CẤU TRÚC PHẦN TỬ IRS ---
struct IRSElement {
    double C;
    double R;
    double L1;
    double L2;
};

// --- CẤU TRÚC THUẬT TOÁN TIẾN HÓA ---
struct Individual {
    vector<IRSElement> elements; // Bộ tham số vật lý của N phần tử
    double fitness;              // Achievable Rate (R_SE)

    bool operator>(const Individual& other) const {
        return fitness > other.fitness;
    }
};

// Hàm tính v_n trực tiếp từ phương trình mạch cộng hưởng song song
comp calculate_vn_from_physics(const IRSElement& el) {
    complex<double> j(0.0, 1.0);
    complex<double> Z_L1 = j * W * el.L1;
    complex<double> Z_L2 = j * W * el.L2;
    complex<double> Z_C = 1.0 / (j * W * el.C);
    
    // Tính Z_n theo công thức trong bài báo
    complex<double> branch2 = Z_L2 + Z_C + el.R;
    complex<double> Zn = (Z_L1 * branch2) / (Z_L1 + branch2);
    
    // Tính v_n
    complex<double> vn = (Zn - Z0) / (Zn + Z0);
    
    return {abs(vn), arg(vn)};
}

// --- DỮ LIỆU KÊNH TRUYỀN ---
vector<comp> h_d; 
vector<vector<comp>> Phi; 

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
normal_distribution<> norm_dist(0.0, 1.0);

void init_channels(double d_AP_USER) {
    double d_IRS_USER = abs(d_AP_IRS - d_AP_USER);
    if (d_IRS_USER < 1.0) d_IRS_USER = 1.0; 

    double d_AP_USER_true = sqrt(d_AP_USER * d_AP_USER + 2.0 * 2.0);
    double d_IRS_USER_true = sqrt((d_AP_IRS - d_AP_USER) * (d_AP_IRS - d_AP_USER) + 2.0 * 2.0);

    double path_loss_direct = sqrt(C0 * pow(d_AP_USER_true, -ALPHA_AP_USER));
    double path_loss_cascaded = sqrt(C0 * pow(d_AP_IRS, -ALPHA_AP_IRS)) * sqrt(C0 * pow(d_IRS_USER_true, -ALPHA_IRS_USER));
    
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

// Đánh giá Fitness cho ES
void evaluate(Individual& ind) {
    double total_norm_sq = 0;
    for (int m = 0; m < M_ANTENNAS; ++m) {
        comp sum_m = {0, 0};
        for (int n = 0; n < N_ELEMENTS; ++n) {
            comp v_n = calculate_vn_from_physics(ind.elements[n]);
            comp v_n_h = {v_n.A, -v_n.theta}; 
            sum_m = sumComplex(sum_m, multComplex(v_n_h, Phi[n][m]));
        }
        comp h_d_m_h = {h_d[m].A, -h_d[m].theta};
        sum_m = sumComplex(sum_m, h_d_m_h);
        total_norm_sq += (sum_m.A * sum_m.A);
    }
    double snr = (total_norm_sq * PT) / NOISE_POWER;
    ind.fitness = log2(1.0 + snr);
}

// Khởi tạo phần tử, cân nhắc xem có đang cố định thông số nào không
IRSElement random_element() {
    uniform_real_distribution<> C_dist(C_MIN, C_MAX);
    uniform_real_distribution<> R_dist(R_MIN, R_MAX);
    uniform_real_distribution<> L1_dist(L1_MIN, L1_MAX);
    uniform_real_distribution<> L2_dist(L2_MIN, L2_MAX);
    
    return {
        FIX_C  ? FIXED_C_VAL  : C_dist(gen),
        FIX_R  ? FIXED_R_VAL  : R_dist(gen),
        FIX_L1 ? FIXED_L1_VAL : L1_dist(gen),
        FIX_L2 ? FIXED_L2_VAL : L2_dist(gen)
    };
}

double clamp(double val, double min_val, double max_val) {
    if (val < min_val) return min_val;
    if (val > max_val) return max_val;
    return val;
}

// --- THUẬT TOÁN 1: EVOLUTION STRATEGY (ES) ---
double run_evolution_strategy() {
    vector<Individual> population(POPULATION_SIZE);
    for (int i = 0; i < POPULATION_SIZE; ++i) {
        population[i].elements.resize(N_ELEMENTS);
        for (int n = 0; n < N_ELEMENTS; ++n) 
            population[i].elements[n] = random_element();
        evaluate(population[i]);
    }

    Individual best_overall = population[0];
    double sigma = 0.2; 

    for (int g = 1; g <= GENERATIONS; ++g) {
        sort(population.begin(), population.end(), greater<Individual>());
        if (population[0].fitness > best_overall.fitness) best_overall = population[0];

        vector<Individual> next_gen;
        for (int i = 0; i < 20; ++i) next_gen.push_back(population[i]);

        int success_count = 0;
        int total_children = 0;

        while (next_gen.size() < POPULATION_SIZE) {
            uniform_int_distribution<> parent_dist(0, 49);
            Individual p1 = population[parent_dist(gen)];
            Individual p2 = population[parent_dist(gen)];

            Individual child;
            child.elements.resize(N_ELEMENTS);
            for (int n = 0; n < N_ELEMENTS; ++n) {
                // Lai ghép: Nếu bị cố định thì giữ nguyên giá trị cố định
                child.elements[n].C  = FIX_C  ? FIXED_C_VAL  : ((prob_dist(gen) < 0.5) ? p1.elements[n].C  : p2.elements[n].C);
                child.elements[n].R  = FIX_R  ? FIXED_R_VAL  : ((prob_dist(gen) < 0.5) ? p1.elements[n].R  : p2.elements[n].R);
                child.elements[n].L1 = FIX_L1 ? FIXED_L1_VAL : ((prob_dist(gen) < 0.5) ? p1.elements[n].L1 : p2.elements[n].L1);
                child.elements[n].L2 = FIX_L2 ? FIXED_L2_VAL : ((prob_dist(gen) < 0.5) ? p1.elements[n].L2 : p2.elements[n].L2);
                
                // Đột biến: Chỉ đột biến các thông số KHÔNG bị cố định
                if (prob_dist(gen) < MUTATION_RATE) {
                    if (!FIX_C)  child.elements[n].C  = clamp(child.elements[n].C  + norm_dist(gen) * sigma * (C_MAX - C_MIN), C_MIN, C_MAX);
                    if (!FIX_R)  child.elements[n].R  = clamp(child.elements[n].R  + norm_dist(gen) * sigma * (R_MAX - R_MIN), R_MIN, R_MAX);
                    if (!FIX_L1) child.elements[n].L1 = clamp(child.elements[n].L1 + norm_dist(gen) * sigma * (L1_MAX - L1_MIN), L1_MIN, L1_MAX);
                    if (!FIX_L2) child.elements[n].L2 = clamp(child.elements[n].L2 + norm_dist(gen) * sigma * (L2_MAX - L2_MIN), L2_MIN, L2_MAX);
                }
            }
            evaluate(child);
            total_children++;
            if (child.fitness > max(p1.fitness, p2.fitness)) success_count++;
            next_gen.push_back(child);
        }
        population = next_gen;

        if (total_children > 0) {
            double success_rate = (double)success_count / total_children;
            if (success_rate > 0.2) sigma /= 0.85; 
            else if (success_rate < 0.2) sigma *= 0.85; 
            sigma = max(0.001, min(sigma, 1.0)); 
        }
    }
    return best_overall.fitness;
}

// --- HÀM PHỤ TRỢ CHO 5 BASELINES BÀI BÁO ---
double beta_model(double theta) {
    return (1.0 - BETA_MIN) * pow((sin(theta - MODEL_PHI) + 1.0) / 2.0, MODEL_K) + BETA_MIN;
}

// Tính Rate dựa trên vector v (Dùng chung cho 5 baselines)
double calculate_baseline_rate(const vector<complex<double>>& v, const vector<vector<complex<double>>>& Phi_c, const vector<complex<double>>& hd_c) {
    double total_norm_sq = 0;
    for (int m = 0; m < M_ANTENNAS; ++m) {
        complex<double> sum_m = 0;
        for (int n = 0; n < N_ELEMENTS; ++n) {
            sum_m += conj(v[n]) * Phi_c[n][m];
        }
        sum_m += conj(hd_c[m]);
        total_norm_sq += norm(sum_m);
    }
    double snr = (total_norm_sq * PT) / NOISE_POWER;
    return log2(1.0 + snr);
}

// --- THUẬT TOÁN 2-6: 5 ĐƯỜNG BASELINE CỦA BÀI BÁO ---
double run_baseline(int method_id) {
    vector<vector<complex<double>>> Phi_c(N_ELEMENTS, vector<complex<double>>(M_ANTENNAS));
    vector<complex<double>> hd_c(M_ANTENNAS);

    // Convert custom struct to std::complex for matrix operations
    for (int n = 0; n < N_ELEMENTS; n++) {
        for (int m = 0; m < M_ANTENNAS; m++) Phi_c[n][m] = polar(Phi[n][m].A, Phi[n][m].theta);
    }
    for (int m = 0; m < M_ANTENNAS; m++) hd_c[m] = polar(h_d[m].A, h_d[m].theta);

    if (method_id == 5) {
        // Lower Bound: No IRS (v = 0)
        vector<complex<double>> v(N_ELEMENTS, 0.0);
        return calculate_baseline_rate(v, Phi_c, hd_c);
    }

    // Tính ma trận Psi (N x N) và hat_hd (N x 1)
    vector<vector<complex<double>>> Psi(N_ELEMENTS, vector<complex<double>>(N_ELEMENTS, 0.0));
    for (int i = 0; i < N_ELEMENTS; ++i) {
        for (int j = 0; j < N_ELEMENTS; ++j) {
            for (int m = 0; m < M_ANTENNAS; ++m) Psi[i][j] += Phi_c[i][m] * conj(Phi_c[j][m]);
        }
    }

    vector<complex<double>> hat_hd(N_ELEMENTS, 0.0);
    for (int i = 0; i < N_ELEMENTS; ++i) {
        for (int m = 0; m < M_ANTENNAS; ++m) hat_hd[i] += Phi_c[i][m] * hd_c[m];
    }

    // Khởi tạo vector pha theta
    vector<double> theta(N_ELEMENTS);
    vector<complex<double>> v(N_ELEMENTS);
    for (int n = 0; n < N_ELEMENTS; n++) {
        theta[n] = (prob_dist(gen) > 0.5) ? PI : -PI;
        v[n] = polar(beta_model(theta[n]), theta[n]);
    }

    // Vòng lặp Alternating Optimization (AO)
    for (int iter = 0; iter < 30; iter++) {
        for (int n = 0; n < N_ELEMENTS; n++) {
            complex<double> Y = 0;
            for (int m = 0; m < N_ELEMENTS; m++) {
                if (m != n) Y += Psi[n][m] * v[m];
            }
            Y += hat_hd[n];
            complex<double> varphi_n = 2.0 * Y;

            double arg_phi = arg(varphi_n);
            double abs_phi = abs(varphi_n);
            double psi_nn = real(Psi[n][n]);

            if (method_id == 1 || method_id == 4) { 
                // Ideal IRS hoặc Ideal IRS Assumption
                theta[n] = arg_phi;
                v[n] = polar(1.0, theta[n]); // beta = 1
            } 
            else if (method_id == 2) { 
                // Practical IRS (AO with Proposition 1)
                double lambda = (arg_phi >= 0) ? 0 : 1;
                double theta_A = arg_phi;
                double theta_C = (lambda == 0) ? PI : -PI;
                double theta_B = (theta_A + theta_C) / 2.0;

                auto obj_func = [&](double th) {
                    double b = beta_model(th);
                    return b * b * psi_nn + b * abs_phi * cos(arg_phi - th);
                };

                double f1 = obj_func(theta_A), f2 = obj_func(theta_B), f3 = obj_func(theta_C);
                double num = theta_C * (3*f1 - 4*f2 + f3) + theta_A * (f1 - 4*f2 + 3*f3);
                double den = 4.0 * (f1 - 2*f2 + f3);
                
                double theta_opt = (abs(den) < 1e-9) ? ((f1 > f3) ? theta_A : theta_C) : (num / den);
                
                double min_th = min(theta_A, theta_C);
                double max_th = max(theta_A, theta_C);
                theta_opt = clamp(theta_opt, min_th, max_th);

                double f_opt = obj_func(theta_opt);
                if (f1 >= f_opt && f1 >= f3) theta_opt = theta_A;
                else if (f3 >= f_opt && f3 >= f1) theta_opt = theta_C;

                theta[n] = theta_opt;
                v[n] = polar(beta_model(theta[n]), theta[n]);
            } 
            else if (method_id == 3) { 
                // Practical IRS (AO with 1D Search)
                double lambda = (arg_phi >= 0) ? 0 : 1;
                double theta_A = arg_phi;
                double theta_C = (lambda == 0) ? PI : -PI;
                double min_th = min(theta_A, theta_C);
                double max_th = max(theta_A, theta_C);

                double best_th = min_th;
                double best_f = -1e9;
                for (int step = 0; step <= 100; step++) {
                    double th = min_th + step * (max_th - min_th) / 100.0;
                    double b = beta_model(th);
                    double val = b * b * psi_nn + b * abs_phi * cos(arg_phi - th);
                    if (val > best_f) { best_f = val; best_th = th; }
                }
                theta[n] = best_th;
                v[n] = polar(beta_model(theta[n]), theta[n]);
            }
        }
    }

    // Evaluate lại cho method 4 (Practical IRS with Ideal IRS Assumption)
    if (method_id == 4) {
        for (int n = 0; n < N_ELEMENTS; n++) {
            v[n] = polar(beta_model(theta[n]), theta[n]);
        }
    }

    return calculate_baseline_rate(v, Phi_c, hd_c);
}

int main() { 
    cout << "====================================================================\n";
    cout << "           SURVEY DISTANCE (COMPARING 6 ALGORITHMS)                 \n";
    cout << "====================================================================\n";
    for (double d_AP_USER = 480.0; d_AP_USER <= 500.0; d_AP_USER += 5.0) { 
        init_channels(d_AP_USER); // KHỞI TẠO CHUNG: Tất cả giải thuật dùng chung kênh
        
        cout << "[ DISTANCE = " << (int)d_AP_USER << "m ]\n";
        
        // 1. ES (Tiến hóa trực tiếp trên linh kiện vật lý L, R, C)
        double rate_es = run_evolution_strategy();
        
        // 2-6. 5 Đường Baseline bài báo gốc
        double rate_ideal          = run_baseline(1); // Upper Bound
        double rate_prop1          = run_baseline(2); // AO Prop 1
        double rate_1d             = run_baseline(3); // AO 1D Search
        double rate_ideal_assum    = run_baseline(4); // Ideal Assumption
        double rate_no_irs         = run_baseline(5); // Lower Bound

        cout << fixed << setprecision(4);
        cout << "  -> (ES) Proposed Physical Opt   : " << rate_es          << " bps/Hz\n";
        cout << "  -> (1)  Upper Bound (Ideal IRS) : " << rate_ideal       << " bps/Hz\n";
        cout << "  -> (2)  AO 1D Search            : " << rate_1d          << " bps/Hz\n";
        cout << "  -> (3)  AO Proposition 1        : " << rate_prop1       << " bps/Hz\n";
        cout << "  -> (4)  Ideal IRS Assumption    : " << rate_ideal_assum << " bps/Hz\n";
        cout << "  -> (5)  Lower Bound (No IRS)    : " << rate_no_irs      << " bps/Hz\n";
        cout << "--------------------------------------------------------------------\n";
    }

    cout << "\n====================================================================\n";
    cout << "           SURVEY N_ELEMENTS (COMPARING 6 ALGORITHMS)               \n";
    cout << "====================================================================\n";
    for (int current_N = 10; current_N <= 50; current_N += 10) { 
        N_ELEMENTS = current_N;
        init_channels(498.0); // KHỞI TẠO CHUNG
        
        cout << "[ N_ELEMENTS = " << N_ELEMENTS << " ]\n";
        
        double rate_es = run_evolution_strategy();
        
        double rate_ideal          = run_baseline(1);
        double rate_prop1          = run_baseline(2);
        double rate_1d             = run_baseline(3);
        double rate_ideal_assum    = run_baseline(4);
        double rate_no_irs         = run_baseline(5);

        cout << fixed << setprecision(4);
        cout << "  -> (ES) Proposed Physical Opt   : " << rate_es          << " bps/Hz\n";
        cout << "  -> (1)  Upper Bound (Ideal IRS) : " << rate_ideal       << " bps/Hz\n";
        cout << "  -> (2)  AO 1D Search            : " << rate_1d          << " bps/Hz\n";
        cout << "  -> (3)  AO Proposition 1        : " << rate_prop1       << " bps/Hz\n";
        cout << "  -> (4)  Ideal IRS Assumption    : " << rate_ideal_assum << " bps/Hz\n";
        cout << "  -> (5)  Lower Bound (No IRS)    : " << rate_no_irs      << " bps/Hz\n";
        cout << "--------------------------------------------------------------------\n";
    }
    
    return 0;
}
