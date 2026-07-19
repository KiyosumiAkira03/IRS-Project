#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <random>
#include <algorithm>
#include <string>
#include <iomanip>
#include <future>
#include <atomic>

using namespace std;

// ==============================================================================
// 1. CẤU HÌNH BÀI TOÁN VÀ THÔNG SỐ VẬT LÝ
// ==============================================================================
int N_ELEMENTS = 40;      // Số phần tử IRS (N) - Sẽ thay đổi ở vòng lặp
const int M_ANTENNAS = 2; // Số anten tại AP (M)
const int POPULATION_SIZE = 100;
const int GENERATIONS = 1000;
const double MUTATION_RATE = 0.15;
const double PI = 3.14159265358979323846;

const double PT = pow(10.0, (36 - 30.0) / 10.0);
const double NOISE_POWER = pow(10.0, (-94 - 30.0) / 10.0);

const double W = 2.0 * PI * 2.4e9; // Tần số góc omega = 2*PI*2.4GHz
const double Z0 = 377.0;           // Trở kháng không gian tự do (Ohm)

// Ràng buộc vật lý của linh kiện
const double C_MIN = 0.466e-12, C_MAX = 2.35e-12; // Farad
const double R_MIN = 0.47,      R_MAX = 30.8;     // Ohm
const double L1_MIN = 0.5e-9,   L1_MAX = 13.0e-9; // Henry
const double L2_MIN = 0.5e-9,   L2_MAX = 13.0e-9; // Henry

// Các giá trị cố định chuẩn mực (Fixed Values)
const double FIXED_C_VAL = 1.408e-12; 
const double FIXED_R_VAL = 2.0;       
const double FIXED_L1_VAL = 2.5e-9;   
const double FIXED_L2_VAL = 0.7e-9;   

// Thông số của Practical Phase Shift Model (Bài báo baseline)
const double BETA_MIN = 0.2;
const double MODEL_K = 1.6;
const double MODEL_PHI = 0.43 * PI;

// Cấu trúc đánh dấu việc cố định linh kiện cho 16 trường hợp
struct FixConfig {
    bool fix_C, fix_R, fix_L1, fix_L2;
    string name;
    int fixed_count;
};

// ==============================================================================
// 2. CÁC CẤU TRÚC DỮ LIỆU & HÀM PHỤ TRỢ SỐ PHỨC
// ==============================================================================
struct comp { 
    double A, theta; 
};

comp sumComplex(comp a, comp b) {
    double real = a.A * cos(a.theta) + b.A * cos(b.theta);
    double imag = a.A * sin(a.theta) + b.A * sin(b.theta);
    return {sqrt(real * real + imag * imag), atan2(imag, real)};
}

comp multComplex(comp a, comp b) {
    return {a.A * b.A, a.theta + b.theta};
}

struct IRSElement { 
    double C, R, L1, L2; 
};

struct Individual {
    vector<IRSElement> elements; 
    double fitness;              
    bool operator>(const Individual& other) const { return fitness > other.fitness; }
};

struct ESResult {
    double rate;
    Individual best_ind;
};

double clamp(double val, double min_val, double max_val) {
    return max(min_val, min(max_val, val));
}

// Tính toán hệ số phản xạ từ mô hình mạch điện tương đương
comp calculate_vn_from_physics(const IRSElement& el) {
    complex<double> j(0.0, 1.0);
    complex<double> Z_L1 = j * W * el.L1;
    complex<double> Z_L2 = j * W * el.L2;
    complex<double> Z_C = 1.0 / (j * W * el.C);
    
    complex<double> branch2 = Z_L2 + Z_C + el.R;
    complex<double> Zn = (Z_L1 * branch2) / (Z_L1 + branch2);
    
    complex<double> vn = (Zn - Z0) / (Zn + Z0);
    return {abs(vn), arg(vn)};
}

// ==============================================================================
// 3. KHỞI TẠO KÊNH TRUYỀN (DÙNG CHUNG)
// ==============================================================================
vector<comp> h_d; 
vector<vector<comp>> Phi; 

const double ALPHA_AP_USER = 3.8, ALPHA_IRS_USER = 2.8, ALPHA_AP_IRS = 2.2;
const double C0 = 1e-4;
const double d_AP_IRS = 500.0;

random_device rd;
mt19937 global_gen(rd());

// Khởi tạo kênh truyền dùng chung cho TẤT CẢ các trường hợp trong 1 bước khảo sát
void init_channels(double d_AP_USER) {
    uniform_real_distribution<> rayleigh(0.5, 1.5);
    uniform_real_distribution<> phase_dist(-PI, PI);
    
    double d_IRS_USER = abs(d_AP_IRS - d_AP_USER);
    if (d_IRS_USER < 1.0) d_IRS_USER = 1.0; 

    double d_AP_USER_true = sqrt(d_AP_USER * d_AP_USER + 4.0);
    double d_IRS_USER_true = sqrt((d_AP_IRS - d_AP_USER) * (d_AP_IRS - d_AP_USER) + 4.0);

    double path_loss_direct = sqrt(C0 * pow(d_AP_USER_true, -ALPHA_AP_USER));
    double path_loss_cascaded = sqrt(C0 * pow(d_AP_IRS, -ALPHA_AP_IRS)) * sqrt(C0 * pow(d_IRS_USER_true, -ALPHA_IRS_USER));
    
    h_d.resize(M_ANTENNAS);
    for (int m = 0; m < M_ANTENNAS; ++m) 
        h_d[m] = {path_loss_direct * rayleigh(global_gen), phase_dist(global_gen)};

    Phi.assign(N_ELEMENTS, vector<comp>(M_ANTENNAS));
    for (int n = 0; n < N_ELEMENTS; ++n) {
        for (int m = 0; m < M_ANTENNAS; ++m) {
            Phi[n][m] = {path_loss_cascaded * rayleigh(global_gen), phase_dist(global_gen)};
        }
    }
}

// Đánh giá Fitness dựa trên kênh truyền dùng chung
void evaluate(Individual& ind) {
    double total_norm_sq = 0;
    for (int m = 0; m < M_ANTENNAS; ++m) {
        comp sum_m = {0, 0};
        for (int n = 0; n < N_ELEMENTS; ++n) {
            comp v_n = calculate_vn_from_physics(ind.elements[n]);
            sum_m = sumComplex(sum_m, multComplex({v_n.A, -v_n.theta}, Phi[n][m]));
        }
        sum_m = sumComplex(sum_m, {h_d[m].A, -h_d[m].theta});
        total_norm_sq += (sum_m.A * sum_m.A);
    }
    ind.fitness = log2(1.0 + (total_norm_sq * PT) / NOISE_POWER);
}

// ==============================================================================
// 4. THUẬT TOÁN EVOLUTION STRATEGY (TIẾN HÓA ĐỘC LẬP TỪ GỐC CHUNG)
// ==============================================================================
ESResult run_evolution_strategy(FixConfig config, vector<Individual> base_population, unsigned int seed) {
    // Sử dụng seed độc lập cho luồng này để quá trình tiến hóa không bị gò bó
    mt19937 local_gen(seed);
    uniform_real_distribution<> prob_dist(0.0, 1.0);
    normal_distribution<> norm_dist(0.0, 1.0);
    uniform_int_distribution<> parent_dist(0, 49);

    vector<Individual> population = base_population;
    
    // ÁP DỤNG RÀNG BUỘC: Đè giá trị cố định lên cá thể, giữ nguyên các thông số tự do
    for (auto& ind : population) {
        for (auto& el : ind.elements) {
            if (config.fix_C)  el.C  = FIXED_C_VAL;
            if (config.fix_R)  el.R  = FIXED_R_VAL;
            if (config.fix_L1) el.L1 = FIXED_L1_VAL;
            if (config.fix_L2) el.L2 = FIXED_L2_VAL;
        }
        evaluate(ind); 
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
            Individual p1 = population[parent_dist(local_gen)];
            Individual p2 = population[parent_dist(local_gen)];

            Individual child;
            child.elements.resize(N_ELEMENTS);
            for (int n = 0; n < N_ELEMENTS; ++n) {
                // Lai ghép
                child.elements[n].C  = config.fix_C  ? FIXED_C_VAL  : ((prob_dist(local_gen) < 0.5) ? p1.elements[n].C  : p2.elements[n].C);
                child.elements[n].R  = config.fix_R  ? FIXED_R_VAL  : ((prob_dist(local_gen) < 0.5) ? p1.elements[n].R  : p2.elements[n].R);
                child.elements[n].L1 = config.fix_L1 ? FIXED_L1_VAL : ((prob_dist(local_gen) < 0.5) ? p1.elements[n].L1 : p2.elements[n].L1);
                child.elements[n].L2 = config.fix_L2 ? FIXED_L2_VAL : ((prob_dist(local_gen) < 0.5) ? p1.elements[n].L2 : p2.elements[n].L2);
                
                // Đột biến
                if (prob_dist(local_gen) < MUTATION_RATE) {
                    if (!config.fix_C)  child.elements[n].C  = clamp(child.elements[n].C  + norm_dist(local_gen) * sigma * (C_MAX - C_MIN), C_MIN, C_MAX);
                    if (!config.fix_R)  child.elements[n].R  = clamp(child.elements[n].R  + norm_dist(local_gen) * sigma * (R_MAX - R_MIN), R_MIN, R_MAX);
                    if (!config.fix_L1) child.elements[n].L1 = clamp(child.elements[n].L1 + norm_dist(local_gen) * sigma * (L1_MAX - L1_MIN), L1_MIN, L1_MAX);
                    if (!config.fix_L2) child.elements[n].L2 = clamp(child.elements[n].L2 + norm_dist(local_gen) * sigma * (L2_MAX - L2_MIN), L2_MIN, L2_MAX);
                }
            }
            evaluate(child);
            total_children++;
            if (child.fitness > max(p1.fitness, p2.fitness)) success_count++;
            next_gen.push_back(child);
        }
        population = next_gen;

        // Cập nhật sigma theo nguyên tắc 1/5 (1/5th success rule)
        if (total_children > 0) {
            double success_rate = (double)success_count / total_children;
            if (success_rate > 0.2) sigma /= 0.85; 
            else if (success_rate < 0.2) sigma *= 0.85; 
            sigma = clamp(sigma, 0.001, 1.0); 
        }
    }
    return {best_overall.fitness, best_overall};
}

// Hàm sinh 16 trường hợp cố định theo Bitmask
vector<FixConfig> generate_16_configs() {
    vector<FixConfig> configs;
    for(int i = 0; i < 16; ++i) {
        FixConfig cfg;
        cfg.fix_C  = (i & 1); cfg.fix_R  = (i & 2); cfg.fix_L1 = (i & 4); cfg.fix_L2 = (i & 8);
        cfg.fixed_count = cfg.fix_C + cfg.fix_R + cfg.fix_L1 + cfg.fix_L2;
        
        if (cfg.fixed_count == 0) cfg.name = "[0 Fix] None";
        else if (cfg.fixed_count == 4) cfg.name = "[4 Fix] All";
        else {
            cfg.name = "[" + to_string(cfg.fixed_count) + " Fix] ";
            if(cfg.fix_C) cfg.name += "C "; if(cfg.fix_R) cfg.name += "R ";
            if(cfg.fix_L1) cfg.name += "L1 "; if(cfg.fix_L2) cfg.name += "L2 ";
        }
        configs.push_back(cfg);
    }
    // Sắp xếp tăng dần theo số lượng thông số bị cố định
    sort(configs.begin(), configs.end(), [](const FixConfig& a, const FixConfig& b) {
        if(a.fixed_count != b.fixed_count) return a.fixed_count < b.fixed_count;
        return a.name < b.name;
    });
    return configs;
}

// ==============================================================================
// 5. CÁC HÀM BASELINES BÀI BÁO 
// ==============================================================================
double beta_model(double theta) {
    return (1.0 - BETA_MIN) * pow((sin(theta - MODEL_PHI) + 1.0) / 2.0, MODEL_K) + BETA_MIN;
}

double calculate_baseline_rate(const vector<complex<double>>& v, const vector<vector<complex<double>>>& Phi_c, const vector<complex<double>>& hd_c) {
    double total_norm_sq = 0;
    for (int m = 0; m < M_ANTENNAS; ++m) {
        complex<double> sum_m = 0;
        for (int n = 0; n < N_ELEMENTS; ++n) sum_m += conj(v[n]) * Phi_c[n][m];
        sum_m += conj(hd_c[m]);
        total_norm_sq += norm(sum_m);
    }
    double snr = (total_norm_sq * PT) / NOISE_POWER;
    return log2(1.0 + snr);
}

double run_baseline(int method_id) {
    uniform_real_distribution<> prob_dist(0.0, 1.0);
    vector<vector<complex<double>>> Phi_c(N_ELEMENTS, vector<complex<double>>(M_ANTENNAS));
    vector<complex<double>> hd_c(M_ANTENNAS);

    for (int n = 0; n < N_ELEMENTS; n++) {
        for (int m = 0; m < M_ANTENNAS; m++) Phi_c[n][m] = polar(Phi[n][m].A, Phi[n][m].theta);
    }
    for (int m = 0; m < M_ANTENNAS; m++) hd_c[m] = polar(h_d[m].A, h_d[m].theta);

    if (method_id == 5) {
        vector<complex<double>> v(N_ELEMENTS, 0.0);
        return calculate_baseline_rate(v, Phi_c, hd_c);
    }

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

    vector<double> theta(N_ELEMENTS);
    vector<complex<double>> v(N_ELEMENTS);
    for (int n = 0; n < N_ELEMENTS; n++) {
        theta[n] = (prob_dist(global_gen) > 0.5) ? PI : -PI;
        v[n] = polar(beta_model(theta[n]), theta[n]);
    }

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
                theta[n] = arg_phi;
                v[n] = polar(1.0, theta[n]); 
            } 
            else if (method_id == 2) { 
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
                theta_opt = clamp(theta_opt, min(theta_A, theta_C), max(theta_A, theta_C));

                double f_opt = obj_func(theta_opt);
                if (f1 >= f_opt && f1 >= f3) theta_opt = theta_A;
                else if (f3 >= f_opt && f3 >= f1) theta_opt = theta_C;

                theta[n] = theta_opt;
                v[n] = polar(beta_model(theta[n]), theta[n]);
            } 
            else if (method_id == 3) { 
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

    if (method_id == 4) {
        for (int n = 0; n < N_ELEMENTS; n++) v[n] = polar(beta_model(theta[n]), theta[n]);
    }
    return calculate_baseline_rate(v, Phi_c, hd_c);
}

// ==============================================================================
// 6. LUỒNG ĐIỀU KHIỂN & IN KẾT QUẢ KHẢO SÁT
// ==============================================================================
void print_survey_step() {
    auto configs = generate_16_configs();
    
    // --- BƯỚC 1: SINH QUẦN THỂ GỐC CHUNG MỘT LẦN DUY NHẤT ---
    uniform_real_distribution<> C_dist(C_MIN, C_MAX);
    uniform_real_distribution<> R_dist(R_MIN, R_MAX);
    uniform_real_distribution<> L1_dist(L1_MIN, L1_MAX);
    uniform_real_distribution<> L2_dist(L2_MIN, L2_MAX);

    vector<Individual> base_population(POPULATION_SIZE);
    for (int i = 0; i < POPULATION_SIZE; ++i) {
        base_population[i].elements.resize(N_ELEMENTS);
        for (int n = 0; n < N_ELEMENTS; ++n) {
            base_population[i].elements[n] = {
                C_dist(global_gen), R_dist(global_gen), 
                L1_dist(global_gen), L2_dist(global_gen)
            };
        }
    }

    // --- BƯỚC 2: CHẠY SONG SONG 16 CẤU HÌNH ---
    vector<future<ESResult>> futures;
    unsigned int base_seed = global_gen(); // Lấy một seed gốc để chia cho các luồng
    
    for(size_t i = 0; i < configs.size(); ++i) {
        // Truyền base_population giống hệt nhau, nhưng seed khác nhau (base_seed + i)
        // để thuật toán tự do tiến hóa theo bề mặt fitness riêng của nó
        futures.push_back(async(launch::async, run_evolution_strategy, configs[i], base_population, base_seed + i));
    }
    
    // --- BƯỚC 3: IN KẾT QUẢ ---
    cout << "  >>> 16 TRƯỜNG HỢP CỐ ĐỊNH ES (Thông số bộ nghiệm của phần tử IRS số 0) <<<\n";
    for(size_t i = 0; i < configs.size(); ++i) {
        ESResult result = futures[i].get(); 
        IRSElement el0 = result.best_ind.elements[0]; // Chỉ lấy phần tử đầu tiên
        
        cout << "    " << left << setw(15) << configs[i].name 
             << " | Rate: " << fixed << setprecision(4) << result.rate << " bps/Hz | "
             << "L1: " << right << setw(5) << setprecision(2) << el0.L1 * 1e9 << " nH, "
             << "L2: " << right << setw(4) << setprecision(2) << el0.L2 * 1e9 << " nH, "
             << "C: "  << right << setw(5) << setprecision(3) << el0.C * 1e12 << " pF, "
             << "R: "  << right << setw(5) << setprecision(2) << el0.R << " Ohm\n";
    }

    cout << "\n  >>> 5 ĐƯỜNG BASELINE BÀI BÁO GỐC <<<\n";
    cout << "    (1) Upper Bound (Ideal IRS) : " << run_baseline(1) << " bps/Hz\n";
    cout << "    (2) AO 1D Search            : " << run_baseline(3) << " bps/Hz\n";
    cout << "    (3) AO Proposition 1        : " << run_baseline(2) << " bps/Hz\n";
    cout << "    (4) Ideal IRS Assumption    : " << run_baseline(4) << " bps/Hz\n";
    cout << "    (5) Lower Bound (No IRS)    : " << run_baseline(5) << " bps/Hz\n";
    cout << "--------------------------------------------------------------------------------------\n";
}

int main() { 
    cout << "======================================================================================\n";
    cout << "               SURVEY DISTANCE (COMPARING MULTIPLE CONFIGURATIONS)                    \n";
    cout << "======================================================================================\n";
    for (double d_AP_USER = 480.0; d_AP_USER <= 500.0; d_AP_USER += 5.0) { 
        init_channels(d_AP_USER); // Khởi tạo kênh chung 1 lần cho cả cụm
        cout << "[ DISTANCE = " << (int)d_AP_USER << "m ]\n";
        print_survey_step();
    }

    cout << "\n======================================================================================\n";
    cout << "               SURVEY N_ELEMENTS (COMPARING MULTIPLE CONFIGURATIONS)                  \n";
    cout << "======================================================================================\n";
    for (int current_N = 10; current_N <= 50; current_N += 10) { 
        N_ELEMENTS = current_N;
        init_channels(498.0); // Khởi tạo kênh chung (khoảng cách cố định ở 498m hoặc tùy bạn)
        cout << "[ N_ELEMENTS = " << N_ELEMENTS << " ]\n";
        print_survey_step();
    }
    
    return 0;
}
