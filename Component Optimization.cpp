#include <bits/stdc++.h>
using namespace std;
#define ll long long
#define ld long double

const ld BETA_MIN = 0.2; // 0.2
const ld K = 1.6; // 1.6
const ld noise = pow(10.0, (-94 - 30.0) / 10.0);
const ld Z0 = 377.0;    // trở kháng không gian tự do
const ld PT = pow(10.0, (36 - 30.0) / 10.0);
const ld omega = 2 * M_PI * 2.4e9;
const ll SEARCH_STEPS = 360;
const ll population = 100;

random_device rd;
mt19937 gen(rd());
uniform_real_distribution<> phase_gacha(-3.14159265358979323846, 3.14159265358979323846);
uniform_real_distribution<> percen_gacha(0.0, 1.0); 
normal_distribution<> rayleigh(0.0, 1.0 / sqrt(2.0));
struct comp // A * e^(i * theta) = A(cos(theta) + i * sin(theta))
{
    ld A; // amplitude
    ld theta; // phase
}; 

struct standardComp // a + bi
{
    ld real;
    ld img;
};

standardComp convertFromEuler(comp a)
{
    return {a.A * cos(a.theta), a.A * sin(a.theta)};
}

comp convertToEuler(standardComp a)
{
    ld A = sqrt(a.real * a.real + a.img * a.img);
    ld theta = atan2(a.img, a.real);
    if (theta < 0) theta += 2 * M_PI; 

    return {A, theta};
}

struct Individual
{
    ld C, R, L1, L2;
    ld minC, maxC;
    ld minR, maxR;
    ld minL1, maxL1;
    ld minL2, maxL2;

    void set_constraint_C(ld mint, ld maxt)
    {   
        minC = mint;
        maxC = maxt;
    }
    void set_constraint_R(ld mint, ld maxt)
    {   
        minR = mint;
        maxR = maxt;
    }
    void set_constraint_L1(ld mint, ld maxt)
    {   
        minL1 = mint;
        maxL1 = maxt;
    }
    void set_constraint_L2(ld mint, ld maxt)
    {   
        minL2 = mint;
        maxL2 = maxt;
    }

    void initC()
    {
        uniform_real_distribution<> tmp(minC, maxC);
        C = tmp(gen);
    }
    void initR()
    {
        uniform_real_distribution<> tmp(minR, maxR);
        R = tmp(gen);
    }
    void initL1()
    {
        uniform_real_distribution<> tmp(minL1, maxL1);
        L1 = tmp(gen);
    }
    void initL2()
    {
        uniform_real_distribution<> tmp(minL2, maxL2);
        L2 = tmp(gen);
    }
};



struct Matrix 
{
    ll row, col;
	vector<vector<comp> > a;
	
	void prep(ll n, ll m)
	{
		a.resize(n, vector<comp>(m));
		row = n; col = m;
		for (ll i = 0; i < n; ++i)
		{
			for (ll j = 0; j < m; ++j)
			{
				a[i][j] = {0, 0};
			}
		}
	}
	
    void show()
    {
        for (ll i = 0; i < row; ++i)
        {
            for (ll j = 0; j < col; ++j)
            {
                cout << "[" << a[i][j].A << "; " << a[i][j].theta << "]";
            }
            cout << '\n';
        }
        
    }

    void diag()
    {
        if (row != 1 && col != 1) return; 
        
        ll size = max(row, col);
        vector<vector<comp> > tmp(size, vector<comp>(size, {0, 0}));

        for (ll i = 0; i < size; ++i)
        {
            tmp[i][i] = (row == 1) ? a[0][i] : a[i][0];
        }
        
        a = tmp;
        row = size;
        col = size;
    }

    void transpose()
    {
        vector<vector<comp> > tmp;
        tmp.resize(col, vector<comp>(row));

        for (ll i = 0; i < row; ++i)
        {
            for (ll j = 0; j < col; ++j)
            {
                tmp[j][i] = a[i][j]; 
            }
        }
        a.clear();
        a.resize(col, vector<comp>(row));
        swap(row, col);
        for (ll i = 0; i < row; ++i)
        {
            for (ll j = 0; j < col; ++j)
            {
                a[i][j] = tmp[i][j];
            }
        }
    }

    void hermitian()
    {
        vector<vector<comp> > tmp;
        tmp.resize(col, vector<comp>(row));

        for (ll i = 0; i < row; ++i)
        {
            for (ll j = 0; j < col; ++j)
            {
                tmp[j][i] = a[i][j]; 
            }
        }
        a.clear();
        a.resize(col, vector<comp>(row));
        swap(row, col);
        for (ll i = 0; i < row; ++i)
        {
            for (ll j = 0; j < col; ++j)
            {
                a[i][j] = tmp[i][j];
                a[i][j].theta *= -1;
            }
        }
    }
};

comp sumComplex(comp a, comp b)
{
    ld real = a.A * cos(a.theta) + b.A * cos(b.theta);
    ld imag = a.A * sin(a.theta) + b.A * sin(b.theta);
    
    ld A = sqrt(real * real + imag * imag);
    ld theta = atan2(imag, real); 
    
    if (theta < 0) theta += 2 * M_PI; 
    
    return {A, theta};
}

comp multComplex(comp a, comp b)
{
    return {a.A * b.A, a.theta + b.theta};
}

comp divComplex(comp a, comp b)
{
    return {a.A / b.A, a.theta - b.theta};
}

Matrix sum(Matrix A, Matrix B)
{
    Matrix tmp2;
    tmp2.prep(A.row, A.col);

    for (ll i = 0; i < A.row; ++i)
    {
        for (ll j = 0; j < A.col; ++j)
        {
            tmp2.a[i][j] = sumComplex(A.a[i][j], B.a[i][j]);
        }
    }

    return tmp2;
}
Matrix mult(Matrix A, Matrix B)
{
    Matrix tmp2;
    tmp2.prep(A.row, B.col);
	for (ll i = 0; i < A.row; ++i)
	{
		for (ll k = 0; k < B.col; ++k)
		{
			for (ll j = 0; j < A.col; ++j)
			{
				tmp2.a[i][k] = sumComplex(tmp2.a[i][k], multComplex(A.a[i][j], B.a[j][k]));
			}
		}
	}
	return tmp2;
}


ld calculateRate(Matrix& v, Matrix& h_r, Matrix& G, Matrix& h_d, ld noise_power)
{
    Matrix Hv = v; Hv.hermitian();
    Matrix Hh_d = h_d; Hh_d.hermitian();

    Matrix phi = h_r;
    phi.hermitian();
    phi.diag();
    phi = mult(phi, G);

    // h_eff^H = v^H * Phi + h_d^H
    Matrix h_eff_H = sum(mult(Hv, phi), Hh_d);
    
    //  P_T * ||h_eff||^2
    ld norm_sq = 0.0L;
    for (ll i = 0; i < h_eff_H.col; ++i) 
    {
        norm_sq += h_eff_H.a[0][i].A * h_eff_H.a[0][i].A;
    }

    return log2(1.0L + ((PT * norm_sq) / noise_power));
}

ld getPathLoss(ld distance, ld exponent) 
{
    ld ref_loss_db = 40.0; 
    ld ref_loss_linear = pow(10, -ref_loss_db / 10.0);
    return ref_loss_linear * pow(distance, -exponent);
}

ld calculateNoIRSRate(Matrix& h_d, ld noise_power)
{
    ld norm_sq = 0.0L;
    for (ll i = 0; i < h_d.row; ++i) 
    {
        norm_sq += h_d.a[i][0].A * h_d.a[i][0].A;
    }
    return log2(1.0L + ((PT * norm_sq) / noise_power));
}

comp calculateReflectionCoefficient(ld C, ld R, ld L1, ld L2)
{
    ld X1 = omega * L1;
    ld X2 = omega * L2 - 1.0L / (omega * C);

    comp Z1 = {abs(X1), (X1 >= 0.0L ? (ld)M_PI / 2.0L : (ld)-M_PI / 2.0L)};
    
    comp Z2 = convertToEuler({R, X2});
    comp Z_num = multComplex(Z1, Z2);
    comp Z_den = convertToEuler({R, X1 + X2});
    comp Z = divComplex(Z_num, Z_den);
    standardComp Z_rect = convertFromEuler(Z);
    
    standardComp v_num_rect = {Z_rect.real - Z0, Z_rect.img};
    standardComp v_den_rect = {Z_rect.real + Z0, Z_rect.img};

    comp v_num = convertToEuler(v_num_rect);
    comp v_den = convertToEuler(v_den_rect);

    return divComplex(v_num, v_den);
}

struct IRSElement 
{
    ld C, R, L1, L2;
};

struct Solution 
{
    vector<IRSElement> elements;
    ld fitness; 
};

ld clamp_val(ld val, ld min_val, ld max_val) 
{
    if (val < min_val) return min_val;
    if (val > max_val) return max_val;
    return val;
}

int main()
{
    ll n = 40; // IRS elements
    ll m = 2;  // AP antennas

    ld C_min, C_max, R_min, R_max, L1_min, L1_max, L2_min, L2_max;
    cout << "Constraint for C: ";
    cin >> C_min >> C_max;
    cout << "Constraint for R: ";
    cin >> R_min >> R_max;
    cout << "Constraint for L1: ";
    cin >> L1_min >> L1_max;
    cout << "Constraint for L2: ";
    cin >> L2_min >> L2_max;

    vector<ld> dist;
    vector<ld> AR;
    for (ld d = 480; d <= 500; d += 5) 
    //for (ll n = 10;  n <= 50; ++n)
    {
        //ld d = 498;
        Matrix h_d, h_r, G; 
        h_d.prep(m, 1); h_r.prep(n, 1); G.prep(n, m);

        
        ld d_ap_irs = 500.0; 
        ld d_ap_user = sqrt(d * d + 2.0 * 2.0);
        ld d_irs_user = sqrt((500.0 - d) * (500.0 - d) + 2.0 * 2.0);
        
        ld h_d_gain = sqrt(getPathLoss(d_ap_user, 3.8));
        ld G_gain = sqrt(getPathLoss(d_ap_irs, 2.2));
        ld h_r_gain = sqrt(getPathLoss(d_irs_user, 2.8));

        ld total_rate = 0.0L;
        ll num_trials = 1; // monte-carlo

        for (ll trial = 0; trial < num_trials; ++trial)
        {
            for (ll i = 0; i < h_d.row; ++i) 
            {
                ld real = rayleigh(gen) * h_d_gain;
                ld imag = rayleigh(gen) * h_d_gain;
                h_d.a[i][0] = {sqrt(real*real + imag*imag), atan2(imag, real)};
            }

            for (ll i = 0; i < h_r.row; ++i) 
            {
                ld real = rayleigh(gen) * h_r_gain;
                ld imag = rayleigh(gen) * h_r_gain;
                h_r.a[i][0] = {sqrt(real*real + imag*imag), atan2(imag, real)};
            }

            for (ll i = 0; i < G.row; ++i) 
            {
                for (ll j = 0; j < G.col; ++j) 
                {
                    ld real = rayleigh(gen) * G_gain;
                    ld imag = rayleigh(gen) * G_gain;
                    G.a[i][j] = {sqrt(real*real + imag*imag), atan2(imag, real)};
                }
            }

            // end setup

            vector<Individual> component(population);
            
            

            for (ll i = 0; i < population; ++i)
            {
                Individual tmp;
                tmp.set_constraint_C(C_min, C_max);
                tmp.set_constraint_R(R_min, R_max);
                tmp.set_constraint_L1(L1_min, L1_max);
                tmp.set_constraint_L2(L2_min, L2_max);

                tmp.initC();
                tmp.initR();
                tmp.initL1();
                tmp.initL2();
                component.push_back(tmp);
                
            }   

            for (ll i = 0; i < population; ++i)
            {
                Matrix v; 
                v.prep(n, 1);
                for (ll i = 0; i < n; ++i)

                v.a[i][0] = calculateReflectionCoefficient(component[i].C, 
                                                            component[i].R, 
                                                            component[i].L1, 
                                                            component[i].L2);
                
            }

            ll mu = population / 4; 
            
            // standard mutation, can change
            ld std_C = (C_max - C_min) * 0.05;
            ld std_R = (R_max - R_min) * 0.05;
            ld std_L1 = (L1_max - L1_min) * 0.05;
            ld std_L2 = (L2_max - L2_min) * 0.05;

            uniform_real_distribution<ld> initC(C_min, C_max);
            uniform_real_distribution<ld> initR(R_min, R_max);
            uniform_real_distribution<ld> initL1(L1_min, L1_max);
            uniform_real_distribution<ld> initL2(L2_min, L2_max);
            
            uniform_int_distribution<ll> parent_selector(0, mu - 1);

            vector<Solution> pop(population);

            for (ll p = 0; p < population; ++p) 
            {
                pop[p].elements.resize(n);
                Matrix v; 
                v.prep(n, 1);

                for (ll i = 0; i < n; ++i) 
                {
                    pop[p].elements[i] = {initC(gen), initR(gen), initL1(gen), initL2(gen)};
                    v.a[i][0] = calculateReflectionCoefficient(pop[p].elements[i].C, 
                                                               pop[p].elements[i].R, 
                                                               pop[p].elements[i].L1, 
                                                               pop[p].elements[i].L2);
                }
                pop[p].fitness = calculateRate(v, h_r, G, h_d, noise);
            }

            for (ll gen_idx = 0; gen_idx < SEARCH_STEPS; ++gen_idx) 
            {
                sort(pop.begin(), pop.end(), [](const Solution& a, const Solution& b) {
                    return a.fitness > b.fitness;
                });

                for (ll p = mu; p < population; ++p) 
                {
                    ll p_idx = parent_selector(gen);
                    Solution offspring = pop[p_idx];

                    Matrix v; 
                    v.prep(n, 1);

                    for (ll i = 0; i < n; ++i) 
                    {
                        normal_distribution<ld> mutC(0, std_C);
                        normal_distribution<ld> mutR(0, std_R);
                        normal_distribution<ld> mutL1(0, std_L1);
                        normal_distribution<ld> mutL2(0, std_L2);

                        offspring.elements[i].C  = clamp_val(offspring.elements[i].C + mutC(gen), C_min, C_max);
                        offspring.elements[i].R  = clamp_val(offspring.elements[i].R + mutR(gen), R_min, R_max);
                        offspring.elements[i].L1 = clamp_val(offspring.elements[i].L1 + mutL1(gen), L1_min, L1_max);
                        offspring.elements[i].L2 = clamp_val(offspring.elements[i].L2 + mutL2(gen), L2_min, L2_max);

                        v.a[i][0] = calculateReflectionCoefficient(offspring.elements[i].C, 
                                                                   offspring.elements[i].R, 
                                                                   offspring.elements[i].L1, 
                                                                   offspring.elements[i].L2);
                    }
                    
                    offspring.fitness = calculateRate(v, h_r, G, h_d, noise);
                    pop[p] = offspring;
                }
            }

            sort(pop.begin(), pop.end(), [](const Solution& a, const Solution& b) {
                return a.fitness > b.fitness;
            });

            cout << "Distance " << d << ": " << pop[0].fitness << '\n';
        }

    

}