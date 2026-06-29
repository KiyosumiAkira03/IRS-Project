#include <bits/stdc++.h>
using namespace std;
#define ll long long
#define ld long double

const ld BETA_MIN = 0.2; // 0.2
const ld K = 1.6; // 1.6
const ld noise = pow(10.0, (-94 - 30.0) / 10.0);
const ld L1 = 2.5e-9;   // nH
const ld L2 = 0.7e-9;   // nH 
const ld Z0 = 377.0;    // ohm
const ld freq = 2.4e9;  // GHz
const ld PT = pow(10.0, (36 - 30.0) / 10.0);
const int SEARCH_STEPS = 360;

random_device rd;
mt19937 gen(rd());
uniform_real_distribution<> phase_gacha(-3.14159265358979323846, 3.14159265358979323846);
uniform_real_distribution<> percen_gacha(0.0, 1.0); 
uniform_real_distribution<> rayleigh(0.5, 1.5); 

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


ld beta(ld theta) 
{
    // beta(theta) = (1 - beta_min) * ((sin(theta - phi) + 1) / 2)^k + beta_min
    ld base = (sin(theta - 0.43 * M_PI) + 1.0) / 2.0; // 0.43 * M_PI
    return (1.0 - BETA_MIN) * pow(base, K) + BETA_MIN;
}

ld f(comp phi, Matrix Psi, ld theta, ll n)
{
    ld beta_term = beta(theta);
    return beta_term * beta_term * Psi.a[n][n].A + beta_term * phi.A * cos(phi.theta - theta);
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

ld getPathLoss(ld distance, ld exponent) {
    ld ref_loss_db = 40.0; 
    ld ref_loss_linear = pow(10, -ref_loss_db / 10.0);
    return ref_loss_linear * pow(distance, -exponent);
}
int main()
{
    ll n = 40; // IRS elements
    ll m = 2;  // AP antennas

    Matrix h_d, h_r, G; 
    h_d.prep(m, 1); h_r.prep(n, 1); G.prep(n, m);

    vector<ld> dist;
    vector<ld> AR;
    for (ld d = 480; d <= 500; d += 1) 
    {
        ld d_ap_irs = 500.0; 
        ld d_ap_user = sqrt(d * d + 2.0 * 2.0);
        ld d_irs_user = sqrt((500.0 - d) * (500.0 - d) + 2.0 * 2.0);
        
        ld h_d_gain = sqrt(getPathLoss(d_ap_user, 3.8));
        ld G_gain = sqrt(getPathLoss(d_ap_irs, 2.2));
        ld h_r_gain = sqrt(getPathLoss(d_irs_user, 2.8));

        ld total_rate = 0.0L;
        ll num_trials = 10; // monte-carlo

        for (ll trial = 0; trial < num_trials; ++trial)
        {
            for (ll i = 0; i < h_d.row; ++i)
                for (ll j = 0; j < h_d.col; ++j)
                    h_d.a[i][j] = {h_d_gain * rayleigh(gen), phase_gacha(gen)};

            for (ll i = 0; i < h_r.row; ++i)
                for (ll j = 0; j < h_r.col; ++j)
                    h_r.a[i][j] = {h_r_gain * rayleigh(gen), phase_gacha(gen)};

            for (ll i = 0; i < G.row; ++i)
                for (ll j = 0; j < G.col; ++j)
                    G.a[i][j] = {G_gain * rayleigh(gen), phase_gacha(gen)};


            Matrix Hh_r = h_r, diagh_r = h_r, Hdiagh_r = h_r;
            diagh_r.diag();
            Hdiagh_r.hermitian();
            Hdiagh_r.diag();
            Hh_r.hermitian();

            Matrix HG = G;
            HG.hermitian();

            Matrix Psi = mult(mult(mult(Hdiagh_r, G), HG), diagh_r);
            Matrix Hat_h_d = mult(mult(Hdiagh_r, G), h_d);

            Matrix v;
            v.prep(n, 1);

            // Ideal IRS: start amplitude = 1; end amplitutde = 1;
            // Practical IRS with Ideal Assumption: start amplitude = 1; end amplitutde = beta(theta);
            for (ll i = 0; i < n; ++i)
            {
                ld init_theta = (percen_gacha(gen) > 0.5) ? M_PI : -M_PI;
                // Ideal IRS - Practical IRS with Ideal Assumption
                v.a[i][0] = {1.0L, init_theta};

                // AO with Proposition 1 - AO with 1D search
                //v.a[i][0] = {beta(init_theta), init_theta};  
            }

            for (ll cnt = 0; cnt < 15; ++cnt) 
            {
                for (ll i = 0; i < n; ++i)
                {
                    comp phi = Hat_h_d.a[i][0];
                    comp tmp = {0, 0};

                    for (ll j = 0; j < n; ++j)
                    {
                        if (j == i) continue;
                        tmp = sumComplex(tmp, multComplex(Psi.a[i][j], v.a[j][0]));
                    }
                    phi = sumComplex(phi, tmp);
                    phi.A *= 2.0L;
                    
                    while (phi.theta > M_PI) phi.theta -= 2.0L * M_PI;
                    while (phi.theta < -M_PI) phi.theta += 2.0L * M_PI;

                    // Upper Bound: Ideal IRS
                    ld max_theta = phi.theta;
                    
                    // 1D Search from -PI to PI
                    // ld max_val = -1e18L;
                    // ld max_theta = v.a[i][0].theta;
 
                    // for (ll step = 0; step <= SEARCH_STEPS; ++step) 
                    // {
                    //     ld test_theta = -M_PI + (2.0L * M_PI * step / SEARCH_STEPS);
                    //     ld current_f = f(phi, Psi, test_theta, i);
                        
                    //     if (current_f > max_val) 
                    //     {
                    //         max_val = current_f;
                    //         max_theta = test_theta;
                    //     }
                    // }

                    // Proposition 1
                    // ld f1 = f(phi, Psi, phi.theta, i);
                    // ll lambda = (phi.theta >= 0) ? 0 : 1;
                    // ld f3 = f(phi, Psi, pow(-1, lambda) * M_PI, i);
                    // ld f2 = f(phi, Psi, (phi.theta + pow(-1, lambda) * M_PI) / 2.0, i);

                    // ld den = 4.0L * (f1 - 2.0L*f2 + f3);
                    // ld max_theta = phi.theta;
                    // if (abs(den) > 1e-9L) {
                    //     max_theta = pow(-1, lambda) * M_PI * (3.0L*f1 - 4.0L*f2 + f3) + phi.theta * (f1 - 4.0L*f2 + 3.0L*f3);
                    //     max_theta /= den;
                    // }

                    while (max_theta > M_PI) max_theta -= 2.0L * M_PI;
                    while (max_theta < -M_PI) max_theta += 2.0L * M_PI;

                    v.a[i][0].theta = max_theta;
                    v.a[i][0].A = 1.0L;
                    // Practical IRS under Ideal Assumption
                    // v.a[i][0].A = beta(max_theta);
                }
            }
            for (ll i = 0; i < n; ++i) {
                    v.a[i][0].A = beta(v.a[i][0].theta); 
                }
            total_rate += calculateRate(v, h_r, G, h_d, noise);
        }

        dist.push_back(d);
        AR.push_back((total_rate / num_trials));
        //cout << "Done: d = " << d << '\n';
    }


    // for (auto i: dist) cout << i << ", ";
    // cout << '\n';
    for (auto i: AR) cout << i << ", ";
    return 0;
}