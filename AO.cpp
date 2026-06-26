#include <bits/stdc++.h>
using namespace std;
#define ll long long
#define ld long double

const ll mod = 1e9 + 7;
const ld BETA_MIN = 0.2; // 0.2
const ld K = 1.6; // 1.6
const ld noise = 1e-12;

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

    // 1 1 1 --> 1 0 0
    //           0 1 0 
    //           0 0 1
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
    comp total_y = {0, 0}; 
    total_y = h_d.a[0][0]; 

    for (ll n = 0; n < v.row; ++n) 
    {
        ld theta_n = v.a[n][0].theta;
        ld beta_n = beta(theta_n); 
        comp reflection = {beta_n, theta_n};
        comp channel_gain = multComplex(h_r.a[n][0], G.a[n][0]); 
        total_y = sumComplex(total_y, multComplex(reflection, channel_gain));
    }
    ld signal_power = total_y.A * total_y.A; 
    ld sinr = signal_power / noise_power;

    return log2(1.0 + sinr);
}
int main()
{
	//ios_base::sync_with_stdio(0); cin.tie(0); cout.tie(0);
    ll n, m; // m: AP element n: IRS element
    cin >> n >> m;

    Matrix h_d, h_r, G; // h_d: AP-to-user (Mx1) h_r: IRS-to-user (Nx1) G: AP-to-IRS (MxN)
    h_d.prep(m, 1); h_r.prep(n, 1); G.prep(n, m);

    
    for (ll i = 0; i < h_d.row; ++i)
    {
        for (ll j = 0; j < h_d.col; ++j)
        {
            h_d.a[i][j] = {2.6 * 1e-7L * rayleigh(gen), phase_gacha(gen)};
        }
    }

    for (ll i = 0; i < h_r.row; ++i)
    {
        for (ll j = 0; j < h_r.col; ++j)
        {
            h_r.a[i][j] = {2.6 * 1e-3L * rayleigh(gen), phase_gacha(gen)};
        }
    }

    for (ll i = 0; i < G.row; ++i)
    {
        for (ll j = 0; j < G.col; ++j)
        {
            G.a[i][j] = {2.6 * 1e-6L * rayleigh(gen), phase_gacha(gen)};
        }
    }

    Matrix Hh_r = h_r, diagh_r = h_r, Hdiagh_r = h_r;
    diagh_r.diag();

    Hdiagh_r.hermitian();
    Hdiagh_r.diag();

    Hh_r.hermitian();

    Matrix HG = G;
    HG.hermitian();

    Matrix Psi = mult(mult(mult(Hdiagh_r, G), HG), diagh_r);
    Matrix Hat_h_d = mult(mult(Hdiagh_r, G), h_d);
    //Hat_h_d.show();
    //Psi.show();

    Matrix v;
    v.prep(n, 1);
    for (ll i = 0; i < n; ++i)
    {
        // standard IRS setup, initial amplitude beta = 1.0, random phase
        v.a[i][0] = {1.0L, phase_gacha(gen)}; 
    }
    for (ll cnt = 0; cnt < 10; ++cnt)
    {
        for (ll i = 0; i < n; ++i)
        {
            // Calculate phi_n
            comp phi = Hat_h_d.a[i][0];
            phi.A *= 2;
            comp tmp = {0, 0};

            for (ll j = 0; j < n; ++j)
            {
                if (j == i) continue;
                tmp = sumComplex(tmp, multComplex(Psi.a[i][j], v.a[j][0]));
            }
            phi = sumComplex(phi, tmp);
            if (phi.theta > M_PI) phi.theta -= 2 * M_PI;
            // 3 points: arg(phi_n) / arg(phi_n) + (-1)^lambda*pi :2 / (-1)^lambda*pi
            // lambda = 0 if arg(phi_n) >= 0, = 1 otherwise
            // f(comp phi, Matrix Psi, ld theta, ll n)

            ld f1 = f(phi, Psi, phi.theta, i);
            ll lambda = (phi.theta >= 0) ? 0 : 1;
            ld f3 = f(phi, Psi, pow(-1, lambda) * M_PI, i);
            ld f2 = f(phi, Psi, (phi.theta + pow(-1, lambda) * M_PI) / 2.0, i);

            ld max_theta = pow(-1, lambda) * M_PI * (3*f1 - 4*f2 + f3) + phi.theta * (f1 - 4*f2 + 3*f3);
            max_theta /= 4*(f1 - 2*f2 + f3);

            //cout << max_theta << ' ' << phi.theta << '\n';
            v.a[i][0].theta = max_theta;
        }
        cout << calculateRate(v, h_r, G, h_d, noise) << ' ';
    }
    
}
