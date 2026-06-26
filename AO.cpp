#include <bits/stdc++.h>
using namespace std;
#define ll long long
#define ld long double
const ll mod = 1e9 + 7;

random_device rd;
mt19937 gen(rd());
uniform_real_distribution<> phase_gacha(-M_PI, M_PI);
uniform_real_distribution<> percen_gacha(0.0, 1.0); 
uniform_real_distribution<> rayleigh(0.5, 1.5); 

struct comp // A * e^(i * phi) = A(cos(phi) + i * sin(phi))
{
    ld A; // amplitude
    ld phi; // phase
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
                cout << "[" << a[i][j].A << "; " << a[i][j].phi << "]";
            }
            cout << '\n';
        }
        
    }

    // 1 1 1 --> 1 0 0
    //           0 1 0 
    //           0 0 1
    void diag()
    {
        if (row != 1) return;
        vector<vector<comp> > tmp;
        tmp.resize(col, vector<comp>(a[0].size()));

        for (ll i = 0; i < col; ++i)
        {
            for (ll j = 0; j < col; ++j)
            {
                if (i == j) tmp[i][j] = a[0][i];
                else tmp[i][j] = {0, 0};
            }
        }
        
        a.resize(col, vector<comp>(col));
        for (ll i = 0; i < col; ++i)
        {
            for (ll j = 0; j < col; ++j)
            {
                a[i][j] = tmp[i][j];
            }
        }
        row = col;
    }

    void transpose()
    {
        vector<vector<comp>> tmp;
        tmp.resize(col, vector<comp>(row));

        for (ll i = 0; i < row; ++i)
        {
            for (ll j = 0; j < col; ++j)
            {
                tmp[j][i] = a[i][j]; 
            }
        }
        
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
        vector<vector<comp>> tmp;
        tmp.resize(col, vector<comp>(row));

        for (ll i = 0; i < row; ++i)
        {
            for (ll j = 0; j < col; ++j)
            {
                tmp[j][i] = a[i][j]; 
            }
        }
        
        a.resize(col, vector<comp>(row));
        swap(row, col);
        for (ll i = 0; i < row; ++i)
        {
            for (ll j = 0; j < col; ++j)
            {
                a[i][j] = tmp[i][j];
                a[i][j].phi *= -1;
            }
        }
    }
};

comp sumComplex(comp a, comp b)
{
    ld real = a.A * cos(a.phi) + b.A * cos(b.phi);
    ld imag = a.A * sin(a.phi) + b.A * sin(b.phi);
    
    ld A = sqrt(real * real + imag * imag);
    ld phi = atan2(imag, real); 
    
    //chuẩn hóa
    if (phi < 0) phi += 2 * M_PI; 
    
    return {A, phi};
}

comp multComplex(comp a, comp b)
{
    return {a.A * b.A, a.phi + b.phi};
}

Matrix mult(Matrix A, Matrix B)
{
    Matrix tmp2;
    tmp2.prep(A.row, B.col);
	for (ll i = 0; i < A.row; ++i)
	{
		for (ll j = 0; j < A.col; ++j)
		{
			for (ll k = 0; k < B.col; ++k)
			{
				tmp2.a[i][k] = sumComplex(tmp2.a[i][k], multComplex(A.a[i][j], B.a[j][k]));
			}
		}
	}
	return tmp2;
}

int main()
{
	ios_base::sync_with_stdio(0); cin.tie(0); cout.tie(0);
    ll n, m; // m: AP element n: IRS element
    cin >> n >> m;

    Matrix h_d, h_r, G; // h_d: AP-to-user (Mx1) h_r: IRS-to-user (Nx1) G: AP-to-IRS (MxN)
    h_d.prep(m, 1); h_r.prep(n, 1); G.prep(m, n);

    for (ll i = 0; i < h_d.row; ++i)
    {
        for (ll j = 0; j < h_d.col; ++j)
        {
            h_d.a[i][j] = {2.6 * 1e-7 * rayleigh(gen), phase_gacha(gen)};
        }
    }

    for (ll i = 0; i < h_r.row; ++i)
    {
        for (ll j = 0; j < h_r.col; ++j)
        {
            h_r.a[i][j] = {2.6 * 1e-3 * rayleigh(gen), phase_gacha(gen)};
        }
    }

    for (ll i = 0; i < G.row; ++i)
    {
        for (ll j = 0; j < G.col; ++j)
        {
            G.a[i][j] = {2.6 * 1e-6 * rayleigh(gen), phase_gacha(gen)};
        }
    }
    
}
