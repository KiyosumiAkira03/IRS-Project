# Intelligent Reflecting Surface: Practical Phase Shift Model and Beamforming Optimization

Nguyen Hoang Minh - 202517030 - Minh.NH2517030@sis.hust.edu.vn

## 1. Paper summarize

An Intelligent Reflecting Surface (IRS) is an array of passive electromagnetic elements capable of dynamically altering the phase shifts of incoming radio waves. By automatically reconfiguring these reflections, an IRS can control the wireless propagation environment to boost signal strength at a receiver. 
Prior academic research overwhelmingly relied on an ideal phase shift model, which assumes that each IRS element can shift a wave's phase freely while maintaining perfect, full reflection (unity amplitude, or $100\%$ reflection efficiency).  
This violates practical hardware constraints and physics. When a passive element shifts a wave's phase, its internal electrical impedance changes. This change in impedance inevitably alters how much energy the element absorbs versus how much it reflects. Consequently, the reflection amplitude is non-uniform and coupled directly with the phase shift.  

To solve this, the authors:
- Propose a **practical phase shift model** that mathematically captures this phase-dependent amplitude variation.  
- Formulate a **joint optimization** problem to properly balance phase alignment and amplitude loss to maximize the system's data transmission rate. 


## 2. System Model

The goal here is to maximize the **Achievable rate** by:
- Optimizing transmit beamforming
- "Min-max"ing IRS reflect beamforming 

The **Achievable rate** can be calculated as:
		$$R = log_2(1 + SNR)$$ 
in which:
$R$ is measured by bps - bits per second.
$SNR$: Signal-plus-Noise Ratio, used to measure the quality of received signal.
In other words, the higher $SNR$ is, the better the output signal

$$SNR = {P_{signal} \over P_{noise}}$$

$P_{signal}$: The power of signal
$P_{noise}$: The power of noise

First, the formula for $P_{signal}$ must be identified

Transmit beamforming: 
- Using multiple antennas to focus signal towards a specific receiver
- Optimize: By changing the position, angle of those multiple antennas (i.o.w: “phase shifts”)

The transmitted vector of the “kth” user equipment (UE):
		$$x_k = b_k * d_k$$
$x_k$: beamformed sample (transmitted vector) for “kth UE”
$b_k$: beamforming vector (precoding weight) for “kth UE”
$d_k$: data symbol of “kth UE” (it is just the data to send)

Total transmitted vector: $w$ = $\sum(x_k)$ = $\sum(b_k * d_k)$

Access Point (AP)
Let:
$h_d$ : baseband equivalent channels from AP to user
$h_r$  : baseband equivalent channels from IRS to user
$G$    :  baseband equivalent channels from AP to IRS
$h_d$ : $(N \times 1 )$ 
$h_r: (M \times 1)$
$G$  : $(N \times M)$ 

$v$: reflection coefficient vector of the IRS
$v_n = [0, 1]$
$arg(v_n)$: The angle of each IRS element = $[-\pi, \pi]$

The **received baseband signal** at the user is thus given by:
$$V_{signal} = (v^H Φ + h^H_d)w*s + z$$ 
in which:
$s$ stands for the data symbol to be tranferred
$h^H_d$ is the direct path from AP to user
$\rightarrow$ $h^H_d \times w$ represents the actual physical wave that survives in the direct path from AP to user 

$Φ = diag(h^H_r)G$
$G$ is the path from AP to the IRS
$h^H_r$ is the path from the IRS to user
Therefore the combined path from AP to user is:
$$h^H_r * G$$
But since $h^H_r$ is a $(M \times 1)$ matrix, $G$ is a $(N \times M)$ matrix, the multiplication is not defined. Therefore we must represent $h^H_r$ in a different way.
Notice that the path of each element of the AP to the IRS is represented as a row in $G$ (total of $N$ rows), and each column represent the path from each element of the IRS to the AP (total of $M$ columns). Therefore we only need to combine each element of $h^H_r$ with the corresponding row of $G$:

$$
\begin{bmatrix}
h_1 & 0 & 0 & ... & 0 \\
0 & h_2 & 0 & ... & 0 \\
... \\
0 & 0 & 0 & ... & h_n 
\end{bmatrix}
$$
$$
\begin{bmatrix}
G_{11} & G_{12} & ... & G1N \\
G_{21} & G_{22} & ... & G2N \\
...  \\
G_{M1} & G_{M2} & ... & GMN 
\end{bmatrix}
$$   
$\rightarrow$ $diag(h^H_r) * G * w$ represents the actual physical wave that survives in the reflection process through the IRS to user.

$z$ denotes the additive white Gaussian noise (AWGN) at the receiver with zero mean and variance $V_{noise} = σ^2$.

Accordingly, **the achievable rate** can be rewritten into:
$$R_{SE} = log_2(1 + {V_{signal}^2 / R \over V_{noise}^2 / R } ) = log_2(1 + {|(v^H Φ + h^H_d)w|^2 \over σ^2})$$
$\rightarrow$ The ultimate goal is to maximize $R_{SE}$

## 3. Practical Phase Shift Model

According to the authors, the total electrical impedance ($Z_n$) opposes to an incoming wave at a given radio frequency ($\omega$) is given by:

$$Z_{n}(C_{n},R_{n})=\frac{j\omega L_{1}(j\omega L_{2}+\frac{1}{j\omega C_{n}}+R_{n})}{j\omega L_{1}+(j\omega L_{2}+\frac{1}{j\omega C_{n}}+R_{n})}$$

in which: 
$L_1$ denote the bottom layer inductance
$L_2$ denote the top layer inductance
$C_n$ denote the effective capacitance
$R_n$ denote the effective resistance


Therefore, the reflection coefficient of the nth element ($v_n$) can be determined by:
$$v_n = \frac{Z_n(C_n, R_n) - Z_0}{Z_n(C_n, R_n) + Z_0}$$

Since $v_n$ is a function of $C_n$ and $R_n$, the reflected electromagnetic waves can be manipulated in a controllable and programmable manner by varying $C_n$’s and $R_n$’s.

However, because both the formula for $Z_n$ and $v_n$ are nested and non-linear. Attempting to use them directly to optimize an array of hundreds of IRS elements simultaneously would require immense computational power, making real-time adjustments impossible for a fast-moving cellular network.
## 4. Optimization (AO Algorithm)

The ultimate goal still remains the same:

$$\max_{\mathbf{w}, \mathbf{v}, \{\theta_{n}\}} |(\mathbf{v}^{H}\mathbf{\Psi}+\mathbf{h}_d^H)\mathbf{w}|^{2}$$$$\text{s.t.} \quad \|\mathbf{w}\|_2^2 \le P_T \quad \text{(7)}$$$$v_n = \beta_n(\theta_n)e^{j\theta_n}, \quad \forall n \quad \text{(8)}$$$$-\pi \le \theta_n \le \pi, \quad \forall n \quad \text{(9)}$$

in which:
$$ \quad \beta_{n}(\theta_{n})=(1-\beta_{min})\left(\frac{sin(\theta_{n}-\phi)+1}{2}\right)^{k}+\beta_{min}$$
is the amplitude of each element, calculated directly with a phase shift $\theta_n$, which act as a "cost" for each phase shift changes.

Because optimizing all variables at once is mathematically non-convex and highly complex, the authors use **Alternating Optimization (AO)**. They hold the Access Point beamforming ($\mathbf{w}$) and $N-1$ IRS elements fixed, isolating just one single IRS element ($n$) at a time, solve for the optimal phase shift for each, rinse and repeat.

When an element is isolated, the problem shrinks down into:$$\max_{\theta_{n}} \beta_{n}^{2}(\theta_{n})\Psi_{n,n} + \beta_{n}(\theta_{n})|\varphi_{n}|\cos(\arg(\varphi_{n})-\theta_{n})$$$$\text{s.t.} \quad -\pi \le \theta_n \le \pi$$

