function [sysd, K] = c2d_equiv2(sysc, h)

A = sysc.A;
B = sysc.B;
C = sysc.C;
D = sysc.D;

nx = size(A,1);
nu = size(B,2);

Abar = [A B;
        zeros(nu,nx+nu)];

Cbar = [C D];

%% Ad, Bd (igual ZOH)
Ah = expm(A*h);

Md = expm(Abar*h);
Bh = Md(1:nx, nx+1:end);

%% Integral numérica direta (sem truque)
f = @(tau) expm(Abar'*tau) * (Cbar'*Cbar) * expm(Abar*tau);

Qh = integral(f, 0, h, 'ArrayValued', true);

%% Cholesky
ChDh = chol(Qh);

Ch = ChDh(:, 1:nx);
Dh = ChDh(:, nx+1:end);

%% Ganho LQR

Q = Qh(1:nx,1:nx);
S = Qh(1:nx,nx+1:end);
R = Qh(nx+1:end,nx+1:end);

K = -dlqr(Ah, Bh, Q, R, S);

%% Sistema Discreto

sysd = ss(Ah,Bh,Ch,Dh,h);

end