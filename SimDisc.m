function [z, J] = SimDisc(sysd, K, N, x0)

Ah = sysd.A;
Bh = sysd.B;
Ch = sysd.C;
Dh = sysd.D;

% Dimensões
nx = size(Ah,1);
nz = size(Ch,1);
nu = size(Bh,2);

x = zeros(nx, N+1);
z = zeros(nz, N);
u = zeros(nu, N);

x(:,1) = x0;

% Simu
for k = 1:N
    u(:,k) = K * x(:,k);              
    z(:,k) = Ch*x(:,k) + Dh*u(:,k);
    x(:,k+1) = Ah*x(:,k) + Bh*u(:,k);
end

% Custo
J = sum(sum(z.^2, 1));

end