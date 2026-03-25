function [z, J] = SimHibr(sysc, H, Ks, tf, x0, Ps)

Ac = sysc.A;
Bc = sysc.B;
Cc = sysc.C;
Dc = sysc.D;

% Dimensões
nx = size(Ac,1);
nz = size(Cc,1);
nu = size(Bc,2);
nH = size(H,2);

x(:,1) = x0;
H_atual = H(1);
K_atual = Ks{1};

% Controle de tempo
ta = 0;
td = H_atual;

while td < tf
    tspan = [ta td];
    [t, x] = ode45(@(t, x) sis(x, Ac, Bc, K_atual), tspan, x0);

    % Estado final da janela
    xf = x(end, :)';

    % Checar qual será o próximo subsistema ativo
    V = zeros(1, nH);
    
    for i = 1:nH
        V(i) = xf' * Ps{i} * xf;
    end
    
    [~, idx] = min(V);

    H_atual = H(idx);
    K_atual = Ks{idx};
    
    % Avança o tempo
    ta = td;
    td = td + H_atual;
    x0 = xf;

end

%J = sum(sum(z.^2, 1));

% --- Dinâmica interna ---
function dxdt = sis(x, A, B, K)
    dxdt = (A+B*K)*x;
end

end