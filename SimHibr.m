function [t_all, x_all, z_all, J] = SimHibr(sysc, H, Ks, tf, x0, Ps)

    A = sysc.A;
    B = sysc.B;
    C = sysc.C;
    D = sysc.D;
    
    nx = size(A,1);
    nz = size(C,1);
    nH = numel(H);
    
    ta = 0;
    idx = 2; % Escolha é com base no prob de otimização? O subsist que dá a "melhor" P?
    H_atual = H(idx);
    K_atual = Ks{idx};
    
    t_all = [];
    x_all = [];
    z_all = [];
    
    x0_atual = x0;

    J = 0;
    
    while ta < tf
        td = min(ta + H_atual, tf);
        
        uk = K_atual * x0_atual;
    
        [t, x] = ode45(@(t,x) (A*x + B*uk), [ta td], x0_atual);
        
        z = (C*x.' + D*uk).';
    
        if isempty(t_all)
            t_all = t;
            x_all = x;
            z_all = z;
        else
            t_all = [t_all; t(2:end)];
            x_all = [x_all; x(2:end,:)];
            z_all = [z_all; z(2:end,:)];
        end
    
        xf = x(end,:)';
    
        V = zeros(1, nH);
        for i = 1:nH
            V(i) = xf' * Ps{i} * xf;
        end
    
        [~, idx] = min(V);
        H_atual = H(idx);
        K_atual = Ks{idx};
        idx

        J = J + trapz(t, sum(z.^2,2));
    
        x0_atual = xf;
        ta = td;
    end

    %J_all = trapz(t_all, sum(z_all.^2,2));

end