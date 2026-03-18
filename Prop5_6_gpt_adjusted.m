function [Ks, P] = Prop5_6_gpt_adjusted(sysc, E, Pi, H)

n_H = size(H,2);
n_A = size(sysc.A,1);

Ks = zeros(n_H, n_A);

% calcular os ganhos uma vez
for i = 1:n_H
    [Ad,Bd,Cd,Dd] = ssdata(c2d(sysc, H(i), 'zoh'));
    
    Q = Cd'*Cd;
    R = Dd'*Dd;
    S = Cd'*Dd;
    
    K = -dlqr(Ad, Bd, Q, R, S);
    Ks(i,:) = K;
end

best_cost = inf;
best_P = [];

for k = 1:n_H

    cvx_begin sdp quiet

    variable Ps(n_A,n_A,n_H) symmetric

    for i = 1:n_H
        
        [Ad,Bd,Cd,Dd] = ssdata(c2d(sysc, H(i), 'zoh'));
        K = Ks(i,:);

        Acl = Ad + Bd*K;
        Ccl = Cd + Dd*K;

        % construção de Ppi
        Ppi = 0;
        for j = 1:n_H
            Ppi = Ppi + Pi(j,i)*Ps(:,:,j);
        end
        
        des5_59{i} = Acl'*Ppi*Acl - Ps(:,:,i) + Ccl'*Ccl;
    end

    % custo que será minimizado neste SDP
    cost = trace(E'*Ps(:,:,k)*E);

    minimize cost

    subject to
        for i = 1:n_H
            des5_59{i} <= 0;
            Ps(:,:,i) >= 0;
        end

    cvx_end

    if cost < best_cost
        best_cost = cost;
        best_P = Ps;
    end

end

P = best_P;
best_cost
end
