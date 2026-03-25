function [Ps, mi] = Prop5_6(sysds, Ks, E, Pi)

n_H = size(sysds,1);
n_x = size(Ks,1);

melhor_mi = inf;
melhor_P = [];

for k = 1:n_H

    cvx_begin sdp quiet

    variable Ps(n_x,n_x,n_H) semidefinite
    variable mi nonnegative

    for i = 1:n_H
        Ah = sysds{i}.A;
        Bh = sysds{i}.B;
        Ch = sysds{i}.C;
        Dh = sysds{i}.D;
        K = Ks{i};

        Acl = Ah + Bh*K;
        Ccl = Ch + Dh*K;

        % construção de Ppi
        Ppi = 0;
        for j = 1:n_H
            Ppi = Ppi + Pi(j,i)*Ps(:,:,j);
        end
        
        des5_59{i} = Acl'*Ppi*Acl - Ps(:,:,i) + Ccl'*Ccl;
    end

    des5_60 = mi - trace(E'*Ps(:,:,k)*E);
    
    minimize mi
    
    subject to
        des5_60 >= 0;
        for i = 1:n_H 
            des5_59{i} <= 0;
            Ps(:,:,i) >= 0;
        end

    cvx_end

    if mi < melhor_mi
        melhor_mi = mi;
        melhor_P = Ps;
        %k
    end

end

Ps = cell(n_x, 1);
for i = 1:n_x
    Ps{i} = melhor_P(:,:,i);
end

end
