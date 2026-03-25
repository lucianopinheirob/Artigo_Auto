function [Ks, P] = Prop5_6(sysc, E, Pi, H)

n_H = size(H,2);
n_A = size(sysc.A,2)
Ks = zeros(n_H, n_A);

cvx_begin sdp quiet

variable mi nonnegative
variable Ps(n_A,n_A,n_H) semidefinite

for i = 1:n_H
    [Ad,Bd,Cd,Dd] = ssdata(c2d(sysc, H(i), 'zoh'));
    Q = Cd'*Cd;
    R = Dd'*Dd;
    S = Cd'*Dd;
    Ks(i, :) = dlqr(Ad, Bd, Q, R, S);
    %% Escrever desigualdades 5.59 e 5.60
end

minimize mi
subject to
    for i = 1:Ns
        des5_59{i} <= 0;
        des5_60{i} <= 0;
    end
cvx_end

end