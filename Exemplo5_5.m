clc

% Sistema Contínuo
A = [0 1; -6 1];
B = [0; 1];
C = [1 0; 0 0];
D = [0; 1];
E = [1; 1];

sysc = ss(A,B,C,D);
n_x = size(sysc.A,1);

% Períodos de Amostragem
H = [1.35 2.4];
n_H = size(H,2);

% Sistemas a Tempo Discreto
sysds = cell(n_H, 1);
Ks = cell(n_H, 1);

for i = 1:n_H
    [sysds{i}, Ks{i}] = c2d_equiv2(sysc, H(i));
end

% Matriz de Metzler
p = 1;
q = 0;

Pi = [p 1-q;
      1-p q];

[Ps, mi] = Prop5_6(sysds, Ks, E, Pi);

tf = 10;
[t, x, z, J] = SimHibr(sysc, H, Ks, tf, E, Ps);

J
mi
plot(t, sum(z.^2,2))

