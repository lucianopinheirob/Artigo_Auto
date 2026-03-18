clc
close all

% Parâmetros do Sistema
A = [0 1; -6 1];
B = [0; 1];
C = [1 0; 0 0];
D = [0; 1];
E = [1; 1];

sysc = ss(A,B,C,D);

% Parâmetros Gerais
H = [1.35 2.4];
p = 1;
q = 0;

Pi = [p 1-q;
      1-p q];

[Ks, P] = Prop5_6_gpt(sysc, E, Pi, H)
%[Ks, P] = Prop5_6_gpt_adjusted(sysc, E, Pi, H)


