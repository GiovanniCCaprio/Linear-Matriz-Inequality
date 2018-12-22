function [feas_2,Ac,Bc,linhas,V] = segundoEstagio_Cc(A,B,C,K,x,vertices) 

% Inicialização
linhas = 0;   % LMI contagem de linhas
n = size(A{1},1);
p = size(C{1},1);
Ai = [];
Bi = [];
Ci = [];
LMIs = [];

for i = 1:vertices
    Ai = [Ai A{i}];
    Bi = [Bi B{i}];
    Ci = [Ci C{i}];
end 

%Criação das Variáveis
polyA = rolmipvar(Ai,'A(a)',vertices,1);
polyB = rolmipvar(Bi,'B(a)',vertices,1);
polyC = rolmipvar(Ci,'C(a)',vertices,1);
polyP = rolmipvar(2*n,2*n,'P(a)','symmetric',vertices,1);

G = rolmipvar(n,p,'G','full',vertices,0);
V = rolmipvar(n,n,'V','full',vertices,0);
H = rolmipvar(n,n,'H','full',vertices,0);
Q = rolmipvar(n,n,'Q','full',vertices,0);
Y = rolmipvar(n,n,'Y','full',vertices,0);

% LMIs
LMIs = [LMIs, polyP >= 0];
linhas = linhas + 2*n; % Soma grau de P

J = [  Q   Q ;
     (Y+V) Y ];

T11 = Q*(polyA + polyB*K);
T12 = Q* polyA;
T21 = Y*(polyA + polyB*K) + G*polyC + H;
T22 = Y*polyA+G*polyC;
T   = [ T11 T12;
        T21 T22];

F11 = T + T';
F12 = polyP - J' + x*T;
F22 = -x*(J + J');
F = [ F11  F12;
      F12' F22];
LMIs = [LMIs, F <= 0];
linhas = linhas + 4*n; % Soma grau de F

% Resolução
optimize(LMIs,[],sdpsettings('verbose',0,'solver','sedumi'));
res = min(checkset(LMIs));

% Número de Variáveis e Linhas
V = size(getvariables(LMIs),2);
linhas = vertices*linhas; % LMIs totais

if res > 0
    feas_2 = 1;
    Ac = inv(double(V))*double(H);
    Bc = inv(double(V))*double(G);
    else 
        feas_2 = 0;
        Ac = NaN;
        Bc = NaN;
end


end

