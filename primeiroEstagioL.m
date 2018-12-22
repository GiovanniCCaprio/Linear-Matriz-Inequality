function [feas,L,linhas,V] = primeiroEstagioL(A,C,x,v)

% Inicialização
linhas = 0;   % LMI contagem de linhas
n = size(A{1},1);
Ai= [];
Ci= [];
LMIs = [];

for i = 1:v
    Ai = [Ai A{i}];
    Ci = [Ci C{i}];
end 

%Criação das Variáveis
polyA = rolmipvar(Ai,'A(a)',v,1);
polyC = rolmipvar(Ci,'C(a)',v,1);
polyP = rolmipvar(n,n,'P(a)','symmetric',v,1);

Z = rolmipvar(n,1,'Z','full',v,0);
X = rolmipvar(n,n,'X','full',v,0);

% LMIs
LMIs = [LMIs, polyP >= 0];
linhas = linhas + n; % Soma grau de P

T11 = X*polyA + polyA'*X' + Z*polyC + polyC'*Z';
T12 = polyP - X + x*polyA'*X' + x*polyC'*Z';
T22 = -x*(X + X');
T = [ T11  T12;
      T12' T22];
LMIs = [LMIs, T <= 0];

linhas = linhas + 2*n; % Soma grau de T

% Resolução
optimize(LMIs,[],sdpsettings('verbose',0,'solver','sedumi'));
res = min(checkset(LMIs));

% Número de Variáveis
V = size(getvariables(LMIs),2);
linhas = v*linhas;

if res > 0
    feas = 1;
    L = inv(double(X))*double(Z);
    else 
        feas = 0;
        L = NaN;
end

end