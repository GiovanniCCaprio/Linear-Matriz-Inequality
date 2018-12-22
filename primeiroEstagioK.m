function [feas,K,linhas,V] = primeiroEstagioK(A,B,x,vertices)

% Inicialização
linhas = 0;   % LMI contagem de linhas
n = size(A{1},1);
Ai= [];
Bi= [];
LMIs = [];

for i = 1:vertices
    Ai = [Ai A{i}];
    Bi = [Bi B{i}];
end 


%Criação das Variáveis
polyA = rolmipvar(Ai,'A(a)',vertices,1);
polyB = rolmipvar(Bi,'B(a)',vertices,1);
polyP = rolmipvar(n,n,'P(a)','symmetric',vertices,1);

Z = rolmipvar(1,n,'Z(a)','full',vertices,1); % TORNA O Z em Z(a)
X = rolmipvar(n,n,'X','full',vertices,0);

% LMIs
LMIs = [LMIs, polyP >= 0];
linhas = linhas + n; % Soma grau de P

T11 = polyA*X + X'*polyA' + polyB*Z + Z'*polyB';
T12 = polyP - X' + x*polyA*X + x*polyB*Z;
T22 = -x*(X + X');
T = [ T11  T12;
      T12' T22];
LMIs = [LMIs, T <= 0];

linhas = linhas + 2*n; % Soma grau de T

% Resolução
optimize(LMIs,[],sdpsettings('verbose',0,'solver','sedumi'));
res = min(checkset(LMIs));

% Número de Variáveis e Linhas
V = size(getvariables(LMIs),2);
linhas = vertices*linhas; % LMIs totais

if res > 0
    feas = 1;
    K = double(Z)*inv(double(X)); %Acha um K(a) com vértices
    else 
        feas = 0;
        K = NaN;
end

end

