%% Trabalho Final - IA 892
clc, clear all, close all

%% Inicialização

xs = [10^(-5) 10^(-3) 10^(-1) 1];
load('DB_dof.mat');
output.tabela1 = [];
output.tabela2 = [];


%% Lema 1 - Teorema 1  

display('Lema 1 - Teorema 1')

for d=1:size(dimensoes,1)
    ordem    = dimensoes(d,1);
    entradas = dimensoes(d,2);
    saidas   = dimensoes(d,3);
    vertices = dimensoes(d,4);
    placar1  = [0 0]; % Soma Estaveis
    placar_v = [0 0]; % Variáveis
    placar_l = [0 0]; % Linhas de LMI
    for i= 1:totalSistemas
            A = BASE{ordem,entradas,saidas,vertices,i}.A;
            B = BASE{ordem,entradas,saidas,vertices,i}.B;
            C = BASE{ordem,entradas,saidas,vertices,i}.C;
            feas = [0 0];
            feas_1 = 0;
            feas_2 = 0;
        for t = 1:4
            x = xs(t);
            [feas_1,K,L1,V1] = primeiroEstagioK(A,B,x,vertices); % 1 Lema
            placar_v(1) = V1;
            placar_l(1) = L1;
            if feas_1 == 1
                    feas = [1 0];
                    [feas_2,Ac,Bc,L2,V2] = segundoEstagio_Cc(A,B,C,K,x,vertices); % 1 Teorema
                    placar_v(2) = V2;
                    placar_l(2) = L2;
                if feas_2 == 1
                    feas = [1 1];
                    break;
                end
            end
        end
        placar1 = placar1 + feas;
    end
    
    % Printar na tela e armazenar iteração
    
    fprintf('terminei [%d %d %d %d]-[%d %d] [%d %d] [%d %d] \n',ordem,entradas,saidas,vertices...
        ,placar1(1),placar1(2),placar_v(1),placar_v(2),placar_l(1),placar_l(2));
    output.tabela1 = [output.tabela1; [ordem,entradas,saidas,vertices,placar1(1),placar1(2)...
        ,placar_v(1),placar_v(2),placar_l(1),placar_l(2)]];
end


%% Lema 2 - Teorema 2

display('Lema 2 - Teorema 2')

for d=1:size(dimensoes,1)
    ordem = dimensoes(d,1);
    entradas = dimensoes(d,2);
    saidas = dimensoes(d,3);
    vertices = dimensoes(d,4);
    placar2 = [0 0]; % Soma Estaveis
    placar_v = [0 0]; % Variáveis
    placar_l = [0 0]; % Linhas de LMI
    for i= 1:totalSistemas
            A = BASE{ordem,entradas,saidas,vertices,i}.A;
            B = BASE{ordem,entradas,saidas,vertices,i}.B;
            C = BASE{ordem,entradas,saidas,vertices,i}.C;
            feas = [0 0];
            feas_1 = 0;
            feas_2 = 0;
        for t = 1:4
            x = xs(t);
            [feas_1,L,L1,V1] = primeiroEstagioL(A,C,x,vertices);  % 2 Lema
            placar_v(1) = V1;
            placar_l(1) = L1;
            if feas_1 == 1
                    feas = [1 0];
                    [feas_2,Ac,Cc,L2,V2] = segundoEstagio_Bc(A,B,C,L,x,vertices); % 2 Teorema
                    placar_v(2) = V2;
                    placar_l(2) = L2;
                if feas_2 == 1
                    feas = [1 1];
                    break;
                end
            end
        end
        placar2 = placar2 + feas;
    end
    
    
    % Printar na tela e armazenar iteração
    
    fprintf('terminei [%d %d %d %d]-[%d %d] [%d %d] [%d %d] \n',ordem,entradas,saidas,vertices...
        ,placar2(1),placar2(2),placar_v(1),placar_v(2),placar_l(1),placar_l(2));
    output.tabela2 = [output.tabela2; [ordem,entradas,saidas,vertices,placar2(1),placar2(2)...
        ,placar_v(1),placar_v(2),placar_l(1),placar_l(2)]];
end


