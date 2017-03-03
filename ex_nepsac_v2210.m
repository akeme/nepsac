%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EPSAC FOR MOBILE ROBOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
format long g

Ts = 0.1;             % sampling time
Tsim = 85;
h = 0.1;           % integrating time


x0 = [0; 0; 0];  % x,y,theta
X = x0;


xa = [0; 0; 0];  % x,y,theta
X1 = xa;

% Runge kutta

% Y(n+1) = Y(n) + (1/6)(k1 +2 k2 +2 k3 +k4);
% k1 = h f(X(n),Y(n))
% k2 = h f(X(n)+h/2, Y(n)+k1/2)
% k3 = h f(X(n)+h/2, Y(n) +k2/2)
% k4 = h f(X(n)+h, Y(n)+k3)

% The next loop difines the reference

for k = 1:round(Tsim/h)
    v(k) = 0.3;
    if k<= 419;
        w(k) = 0.15;
    else
        w(k) = -0.15;
    end
    
    k1 = h*[v(k)*cos(x0(3)); v(k)*sin(x0(3));w(k)];
    k2 = h*[v(k)*cos(x0(3)+k1(3)/2); v(k)*sin(x0(3)+k1(3)/2);w(k)];
    k3 = h*[v(k)*cos(x0(3)+k2(3)/2); v(k)*sin(x0(3)+k2(3)/2);w(k)];
    k4 = h*[v(k)*cos(x0(3)+k3(3)); v(k)*sin(x0(3)+k3(3));w(k)];
    
    X = [X, X(:,k) + (1/6)*(k1+2*k2+2*k3+k4)];
    x0 = X(:,k+1);
    
%     X1 = [X1, X1(:,k)+ h*([v(k)*cos(xa(3));v(k)*sin(xa(3));w(k)])];
%     xa = X1(:,k+1);
    
end

Xr = X;
Ur = [v;w];
plot(Xr(1,:),Xr(2,:),'b--')
hold on
% plot(X1(1,:),X1(2,:),':')
% legend('runge_kutta','euler')
grid

vr = v;
wr = w;

tic

%%
% inicialização do controlador

N = 5;
Q = [1 0 0;
    0 1 0;
    0 0 .1];   %no segundo exemplo foi usado 0.1
R = 0.1*eye(2);
%restrições

vmin = -0.4;
vmax = 0.4;
wmin = -0.4;
wmax = 0.4;
%  vmin = -10;
%  vmax = 10;
%  wmin = -10;
%  wmax = 10;

%lb e ub utilizads nas ferramentas de otmização

% lb = [];
% ub = [];


for k = 1:N
    Qb(3*k-2:3*k,3*k-2:3*k) =Q;
    Rb(2*k-1:2*k,2*k-1:2*k)=R;
 %  lb = [lb;[vmin;wmin]];
 %  ub = [ub;[vmax;wmax]]
end

%saída Qb , Rb 
% lb e ub -> só se for utilizado ferramentas de otimização


%%
% matrizes de ponderação do Nepsac
Qn = Qb;
%Rn = [R(1,1)*eye(N),zeros(N);zeros(N),R(2,2)*eye(N)];



x0 = [0;-1; pi/2];  % x,y,theta , posição inicial do robô
X = x0;
v = [];
w = [];

Xm = x0;

%%

% determinação do controlador
%entradas : Ur, Qb, Rb, N 

for k = 1:round(Tsim/h)-N;
    % Algoritmo do nepsac
    if k ==1
        Ubase = Ur(:,k:N+k-1);    % [v;w]
    else
        Ubase = [Ubase(:,2:end),Ur(:,k+N-1)];    % [v;w]
    end
    Uopt = 1;
    contador = 0;
    while sum(abs(Uopt)) > 10e-4;
        contador = contador +1;
        % cálculo da saída base
        pert = (x0-Xm);
        %pert = [0;0;0];
        
        for j = 1:N;
            if j ==1
                Xbase(1,j) = x0(1) + Ts*Ubase(1,j)*cos(x0(3));%+Ts*pert(1);
                Xbase(2,j) = x0(2) + Ts*Ubase(1,j)*sin(x0(3));%+Ts*pert(2);
                Xbase(3,j) = x0(3) + Ts*Ubase(2,j);%+Ts*pert(3);
            else
                Xbase(1,j) = Xbase(1,j-1) + Ts*Ubase(1,j)*cos(Xbase(3,j-1));%+pert(1);
                Xbase(2,j) = Xbase(2,j-1) + Ts*Ubase(1,j)*sin(Xbase(3,j-1));%+pert(2);
                Xbase(3,j) = Xbase(3,j-1) + Ts*Ubase(2,j);%+pert(3);
            end
        end
        
        dUv = 0.00001;
        dUw = 0.00001;
        
        for j = 1:N;
            if j ==1
                %  V
                Xv(1,j) = x0(1) + Ts*(Ubase(1,j)+dUv)*cos(x0(3));
                Xv(2,j) = x0(2) + Ts*(Ubase(1,j)+dUv)*sin(x0(3));
                Xv(3,j) = x0(3) + Ts*Ubase(2,j);
                % W
                Xw(1,j) = x0(1) + Ts*(Ubase(1,j))*cos(x0(3));
                Xw(2,j) = x0(2) + Ts*(Ubase(1,j))*sin(x0(3));
                Xw(3,j) = x0(3) + Ts*(Ubase(2,j)+dUw);
                
            else
                %  V
                Xv(1,j) = Xv(1,j-1) + Ts*Ubase(1,j)*cos(Xv(3,j-1));
                Xv(2,j) = Xv(2,j-1) + Ts*Ubase(1,j)*sin(Xv(3,j-1));
                Xv(3,j) = Xv(3,j-1) + Ts*Ubase(2,j);
                % W
                Xw(1,j) = Xw(1,j-1) + Ts*Ubase(1,j)*cos(Xw(3,j-1));
                Xw(2,j) = Xw(2,j-1) + Ts*Ubase(1,j)*sin(Xw(3,j-1));
                Xw(3,j) = Xw(3,j-1) + Ts*Ubase(2,j);
                
            end
        end
        av = (Xv-Xbase)/dUv;
        aw = (Xw-Xbase)/dUw;
        
        Gv = zeros(3*N,N);
        Gw = zeros(3*N,N);
        
        pcv = [];
        pcw = [];
        Gnep = zeros(3*N,2*N);
        nep = [];
        xbnep = [];
        for j = 1:N
            nep = [nep;[av(:,j),aw(:,j)]];
            Gnep((3*(N-j)+1):end,2*(N-j)+1:2*(N-j+1)) = nep;
            pcv = [pcv;av(:,j)];
            pcw = [pcw;aw(:,j)];
            Gv(3*(N-j)+1:end,N+1-j) = pcv;
            Gw(3*(N-j)+1:end,N+1-j) = pcw;
            xbnep = [xbnep;Ubase(:,j)];
        end
        G = [Gv,Gw];
        
        
        
        % Montar a Matriz Dinâmica de resposta ao impulso
        
        teste = 1;
        
        % Fim do nepsac
        
        
        
        %         A = [];
        %         xb = [];
        %         aux = eye(3);
        %         for j = 0:N-1
        %             Aj{j+1} = [1 0 -Ts*vr(k+j)*sin(Xr(3,k+j));
        %                 0 1  Ts*vr(k+j)*cos(Xr(3,k+j));
        %                 0 0  1];
        %             Bj{j+1} = [cos(Xr(3,k+j))*Ts  0;
        %                 sin(Xr(3,k+j))*Ts  0;
        %                 0                Ts];
        %             aux = Aj{j+1}*aux;
        %             A = [A;aux];
        %             xb = [xb;Ur(:,k+j)];
        %         end
        %
        %         B = [];
        %
        %         for j = 1:N
        %             aux2  = eye(3);
        %             Baux = [];
        %             for i = 1:j
        %                 if i<j
        %                     Baux = [aux2*Bj{j-i+1},Baux];
        %                     aux2 = aux2*Aj{j-i+1};
        %                 else
        %                     Baux = [aux2*Bj{j-i+1},Baux,zeros(3,(N-j)*2)];
        %                 end
        %             end
        %             B = [B;Baux];
        %
        %         end
        %
        %         H = 2*(B'*Qb*B+Rb);
        %         f = 2*B'*Qb*A*(X(:,k)-Xr(:,k));
        %Uopt = quadprog(H,f,[],[],[],[],[lb-xb],[ub-xb]);
        
        Hnep = 2*(Gnep'*Qb*Gnep+Rb);
        XbXr = [];
        for j = 1:N
            XbXr = [XbXr;(Xbase(:,j)-Xr(:,k+j))];
        end
        fnep = 2*Gnep'*Qb*XbXr;
 %       options = optimset('Display','off','Largescale','off');
 %       Uopt = quadprog(Hnep,fnep,[],[],[],[],[lb-xbnep],[ub-xbnep],[],options);
       Uopt = -inv(Hnep)*fnep;
        U(:,k) = Ur(:,k)+Uopt(1:2,1);
        u01 = [];
        for j = 1:N
            u01 = [u01, Uopt(2*(j-1)+1:2*j)];
        end
        Ubase= Ubase+u01;
        Uopt =0; %comentando vira o NEPSAC - aviso do programa do Labview
    end
    custo(k) = contador;
    U(:,k) = Ubase(:,1);%+Uopt(1:2,1);
    %U(:,k) = Ur(:,k)+Uopt(1:2,1);
    v_r(k) = Ur(1,k);
    w_r(k) = Ur(2,k);
    v(k) = U(1,k);
    w(k) = U(2,k);
    
     %restrições
    
    if v(k) >=0.4;
        v(k)=0.4;
    end
    
    if v(k)< -0.4;
        v(k) = -0.4;
    end
    
    if w(k)>= 0.4;
        w(k) = 0.4;
    end
    
    if w(k)<-0.4;
        w(k)=-0.4;
    end
    
    %saídas : v(k) e w(k)-> entrada p bloco das velocidades de cada roda
    
    %  k1 = h*[v(k)*cos(x0(3)); v(k)*sin(x0(3));w(k)];
    %  k2 = h*[v(k)*cos(x0(3)+k1(3)/2); v(k)*sin(x0(3)+k1(3)/2);w(k)];
    %  k3 = h*[v(k)*cos(x0(3)+k2(3)/2); v(k)*sin(x0(3)+k2(3)/2);w(k)];
    %  k4 = h*[v(k)*cos(x0(3)+k3(3)); v(k)*sin(x0(3)+k3(3));w(k)];
    
    % X = [X, X(:,k) + (1/6)*(k1+2*k2+2*k3+k4)];
    
    %cálculo da saída do modelo
    % no mathscript é depois do bloco das velocidades do robo
    
    
    Xm(1,1) = x0(1) + Ts*v(k)*cos(x0(3));
    Xm(2,1) = x0(2) + Ts*v(k)*sin(x0(3));
    Xm(3,1) = x0(3) + Ts*w(k);
    X = [X,Xm];
    x0 =Xm;
    %x0 = X(:,k+1);
     
    %no mathscipt calculado pela posição real e a de referência
    
    Erro(:,k) = Xr(:,k+1)- X(:,k+1);
    E1_q(k) = Erro(1,k)*Erro(1,k);
    E2_q(k) = Erro(2,k)*Erro(2,k);
    E3_q(k) = Erro(3,k)*Erro(3,k);
     
    E1_t(k) = k*abs(Erro(1,k));
    E2_t(k) = k*abs(Erro(2,k));
    E3_t(k) = k*abs(Erro(3,k));
    
    E1_tq(k) = k*E1_q(k);
    E2_tq(k) = k*E2_q(k);
    E3_tq(k) = k*E3_q(k);
    %k
end

%%
t_nepsac = toc

plot(X(1,:),X(2,:),'r')
legend('referência','trajetória do robô')
title('EPSAC')
xlabel('x[m]')
ylabel('y[m]')

figure
ne = length(Erro);
tempo = Ts*(0:ne-1);
plot(tempo,v,'r')
hold on
plot(tempo,v_r,'b--')
plot (tempo,w,'r')
plot(tempo,w_r,'b--')
legend('v','w','vr','wr')
title('Velocidades EPSAC')
xlabel('t[s]')
ylabel('v[m/s] w[rad/s]')

figure

plot(tempo,Erro(1,:),':',tempo,Erro(2,:),tempo,Erro(3,:),'-.');
legend('x','y','\theta')
title('Erro EPSAC')
xlabel('t[s]')
ylabel('x[m] y[m] \theta[rad]')


 E1 = Erro(1,:);
 E2 = Erro(2,:);
 E3 = Erro(3,:);
ERRO_TOTAL = E1*E1'*Q(1,1)+E2*E2'*Q(2,2)+E3*E3'*Q(3,3);
grid


%indices de desempenho

% Integral do erro absoluto
IAE1 = trapz(abs(E1))
IAE2 = trapz(abs(E2))
IAE3 = trapz(abs(E3))

% Integral do erro ao quadrado

ISE1 = trapz(E1_q)
ISE2 = trapz(E2_q)
ISE3 = trapz(E3_q)


ITAE1 = trapz(E1_t)
ITAE2 = trapz(E2_t)
ITAE3 = trapz(E3_t)


ITSE1 = trapz(E1_tq)
ITSE2 = trapz(E2_tq)
ITSE3 = trapz(E3_tq)

% desvio1 = std(E1);
% desvio2 = std(E2);
% desvio3 = std(E3);

% erro_max1 = max(abs(E1));
% erro_max2 = max(abs(E2));
% erro_max3 = max(abs(E3));

commandwindow
