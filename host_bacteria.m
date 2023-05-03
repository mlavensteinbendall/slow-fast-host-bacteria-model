% clear workspace 
clear 
clc
close all

% Objective: The goal of this code is to recreate plots from the preprinted
% paper "Slow-fast model of host-bacteria interaction system" by Yvon
% Maday, Darryl Ondoua, Benoit Sarels, which models a host health is
% dependent on the existence of two bacteria.

% This code is modeling the health of the host with relation to two
% colonies of bacteria.


% Variables for Bacteria 1
a_x = 2;
b_x = 1;
c_x = 1.5;
alpha_x = 0.2;
beta_x = 3;

% Variables for Bacteria 2
a_y = 2;
b_y = 1;
c_y = 1.5;
alpha_y = 0.2;
beta_y= 3;


% Variables for relationship between bacteria and host health
gamma_0 = 0.25;
gamma_1 = 0.8;


%Time
eps = 1;            % time scaling T = eps * t (0 < eps << 1)
t_max = 5;          % maximum time
t_interval = [0,7]; % time interval


% Initial conditions (Bacteria 1, Bacteria 2, Host Health)
init_cond = [0.7, 0.2, 0.5]';



% Solve system of ODE -- Slow Model
[t,sol] = ode45(@(t,Y) odefcn(t,Y,a_x,b_x,c_x,alpha_x,beta_x,a_y,b_y,c_y,alpha_y,beta_y,gamma_1,gamma_0,eps,t_max) , t_interval , init_cond);

% Solve system of ODE -- Fast Model
% [t2,sol2] = ode45(@(t2,X) odefcn(t2,X,a_x,b_x,c_x,alpha_x,beta_x,a_y,b_y,c_y,alpha_y,beta_y,gamma_1,gamma_0,eps,t_max) , t_interval2 , init_cond);



% Plot -- Solution of the ODEs for the host health and bacterias over slow time
figure(1)
hold on;
plot(t,sol(:,1),'r');
plot(t,sol(:,2),'b');
plot(t,sol(:,3),'g');
hold off;
title('Values over Slow Time')
ylabel('Values')
xlabel('Slow Time (T)')
legend('Bacteria 1','Bacteria 2','Host')

% Plot -- Solution of the ODEs for the host health and bacterias over fast time
% figure(2)
% hold on;
% plot(t2,sol2(:,1),'r');
% plot(t2,sol2(:,2),'b');
% plot(t2,sol2(:,3),'g');
% hold off;
% title('Values over Fast Time')
% ylabel('Values')
% xlabel('Fast Time (t)')
% legend('Bacteria 1','Bacteria 2','Host')

% Defining Nullclines of the bacteria 
nullcline_x = @(x,z) (-1/alpha_x)*log((1/a_x)*(b_x*exp(-beta_x*x)+c_x + gamma_0 - gamma_1 * z));
nullcline_y = @(x,z) (-1/alpha_y)*log((1/a_y)*(b_y*exp(-beta_y*x)+c_y+ (-gamma_1 * z + gamma_0)));

n=0:.01:3;


figure(2)
plot(n,nullcline_x(n,0.165))
hold on;
plot(nullcline_y(n,0.165),n)
plot(0,0,'o')
title('Nullclines for Bacterias')
xlabel('Bacteria 1')
ylabel('Bacteria 2')
legend('phi(s,n)=0','psi(s,n)=0')



% Y=[ds/dt, dn/dt] -- Slow Model
function dYdt = odefcn(t,Y,a_s,b_s,c_s,alpha_s,beta_s,a_n,b_n,c_n,alpha_n,beta_n,gamma_1,gamma_0,eps,t_max)
 dYdt = [
        Y(1)*(phi(Y,a_s,b_s,c_s,alpha_s,beta_s) + (gamma_1 * Y(3) - gamma_0))/eps;
        Y(2)*(psi(Y,a_n,b_n,c_n,alpha_n,beta_n) +(gamma_1 * Y(3) - gamma_0))/eps;
        (Y(3)*(1-Y(3))*(T(t_max,Y(1),Y(2)) - t))]; 
end

% Y=[ds/dt, dn/dt] -- Fast Model
% function dYdt = odefcn(t,Y,a_s,b_s,c_s,alpha_s,beta_s,a_n,b_n,c_n,alpha_n,beta_n,gamma_1,gamma_0,eps,t_max)
%  dYdt = [
%         Y(1)*(phi(Y,a_s,b_s,c_s,alpha_s,beta_s) + (gamma_1 * Y(3) - gamma_0));
%         Y(2)*(psi(Y,a_n,b_n,c_n,alpha_n,beta_n) +(gamma_1 * Y(3) - gamma_0));
%         eps*(Y(3)*(1-Y(3))*(T(t_max,Y(1),Y(2)) - eps*t))]; 
% end


% Function for Bacteria 1
function phi_sol = phi(Y,a_s,b_s,c_s,alpha_s,beta_s)
    phi_sol = (a_s * exp(-Y(1) * alpha_s) - b_s * exp((-Y(2)) * beta_s) - c_s);
end

% Function for Bacteria 2
function psi_sol = psi(Y,a_n,b_n,c_n,alpha_n,beta_n)
    psi_sol = (a_n * exp(-Y(2) * alpha_n) - b_n * exp((-Y(1)) * beta_n) - c_n);
end

% Function for mature age of the host
function T_sol = T(t_max,x,y)
    T_sol = t_max*(x/(1+x))*(y/(1+y));
end
