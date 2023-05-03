% clear workspace 
clear 
clc
close all

% Objective: The goal of this code is to recreate plots from the preprinted
% paper "Slow-fast model of host-bacteria interaction system" by Yvon
% Maday, Darryl Ondoua, Benoit Sarels, which models a host health is
% dependent on the existence of two bacteria.

% This code is modeling the two colonies of bacteria with the assumption
% the host is absent and provided all the nutrients.


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

% Note: * Currently set to case provided on page 8  
%       * We want a > b + c for non-dead population 


% Initial conditions (Bacteria 1, Bacteria 2)
init_cond = [0.5,0.25]';

% In relation to Figure 6 (page 10), the intial conditions for the bacteria
% did not provided the same plots. Hence, here are the conditions that gave
% the same plots:
% When init_cond = [0.5,0.21995]', we see the extiction graph
% When init_cond = [0.5,0.219954039]', we see an even unstable saddle point
% When init_cond = [0.5,0.21999]', we see a stable coexistence



% Time interval
t_interval = [0,60];



% Solve system of ODE
[t,sol] = ode45(@(t,Y) odefcn(Y,a_x,b_x,c_x,alpha_x,beta_x,a_y,b_y,c_y,alpha_y,beta_y) , t_interval , init_cond);



% Plot -- Solution of the ODEs for the bacterias over time
figure(1)
hold on;
plot(t,sol(:,1),'r');
plot(t,sol(:,2),'b');
hold off;
ylabel('Values')
xlabel('Time')
legend('Bacteria 1','Bacteria 2')

% Plot -- Solution of Bacteria 1 vs. Solution of Bacteria 2 (was not in
% paper)
figure(2)
plot(sol(:,1),sol(:,2),'r');
hold on;
plot(init_cond(1),init_cond(2),'bo');
% plot(s_null(:,1))
hold off;
title('Bacteria 1 vs. Bacteria 2')
xlabel('Bacteria 1')
ylabel('Bacteria 2')


% Defining Nullclines 
nullcline_x = @(y) (-1/alpha_x)*log((1/a_x)*(b_x*exp(-beta_x*y)+c_x));
nullcline_y = @(x) (-1/alpha_y)*log((1/a_y)*(b_y*exp(-beta_y*x)+c_y));

n = linspace(0,3);

P = InterX([n; nullcline_x(n)],[nullcline_y(n);n]) % Prints the equalibrium points

% Plotting the Nullclines
figure(3)
hold on;
plot(n,nullcline_x(n))
plot(nullcline_y(n),n)
plot(0,0,'o')   % orgin
hold off;
title('Nullclines for Bacterias')
xlabel('Bacteria 1')
ylabel('Bacteria 2')
legend('phi(s,n)=0','psi(s,n)=0')



% Y=[dx/dt, dy/dt] -- system given on page 4 and 5
function dYdt = odefcn(Y,a_s,b_s,c_s,alpha_s,beta_s,a_n,b_n,c_n,alpha_n,beta_n)
 dYdt = [
        Y(1)*((a_s * exp(- Y(1) * alpha_s)) - (b_s * exp(- Y(2) * beta_s)) - c_s); %ds/dt
        Y(2)*((a_n * exp(- Y(2) * alpha_n)) - (b_n * exp(- Y(1) * beta_n)) - c_n)]; %dn/dt
end
