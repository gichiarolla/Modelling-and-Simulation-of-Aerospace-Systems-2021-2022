% Modelling and Simulation of Aerospace Systems (2021/2022)
% Assignment # 2
% Author Giovanni Chiarolla matr. 964831
%% Ex1
clear;clc;close all

T0 = 0.1; % Nm
J2 = 0.1; % kgm
J1 = 0.2; % kgm
tspan = [0, 10];
% k and b assumptions
k = 3;
b = 0.8;
x0 = [0, 0, 0, 0];
options = odeset;
[t,x] = ode45(@(t,x) state_ex1(t,x,J1,J2,k,b,T0),tspan,x0,options);
ddx = acc_ex1(x,J1,J2,k,b,T0);
%plots
figure(1)
plot(t,ddx(:,1),'LineWidth',1.2)
hold on 
plot(t,ddx(:,2),'LineWidth',1.2)
grid on
legend('$\ddot{\theta}_{1}$','$\ddot{\theta}_{2}$','Interpreter','latex','Location','best',...
    'FontSize',13)
xlabel('time [s]','Interpreter','latex','FontSize',13)
ylabel('$ \ddot{\theta}_{i}$ [$\frac{rad}{s^{2}}$]','Interpreter','latex','FontSize',13)
axis ([0 10 -0.4 1])
figure(2)
plot(t,ddx(:,1),'LineWidth',1.2)
hold on 
plot(t,ddx(:,2),'LineWidth',1.2)
grid on

file = fopen('sample.txt','r');
sample = fscanf(file,'%f',[3,inf]);
fclose(file);
sample = sample';
hold on
 
plot(sample(:,1), sample(:,2),'--','LineWidth',1.3)
hold on 
plot(sample(:,1), sample(:,3),'--','LineWidth',1.3)
legend('$\ddot{\theta}_{1}$','$\ddot{\theta}_{2}$','$\ddot{\theta}_{1} (Sample)$','$\ddot{\theta}_{2} (Sample)$','Interpreter','latex','Location','best',...
    'FontSize',13)
xlabel('time [s]','Interpreter','latex','FontSize',13)
ylabel('$ \ddot{\theta}_{i}$ [$\frac{rad}{s^{2}}$]','Interpreter','latex','FontSize',13)
axis ([0 10 -0.4 1])
 kb0 = [k b];

kb = lsqcurvefit(@(kb, t_sample) curve_fit(kb,t_sample,J1,J2,T0),kb0,...
    sample(:,1),sample(:,2:3));

k = kb(1);
b = kb(2);
x0 = [0, 0, 0, 0];
options = odeset;
[t,x] = ode45(@(t,x) state_ex1(t,x,J1,J2,k,b,T0),tspan,x0,options);
ddx = acc_ex1(x,J1,J2,k,b,T0);
figure(3)
plot(t,ddx(:,1),'LineWidth',1.3)
hold on 
plot(t,ddx(:,2),'LineWidth',1.3)
hold on 
plot(sample(:,1), sample(:,2),'--','LineWidth',1.3)
hold on 
plot(sample(:,1), sample(:,3),'--','LineWidth',1.3)
grid on
legend('$\ddot{\theta}_{1}$','$\ddot{\theta}_{2}$','$\ddot{\theta}_{1} (Sample)$','$\ddot{\theta}_{2} (Sample)$','Interpreter','latex','Location','best',...
    'FontSize',13)
xlabel('time [s]','Interpreter','latex','FontSize',13)
ylabel('$ \ddot{\theta}_{i}$ [$\frac{rad}{s^{2}}$]','Interpreter','latex','FontSize',13)
axis ([0 10 -0.4 1])
figure(4)

ddxin = interp1(t,ddx,sample(:,1));
loglog(sample(:,1), abs(ddxin(:,1)-sample(:,2))) 
xlabel('time [s]','Interpreter','latex','FontSize',13)
ylabel('$ abs(\ddot{\theta}_{1}-\ddot{\theta}_{1} (Sample))$','Interpreter','latex','FontSize',13)
grid on 
figure(5)

ddxin = interp1(t,ddx,sample(:,1));
loglog(sample(:,1), abs(ddxin(:,2)-sample(:,3))) 
xlabel('time [s]','Interpreter','latex','FontSize',13)
ylabel('$ abs(\ddot{\theta}_{2}-\ddot{\theta}_{2} (Sample))$','Interpreter','latex','FontSize',13)
grid on 

fprintf(' b = %8.8f\n', b)
fprintf(' k = %8.8f\n', k)
%% Ex2
clear; clc; close all
%Initial conditions
t0=0; x0=0; v0=0;                               % initial condition
data.t1=1; data.t2=1.5; data.tf=3;                   % time [s]
F=@(x) 1e3+x*120e3;                                      % f0+x*k [N]
%accumul
data.acc.Vi=10e-3;                                       % [m^3]
data.acc.Pi=2.5e6;                                       % [Pa]
data.acc.P0=21e6;                                        % [Pa]
data.acc.gamma=1.2;
%delivery
data.del.ka=1.12;
data.del.kcv=2;
data.del.D23=18e-3;                                      % [m]
data.del.L23=2;                                          % [m]
data.del.f23=0.032;
%distributor
data.dis.kd=12;
data.dis.d0=5e-3;                                        % [m]
%actuator
data.act.Dc=50e-3;                                       % [m]
data.act.Dr=22e-3;                                       % [m]
data.act.m=2;                                            % [kg]
data.act.x_max=200e-3;                                   % [m]

%return
data.ret.D67=18e-3;                                      % [m]
data.ret.L67=15;                                         % [m]
data.ret.f67=0.035;                                      
%tank
data.tan.pt=1e5;                                         % [Pa]
data.tan.Vt0=1e-3;                                       % [m^3]
data.tan.kt=1.12;

%fluid
data.fluid.rho=890;                                        % [kg/m^3]
data.Vn2=data.acc.Pi*data.acc.Vi/data.acc.P0;
data.Va=data.acc.Vi-data.Vn2;
X0=[x0 v0 data.Va data.tan.Vt0]';


tspan=[t0 data.tf];
[tt,yy] = ode23s(@(t,y) ode_ex2(t,y,data,F),tspan,X0);


par_out = zeros(length(tt),17);
YY = zeros(length(tt),4);
u = zeros(length(tt),1);
for ii = 1:length(tt)
    [YY(ii,:),par_out(ii,:)] = ode_ex2(tt(ii),yy(ii,:),data,F);
    u(ii) = command(tt(ii),data);
end


QA = par_out(:,1);
Q1 = par_out(:,2);
Q2 = par_out(:,3);
Q3 = par_out(:,4);
Q4 = par_out(:,5);
Q5 = par_out(:,6);
Q6 = par_out(:,7);
Q7 = par_out(:,8);
PA = par_out(:,9);
P1 = par_out(:,10);
P2 = par_out(:,11);
P3 = par_out(:,12);
P4 = par_out(:,13);
P5 = par_out(:,14);
P6 = par_out(:,15);
P7 = par_out(:,16);
A0 = par_out(:,17);


% plots

    figure(1)
    hold on
    grid on
    plot(tt,data.act.x_max*ones(length(tt),1)*1e3,'--',tt,yy(:,1)*1e3,'LineWidth',1.3)
    xlabel('Time [s]', 'interpreter', 'latex')
    ylabel('x [mm]', 'interpreter', 'latex')
    legend('Maximum stroke','Piston position over time','location','northwest','interpreter', 'latex')
    
    figure(2)
    hold on
    grid on
    plot(tt,YY(:,1),'LineWidth',1.3);
    xlabel('time [s]', 'interpreter', 'latex')
    ylabel('${\dot{x}\ \left[\frac{m}{s}\right]}$', 'interpreter', 'latex')
    legend('Piston velocity', 'interpreter', 'latex')
    
    
    figure(3)
    hold on
    grid on
    plot(tt,yy(:,3)*1e3,tt,yy(:,4)*1e3,'LineWidth',1.3)
    
    xlabel('Time [s]')
    ylabel('Volume [l]')
    legend('Accumulator','Tank')
    
    figure(4)
    hold on
    grid on
    P=plot(tt,QA.*tt*1e3,tt,Q4.*tt*1e3,'--',tt,Q5.*tt*1e3,tt,Q7.*tt*1e3,'--');
    xlabel('Time [s]', 'interpreter', 'latex')
    ylabel('Voumetric flow rate $\frac{l}{s}$', 'interpreter', 'latex')
    legend('$Q_A$','$Q_4$','$Q_5$','$Q_7$', 'interpreter', 'latex')
    P(1).LineWidth = 1.5; P(2).LineWidth = 1.7; P(3).LineWidth = 1.5; P(4).LineWidth = 1.5; 
    
%     figure(5)
%     hold on
%     grid on
%     [ylab,~,~] = plotyy(tt,u,tt,A0*1e6);
%     
%     xlabel('time [s]')
%     ylabel(ylab(1),'Valve Position')
%     ylabel(ylab(2),'valve orifice opening')
 
        figure(6)
    hold on
    grid on
    P=plot(tt,PA*1e-6,'--',tt,P1*1e-6,'--',tt,P2*1e-6,'--',tt,P3*1e-6,'--');
    xlabel('time [s]', 'interpreter', 'latex')
    ylabel('pressure [MPa]', 'interpreter', 'latex')
    legend('$PA$','$P1$','$P2$','$P3$', 'interpreter', 'latex')
    axis([0 3 14 22])
    P(1).LineWidth = 3.5; P(2).LineWidth = 3; P(3).LineWidth = 2.5; P(4).LineWidth = 2;

    figure(7)
    hold on
    grid on
    P=plot(tt,data.tan.pt*ones(size(tt))*1e-6,tt,P6*1e-6,tt,P7*1e-6);
    xlabel('time [s]', 'interpreter', 'latex')
    ylabel('pressure [MPa]', 'interpreter', 'latex')
    legend('$Pt$','$P6$','$P7$', 'interpreter', 'latex')
    P(1).LineWidth = 1.5; P(2).LineWidth = 1.5; P(3).LineWidth = 1.5;
    
    figure(8)
    hold on
    grid on
    plot(tt,P4*1e-6,tt,P5*1e-6,tt,data.acc.P0*ones(length(tt),1)*1e-6,'--', 'LineWidth',1.5)
    xlabel('time [s]', 'interpreter', 'latex')
    ylabel('pressure [MPa]', 'interpreter', 'latex')
    legend('$P4$','$P5$','$P0$', 'interpreter', 'latex')

%3)
options = odeset('event',@pilot_event);
[T,~] = ode23s(@ode_ex2,tspan,X0,options,data,F);
te = T(end);
fprintf(' te = %8.8f\n', te)
%% Ex3
clear; clc; close all;
R1 = 1000; 
R2 = 100;
L = 1e-3;
C = 1e-3;
f = 5; % Hz
Vc0 = 1; % V
tspan = [0 10];

x0 = [Vc0 0];
options = odeset; % Semi-stiff problem  trapezoidal rule ode23t
[t,x] = ode23t(@(t,x) state_ex3(t,x,C,R1,R2,L),tspan,x0,options);
figure(1)

plot(t,x(:,1),'LineWidth',1.2)
grid on 
xlabel('$time$ $[s]$','Interpreter','latex','FontSize',13)
ylabel('$ {V}_{c}$ [$V$]','Interpreter','latex','FontSize',13)
% pt 2
vf = @(t) sin(2*pi*f*t)*atan(t);
dvf = @(t) 2*pi*f*cos(2*pi*f*t)*atan(t) + sin(2*f*pi*t)*1/(1+t^2);
x0 = [Vc0 0];
options = odeset;
[t,x] = ode23t(@(t,x) state_ex3_pt2(t,x,C,R1,R2,L,vf,dvf),tspan,x0,options);
figure(2)
plot(t,x(:,1),'LineWidth',1.2)
grid on
xlabel('$time$ $[s]$','Interpreter','latex','FontSize',13)
ylabel('$ {V}_{c}$ [$V$]','Interpreter','latex','FontSize',13)

%% Ex4
clear;clc;close all;
% DATA
par.L = 8/100;             %[m]
par.r = 40e-2;             %[m]
%par.A = 2*pi*par.r*par.L; 
par.A = 80/1e4;%[m^2]

% Inner liner
par.l1 = 1/100;             %[m]
par.k1 = 200;               %[W/m*K]
par.R1 = par.l1/(par.k1*par.A);
% Conductor 
par.l2 = 5/100;             %[m]
par.k2 = 400;               %[W/m*K]
par.R2 = par.l2/(par.k2*par.A);
rho2 = 8100;                 %[kg/m3]
V2 = par.A * par.l2;
M2 = V2 * rho2;
par.c2 = 385 *M2;             %[J/K]
% Interface
par.R3 = 0.07;             %[m^2*k/W]
% Insulator
par.l4 = 3/100;             %[m]
par.k4 = 160;               %[W/m*K]
rho4 = 8400;               %[kg/m3] 
V4 = par.A * par.l4;
M4 = V4 *rho4;
par.c4 = 456*M4;              %[J/K]
par.R4 = par.l4/(par.k4*par.A); %[m^2*k/W]
% Outer coating
par.l5 = 1/100;             %[m]
par.k5 = 25;               %[W/m*K]
par.R5 = par.l5/(par.k5*par.A); %[m^2*k/W]

% 4.3 Dynamic simulation to show the temperature profiles


par.T0 = 20 + 273.15;          %[K] inside the nozzle at t = 0
par.Ti_0 = par.T0;             %[K] outside fixed 
par.Ti_f = 1000 + 273.15;      %[K] inside the nozzle at t = 60 
par.t1 =linspace(0,1,100);                    %[s]
par.t2 =linspace(1+eps,60,100);                    %[s]
par.tspan = [par.t1 par.t2];                   %[s]
par.Ti = @(t) par.Ti_0 + (par.Ti_f-par.Ti_0)*t;%[K]
par.Tinn = [par.Ti(par.t1) par.Ti_f*ones(1,length(par.t2))]';
par.Tout = par.Ti_0.*ones(length(par.tspan),1);

% Req of the mononodal case
par.Req1 = 0.5*par.R1;
par.Req2 = 0.5*par.R1 + 0.5*par.R2;
par.Req3 = 0.5*par.R2 + 0.5*par.R3;
par.Req4 = 0.5*par.R3 + 0.5*par.R4;
par.Req5 = 0.5*par.R4 + 0.5*par.R5;
par.Req6 = 0.5*par.R5;

[~, TS_1] = ode45(@(tspan,TS) integrator_ex4(tspan,TS,par), [par.t1 par.t2] , par.Ti_0*ones(1,2));
[T1_1, T3_1, T5_1] = evTemp(TS_1(:,1), TS_1(:,2), par, [par.t1 par.t2]);
% Plot of the results
figure(1)
hold on 
grid on
T = [par.Tinn T1_1 TS_1(:,1) T3_1 TS_1(:,2) T5_1 par.Tout];
set(gca, 'FontSize', 10) ;
xlabel('$ t \ [s]$', 'Interpreter', 'Latex'); 
ylabel('$ T \ [K]  $', 'Interpreter', 'Latex') ;

plot(par.tspan,par.Tinn,'LineWidth',1.2)
hold on
plot(par.tspan,T1_1,'LineWidth',1.2)
hold on
plot(par.tspan,TS_1(:,1),'LineWidth',1.2)
hold on
plot(par.tspan,T3_1,'LineWidth',1.2)
hold on
plot(par.tspan,TS_1(:,2),'LineWidth',1.2)
hold on
plot(par.tspan,T5_1,'LineWidth',1.2)
hold on
plot(par.tspan,par.Tout,'LineWidth',1.2)
legend('$N_{Inn}$','$N_{1}$','$N_{2}$','$N_{3}$','$N_{4}$','$N_{5}$','$N_{Out}$','Interpreter','latex','Location','best',...
    'FontSize',13)

% Req of the mononodal case
par.Req1_2 = 0.5*par.R1;
par.Req2_2 = 0.5*par.R1 + 1/3*par.R2;
par.Req3_2 = 1/3*par.R2; 
par.Req4_2 = 1/3*par.R2 + 0.5*par.R3;
par.Req5_2 = 0.5*par.R3 + 1/3*par.R4;
par.Req6_2 = 1/3*par.R4;
par.Req7_2 = 1/3*par.R4 + 0.5*par.R5;
par.Req8_2 = 0.5*par.R5;

[t, TS_2] = ode45(@(tspan,TS) integrator_ex4_pt2(tspan,TS,par), [par.t1 par.t2] , par.Ti_0*ones(1,4));
[T1_2, T4_2, T7_2] = evTemp_pt2(TS_2(:,1), TS_2(:,2),TS_2(:,3),TS_2(:,4), par, [par.t1 par.t2]);
% Plot of the results
figure(2)
hold on 
grid on
set(gca, 'FontSize', 10) ;
xlabel('$  t \ [s]$', 'Interpreter', 'Latex'); 
ylabel('$ T \ [K]  $', 'Interpreter', 'Latex') ;

plot(par.tspan,par.Tinn,'LineWidth',1.2)
hold on
plot(par.tspan,T1_2,'LineWidth',1.2)
hold on
plot(par.tspan,TS_2(:,1),'LineWidth',1.2)
hold on
plot(par.tspan,TS_2(:,2),'Linewidth',1.2)
hold on 
plot(par.tspan,T4_2,'LineWidth',1.2)
hold on
plot(par.tspan,TS_2(:,3),'LineWidth',1.2)
hold on
plot(par.tspan,TS_2(:,4),'LineWidth',1.2)
hold on
plot(par.tspan,T7_2,'LineWidth',1.2)
hold on
plot(par.tspan,par.Tout,'LineWidth',1.2)
legend('$N_{Inn}$','$N_{1}$','$N_{2}$','$N_{3}$','$N_{4}$','$N_{5}$','$N_{6}$','$N_{7}$','$N_{Out}$','Interpreter','latex',...
    'FontSize',13)
%% Functions
function dx = state_ex1(~,x,J1,J2,k,b,T0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             - State space implementation of exercise 1 ODE. 
%
% INPUT
%    x:      state space                             [4x1]
%   J1:      inertia moment disk 1                   [1x1]
%   J2:      inertia moment disk 1                   [1x1]
%    k:      stiffness                               [1x1]
%    b:      viscous friction coefficient            [1x1]
%   T0:      torque                                  [1x1]
%
%
% OUTPUT
%    dx:    integrated state                         [4x1]
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 = x(3);
x2 = x(4);
th1 = x(1);
th2 = x(2);

dx =[ x1
     x2 
   - k/J1 * (th1- th2 )
    T0/J2 - b/J2 * sign(x2)*x2^2 - k/J2 * (th2-th1) ];
end
function ddx = acc_ex1(x,J1,J2,k,b,T0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             - State space implementation of exercise 1 ODE. 
%
% INPUT
%    x:      state space                             [4x1]
%   J1:      inertia moment disk 1                   [1x1]
%   J2:      inertia moment disk 1                   [1x1]
%    k:      stiffness                               [1x1]
%    b:      viscous friction coefficient            [1x1]
%   T0:      torque                                  [1x1]
%
%
% OUTPUT
%   ddx:    integrated state                         [4x1]
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%x1 = x(:,3);
x2 = x(:,4);
th1 = x(:,1);
th2 = x(:,2);

ddx = [- k/J1 * (th1 - th2 ),T0/J2 - b/J2 * sign(x2).*x2.^2 - k/J2 * (th2-th1) ];
end
function q = curve_fit(kb,t_sample,J1,J2,T0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             - Curve fit function used to retrieve k and b from samples. 
%
% INPUT
%         kb:      vector containing stiffness and viscous [1x2]
%            coefficient                             
%   t_sample:      samples                                 [1001x1]
%         J2:      inertia moment disk 2                   [1x1]
%         J1:      inertia moment disk 1                   [1x1]
%         T0:      torque                                  [1x1]
%
%
% OUTPUT
%         q:       solution of k and b                     [1x2]
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = [0, 0, 0, 0];
options = odeset;
[~,s] = ode45(@(t,x) state_ex1(t,x,J1,J2,kb(1),kb(2),T0),t_sample,x0,options);
q = acc_ex1(s,J1,J2,kb(1),kb(2),T0);
end
function [dV] = state_ex3(~,x,C,R1,R2,L)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             - State space implementation of exercise 3 ODE. 
%
% INPUT
%    x:      state space                             [2x1]
%    C:      capacity                                [1x1]
%   R1:      Resistence 1                            [1x1]
%   R2:      Resistence 2                            [1x1]
%    L:      inductance                              [1x1]
%
%
% OUTPUT
%    dV:    integrated state                         [2x1]
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dV= [x(2); -x(1)*(R1/((R1+R2)*C*L))-x(2)*((L+R1*R2*C)/((R2+R1)*C*L))];

end
function [dV] = state_ex3_pt2(t,x,C,R1,R2,L,vf,dvf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             - State space implementation of exercise 3 pt2 ODE. 
%
% INPUT
%    t:      time                                    [1x1]
%    x:      state space                             [2x1]
%    C:      capacity                                [1x1]
%   R1:      Resistence 1                            [1x1]
%   R2:      Resistence 2                            [1x1]
%    L:      inductance                              [1x1]
%   vf:      external voltage function               [1x1]
%  dvf:      derivative of external voltage func     [1x1]
%
% OUTPUT
%    dV:    integrated state                         [2x1]
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = (R1*vf(t) +L*dvf(t))/((R1+R2)*L*C);
dV= [x(2); -x(1)*(R1/((R1+R2)*C*L))-x(2)*((L+R1*R2*C)/((R2+R1)*C*L))+ a];

end
function TS = integrator_ex4(t, temp, par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             - State space implementation of exercise 4 pt1 ODEs. 
%
% INPUT
%    t:      time                                    [1x1]
% temp:      temperature                             [2x1]
%  par:      various parameters                      [1x1 struct]
%   
%
% OUTPUT
%    TS:    integrated state                         [2x1]
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T2 = temp(1);
T4 = temp(2);
[T1, T3, T5] = evTemp(T2,T4,par,t);
dT2 = (1/par.c2)*((T1-T2)/par.Req2 - (T2-T3)/par.Req3);
dT4 = (1/par.c4)*((T3-T4)/par.Req4 - (T4-T5)/par.Req5);
TS = [dT2 dT4]';
end
function TS = integrator_ex4_pt2(t, temp, par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             - State space implementation of exercise 4 pt2 ODEs. 
%
% INPUT
%    t:      time                                    [1x1]
% temp:      temperature                             [2x1]
%  par:      various parameters                      [1x1 struct]
%   
%
% OUTPUT
%    TS:    integrated state                         [2x1]
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T2 = temp(1);
T3 = temp(2);
T5 = temp(3);
T6 = temp(4);
[T1, T4, T7] = evTemp_pt2(T2,T3,T5,T6,par,t);
dT2 = (1/(0.5*par.c2))*((T1-T2)/par.Req2_2 - (T2-T3)/par.Req3_2);
dT3 = (1/(0.5*par.c2))*((T2-T3)/par.Req3_2 - (T3-T4)/par.Req4_2);
dT5 = (1/(0.5*par.c4))*((T4-T5)/par.Req5_2 - (T5-T6)/par.Req6_2);
dT6 = (1/(0.5*par.c4))*((T5-T6)/par.Req6_2 - (T6-T7)/par.Req7_2);
TS = [dT2 dT3 dT5 dT6]';
end
function [T1,T3,T5] = evTemp(T2,T4,par,time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             - Function representng the algebraic equations
%               of exercise 4 pt1. 
%
% INPUT
%   T2:      temperature T2                          [1x1]
%   T4:      temperature T4                          [1x1]
%  par:      various parameters                      [1x1 struct]
% time:      time                                    [1x1]
%   
%
% OUTPUT
%    T1:    temperature T1                           [1x1]
%    T3:    temperature T3                           [1x1]
%    T5:    temperature T5                           [1x1]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tinn = zeros(length(time),1);
Tout = (20+273.15)*ones(length(time),1);

for i = 1:length(time)
    if (time(1,i) <1)
        Tinn(i,1) = (20+273.15) + (1000-20) * time(1,i);
    else
        Tinn(i,1) = 1000 + 273.15;
    end
end
T1 = (Tinn/par.Req1 + T2/par.Req2)/(1/par.Req1 + 1/par.Req2);
T3 = (T2/par.Req3 + T4/par.Req4)/(1/par.Req3 + 1/par.Req4);
T5 = (T4/par.Req5 + Tout/par.Req6)/(1/par.Req5 + 1/par.Req6);
end
function [T1, T4, T7] = evTemp_pt2(T2,T3,T5,T6,par,time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             - Function representng the algebraic equations
%               of exercise 4 pt2. 
%
% INPUT
%   T2:      temperature T2                          [1x1]
%   T3:      temperature T3                          [1x1]
%   T5:      temperature T5                          [1x1]
%   T6:      temperature T6                          [1x1]
%  par:      various parameters                      [1x1 struct]
% time:      time                                    [1x1]
%   
%
% OUTPUT
%    T1:    temperature T1                           [1x1]
%    T4:    temperature T4                           [1x1]
%    T7:    temperature T7                           [1x1]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tinn = zeros(length(time),1);
Tout = (20+273.15)*ones(length(time),1);

for i = 1:length(time)
    if (time(1,i) <1)
        Tinn(i,1) = (20+273.15) + (1000-20) * time(1,i);
    else
        Tinn(i,1) = 1000 + 273.15;
    end
end
T1 = (Tinn/par.Req1_2 + T2/par.Req2_2)/(1/par.Req1_2 + 1/par.Req2_2);
T4 = (T3/par.Req4_2 + T5/par.Req5_2)/(1/par.Req4_2 + 1/par.Req5_2);
T7 = (T6/par.Req7_2 + Tout/par.Req8_2)/(1/par.Req7_2 + 1/par.Req8_2);
end
function [value, isterminal, direction] = pilot_event(~,x,data,~)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             - Control function preventing singularities 
%
% INPUT
%    x:      input                                   [1x1]
% data:      data                                    [1x1 struct]
%
% OUTPUT
%         value:    value                            [1x1]
%    isterminal:    isterminal                       [1x1]
%     direction:    direction                        [1x1]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xn = data.act.x_max; 
    
    value = xn - x(1) ;
    isterminal = 1 ;
    direction = -1 ;
end
function z = command(t,data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             - Command function of the valve 
%
% INPUT
%    t:      input                                   [1x1]
% data:      data                                    [1x1 struct]
%
% OUTPUT
%         z:    output                               [1x1]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if t <= data.t1
        z = 0;
    elseif t > data.t1 && t < data.t2
        z = (t-data.t1)/(data.t2-data.t1);
    elseif t >= data.t2
        z = 1;
    end
end
function [yy,parout] = ode_ex2(t,xx,data,F)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             - State space implementation of exercise 1 ODE. 
%
% INPUT
%    t:      state space                             [4x1]
%   xx:      position                                [1x1]
% data:      data                                    [1x1 struct]
%    F:      stiffness                               [1x1 function handle]
%
%
% OUTPUT
%   yy:    integrated state                           [4x1]
%   parout: output parameter (pressure and flow rate) [18x1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initial contitions
    x = xx(1);
    v = xx(2);
    Va = xx(3);
   % Vt = xx(4);
    
    A4 = (pi*data.act.Dc^2/4);
    A5 = (pi/4*(data.act.Dc^2-data.act.Dr^2));
    
    % actuator boundary conditions
    if x < 0
        x = 0;
    end
    if x <= 0 && v < 0
        x = 0;
        v = 0;
    end
    if x > data.act.x_max
        x = data.act.x_max;
    end
    if x >= data.act.x_max && v > 0
        x = data.act.x_max;
        v = 0;
    end
    
    % control history
    z = command(t,data);
    alpha = 2*acos(1 - 2*z);
    a_0 = data.dis.d0^2/8*(alpha - sin(alpha));
    
    % volume flow rates
    Q_4 = A4*v;
    Q_5 = A5*v;
    if z == 0
        Q_3 = 0;
        Q_6 = 0;
        Q_4 = 0;
        Q_5 = 0;
    end
    if z > 0
        Q_3 = Q_4;
        Q_6 = Q_5;
    end
    
    Q_2 = Q_3;
    Q_1 = Q_2;
    if Q_1 <= 0
        Q_2 = 0;
        Q_3 = Q_2;
        if z > 0
            Q_4 = Q_3;
        end
    end

    Q_A = Q_1;
    Q_7 = Q_6;
    
    % pressures calculation
    p_A = data.acc.P0*(data.Vn2/(data.acc.Vi-Va))^data.acc.gamma;
    p_1 = p_A - data.del.ka*0.5*data.fluid.rho*(Q_1/(pi*data.del.D23^2/4))^2;
    p_2 = p_1 - data.del.kcv*0.5*data.fluid.rho*(Q_2/(pi*data.del.D23^2/4))^2;
    p_3 = p_2 - data.del.f23*data.del.L23/data.del.D23*0.5*data.fluid.rho*(Q_3/(pi*data.del.D23^2/4))^2;
    p_7 = data.tan.pt + data.tan.kt*0.5*data.fluid.rho*(Q_7/(pi*data.ret.D67^2/4))^2;
    p_6 = p_7 + data.ret.f67*data.ret.L67/data.ret.D67*0.5*data.fluid.rho*(Q_7/(pi*data.ret.D67^2/4))^2;
    p_eq0 = (data.tan.pt*A5 + F(0))/A4;

    % control valve instructions
    if z == 0
        p_4 = p_eq0;
        p_5 = data.tan.pt;
    else
        if a_0 < 1e-9
            p_4 = p_eq0;
            p_5 = data.tan.pt;
        else
            p_4 = p_3 - data.dis.kd*0.5*data.fluid.rho*(Q_3/a_0)^2;
            p_5 = p_6 + data.dis.kd*0.5*data.fluid.rho*(Q_6/a_0)^2;
        end
    end
    
    % ODE system
    Va_dot = -Q_A;
    Vt_dot = Q_7;
    x_dot = v;
    v_dot = (p_4*A4 - p_5*A5 - F(x))/data.act.m;

    % actuator phyhsical constraints
    if (x >= data.act.x_max && v_dot>0) || (x <= 0 && v_dot<0)
        v_dot = 0;
    end
    if Va >= data.acc.Vi
        Va_dot = 0;
    end
    
    yy(1,1) = x_dot;
    yy(2,1) = v_dot;
    yy(3,1) = Va_dot;
    yy(4,1) = Vt_dot;
    
    parout = [Q_A Q_1 Q_2 Q_3 Q_4 Q_5 Q_6 Q_7...
            p_A p_1 p_2 p_3 p_4 p_5 p_6 p_7 a_0];
end


