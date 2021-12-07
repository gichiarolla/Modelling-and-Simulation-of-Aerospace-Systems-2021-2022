% Modelling and Simulation of Aerospace Systems (2021/2022)
% Assignment # 1
% Author Giovanni Chiarolla matr. 964831
%% Ex 1 
clearvars; close all; clc;
f = @(x) cos(x)-x;
figure
grid on
fplot(f)
hold on 
yline(0)
grid on
xlabel('x','Interpreter','latex')
ylabel('f(x)','Interpreter','latex')

a = 0.5;
b = 1.0;
err = 1e-8;

for i = 1:10 % Bisection Method
tic;
[root_BI,feval_BI] = bisection(f,a,b,err);
t_BI = toc();


% Secant method
tic
[root_SE,feval_SEC] = secant(f,a,b,err);
t_SEC = toc();


% Regula-falsi Method
tic
[root_RF,feval_RF] = regulafalsi(f,a,b,err);
t_RF = toc();

end
fprintf(' root_BI = %8.8f\n', root_BI)
fprintf(' root_SE = %8.8f \n', root_SE)
fprintf(' root_RF = %8.8f \n', root_RF)

fprintf(' Bisection function evaluations = %4.0f\n', feval_BI) 
fprintf(' Regula-Falsi function evaluations = %4.0f\n', feval_RF) 
fprintf(' Secant function evaluations = %4.0f\n', feval_SEC) 
%% Ex 2
 clearvars; close all; clc;
% Plot of the vector field (to find initial conditions) 
[X,Y] = meshgrid(-2:0.1:2,-2:0.1:2);
eq1 = X.^2 - X - Y;
eq2 =  X.^2/16 + Y.^2 - 1;
q=quiver(X,Y,eq1,eq2,1.5);

 syms  x1 x2

f =  [x1^2-x1-x2, x1^2/16 + x2^2 - 1];
% Initial conditions 
x0 = [3,2];
x02 = [-2,1];
% Newton's method tolerance and max nmax iterations
tol = 1e-8;
nmax = 6;
% Analytical jacobian 
[x_an(:,1),it_an] = newton_f_4(f,x0, tol,x1,x2,nmax);
[x_an(:,2),it_a] =  newton_f_4(f,x02, tol,x1,x2,nmax);

% function
f = @(x1,x2) [x1^2-x1-x2; x1^2/16 + x2^2 - 1];
% forward
[x_fd(:,1),it_fd] = newton_finite(x0', nmax, tol, f, 1);
[x_fd(:,2),it_f] = newton_finite(x02', nmax, tol, f, 1);
% centred
[x_cd(:,1),it_cd] = newton_finite(x0', nmax, tol, f, 0);
[x_cd(:,2),it_c] = newton_finite(x02', nmax, tol, f, 0);
% Compute exact roots
syms x y real
[x,y] = solve(x^2-x-y==0, x^2/16 +y^2-1==0);
SOL = double([x';y']);
% errors
 error_an = vecnorm(x_an - SOL);
 error_fd = vecnorm(x_fd - SOL);
 error_cd = vecnorm(x_cd - SOL);
 

%% Ex 3
clearvars; close all; clc;
% heun method
x0 = 1/2;
ti = 0;
tf = 2;
h_m = [0.5, 0.2, 0.05, 0.01] ;
f = @(x,t) x - t^2 + 1;
fa  = @( t ) t.^2 +2*t+1-1/2* exp(t);
xsolhe = zeros(length(h_m),1);
xsolrk =  zeros(length(h_m),1);
timehe = zeros(1,length(h_m));
timerk4 = zeros(1,length(h_m));
err_he = zeros(1,length(h_m));
err_rk4 = zeros(1,length(h_m));



figure (1)
hold on

for i = 1: length(h_m)
hi = h_m(i);
tic
[xhe,t_he,feval_he] = heun_m(x0,f,ti,tf,hi);


err_he(i) = max(abs(xhe'-fa(ti:h_m(i):tf)));
w = @() heun_m(x0,f,ti,tf,h_m(i));
    timehe(i) = timeit(w);
plot(t_he,xhe,'--','LineWidth',0.5)
hold on 



end
fplot(fa,'LineWidth',1.5)
grid on

legend('h=0.5','h=0.2','h=0.05','h=0.01','analytical','Interpreter','latex','Location','best',...
    'FontSize',13)
xlabel('$time [s]$','Interpreter','latex','FontSize',13)
ylabel('$x_{Heun}(t)$','Interpreter','latex','FontSize',13)
axis([0 2 0 5.5])
% Runge Kutta 4
figure(2)
hold on
for i = 1 : length(h_m)

[xrk4, t_rk4, feval_rk4] = rk4(x0,f,ti,tf,h_m(i));

xsolhe(i) = xhe(end);
xsolrk(i) = xrk4(end);
err_rk4(i) = max(abs(xrk4'-fa(ti:h_m(i):tf)));
w = @() rk4(x0,f,ti,tf,h_m(i));
    timerk4(i) = timeit(w);
plot(t_rk4,xrk4,'--','LineWidth',0.5)
hold on

end
fplot(fa,'LineWidth',1.5)
grid on

legend('h=0.5','h=0.2','h=0.05','h=0.01','analytical','Interpreter','latex','Location','best',...
    'FontSize',13)
xlabel('$time [s]$','Interpreter','latex','FontSize',13)
ylabel('$x_{RK4}(t)$','Interpreter','latex','FontSize',13)
axis([0 2 0 5.5])


figure(3)
fig_RK2=loglog(timehe,err_he);
fig_RK2.LineWidth=2;
grid on
xlabel('Integration time [s]','Interpreter','latex','FontSize',13)
ylabel('Solution error','Interpreter','latex','FontSize',13)

figure(4)
fig_RK4=loglog(timerk4,err_rk4);
fig_RK4.LineWidth=2;
xlabel('Integration time [s]','Interpreter','latex','FontSize',13)
ylabel('Solution error','Interpreter','latex','FontSize',13)
grid on
%% Ex 4
% (It takes 30 seconds)
clearvars; close all; clc; 
A = @(alpha) [0,1;-1,2*cos(alpha)];
% Operator Frk2
I = eye(2);
Frk2 = @(alpha,h) (I + [0,1;-1,2*cos(alpha)]*h + h^2/2*[0,1;-1,2*cos(alpha)]^2);
h_max = fzero(@(h) (max(abs(eig(Frk2(pi,h))))-1),0);
alpha = pi:-0.01:0;
a = zeros(2,length(alpha));
b = zeros(2,length(alpha));
h_max_col = zeros(1,length(alpha));
figure(1)
hold on
for i = 1: length(alpha)  
f = @(h) max(abs(eig(Frk2(alpha(i),h))));
fp = fplot(f); % plot of the solutions of the problem. Used to find proper initial conditions in following fzeros 
fp.Color = [alpha(i), 0 , 0]/pi; % varia la tonalita di rosso 
if alpha(i) < pi/2 + 0.025 && alpha(i) > pi/2 -0.025
    fp.Color = [0, 0 , 1];
    fp.LineWidth = 2; 
end    
if alpha(i)<pi/2    
h_max_col(i) = fzero(@(h) max(abs(eig(Frk2(alpha(i),h))))-1,0);
else 
 h_max_col(i) = fzero(@(h) max(abs(eig(Frk2(alpha(i),h))))-1,alpha(i));   
a(:,i) = eig(A(alpha(i)));
b(:,i) = a(:,i)*h_max_col(i);
end
end
grid on
% Solution for alpha equal to pi (Heun)
h1 = fzero(@(h) max(abs(eig(Frk2(pi,h))))-1,pi);
xlim([-1 4])
ylim([-1 5])
line(xlim(),[1,1], 'LineWidth',0.8,'Color',[0 0 0 ])
line([0,0],ylim(), 'LineWidth',0.8,'Color',[0 0 0 ])
figure(2)

hold on 

area(real(b(1,:)),imag(b(1,:)),'FaceColor','r','FaceAlpha',.3,'EdgeAlpha',1)
plot(real(b(1,1)),imag(b(1,1)),'.','Color',[1 0 0],'MarkerSize',20)
area(real(b(2,:)),imag(b(2,:)),'FaceColor','r','FaceAlpha',.3,'EdgeAlpha',1)
axis([-3 3 -3 3])
grid on
line(xlim(),[0,0], 'LineWidth',0.5,'Color',[0 0 0 ])
line([0,0],ylim(), 'LineWidth',0.5,'Color',[0 0 0 ])
legend('Heun Stability  Region',['$max(|eig(F(h,\alpha)|) = 1$' newline '  for $\alpha = pi$'],'Interpreter',...
    'latex','Location','best','FontSize',8)
xlabel('$Re\{h\lambda\}$','Interpreter','latex')
ylabel('$Im\{h\lambda\}$','Interpreter','latex')

% Operator Frk4
Frk4 = @(A,h) I + h*A +1/2*h^2*A^2 + 1/6*h^3*A^3 + 1/24*h^4*A^4;
figure(3)
hold on

h_max_col_a = zeros(1,length(alpha));
b_a = zeros(2,length(alpha));  
 for i = 1 : length(alpha)
    f = @(h) max(abs(eig(Frk4(A(alpha(i)),h))));
    fp = fplot(f);
fp.Color = [abs(alpha(i)), 0 , 0]/pi; % varia la tonalita di rosso 
if alpha(i) < pi/2 + 0.025 && alpha(i) > pi/2 -0.025
    fp.Color = [0, 0 , 1];
    fp.LineWidth = 2;
    
end  
h0 = 3;
h_max_col_a(i) = fzero(@(h) max(abs(eig(Frk4(A(alpha(i)),h))))-1,alpha(i));
h_max_col(i) = fzero(@(h) max(abs(eig(Frk4(A(alpha(i)),h))))-1,h0);
a(:,i) = eig(A(alpha(i)));
b(:,i) = a(:,i)*h_max_col(i);
b_a(:,i) = a(:,i)*h_max_col_a(i);
w =abs( b(1,i));
 end
 % Solution for alpha equal to pi (RK4)
 h2 = fzero(@(h) max(abs(eig(Frk4(A(pi),h))))-1,h0);
grid on
xlim([-1 4])
ylim([-1 5])
line(xlim(),[1,1], 'LineWidth',0.8,'Color',[0 0 0 ])
line([0,0],ylim(), 'LineWidth',0.8,'Color',[0 0 0 ])
 figure(4)


hold on
plot(real(b(1,1)),imag(b(1,1)),'.','Color',[0 0 0],'MarkerSize',20);
h_m = [0.5, 0.2, 0.05, 0.01];
hold on
scatter(h_m, zeros(1, length(h_m)),30,[1 0 1],'*')
hold on 
plot(real(b(1,:)),imag(b(1,:)),'LineStyle','none','Marker','o','MarkerSize',...
    1.5,'MarkerFaceColor','r','Color','r')
plot(real(b(2,:)),imag(b(2,:)),'LineStyle','none','Marker','o','MarkerSize',...
    1.5,'MarkerFaceColor','r','Color','r')


hold on
plot(real(b_a(1,:)),imag(b_a(1,:)),'LineStyle','none','Marker','o','MarkerSize',...
    1.5,'MarkerFaceColor','r','Color','r')
plot(real(b_a(2,:)),imag(b_a(2,:)),'LineStyle','none','Marker','o','MarkerSize',...
    1.5,'MarkerFaceColor','r','Color','r')

grid on
text(-2.5,0.2,'RK4 Stability Region','FontSize',10,'Color', 'r')
axis([-4 4 -3.5 3.5])


line(xlim(),[0,0], 'LineWidth',0.5,'Color',[0 0 0 ])
line([0,0],ylim(), 'LineWidth',0.5,'Color',[0 0 0 ])
legend([' $max(|eig(F(h,\alpha)|) = 1$' newline ' for $\alpha = pi$'],...
    '$h_{i}\lambda$  of Exercise 3','Interpreter',...
    'latex','Location','best','FontSize',8)
xlabel('$Re\{h\lambda\}$','Interpreter','latex')
ylabel('$Im\{h\lambda\}$','Interpreter','latex')

fprintf(' h for alpha equal to pi (Heun) = %8.8f\n', h1)
fprintf(' h for alpha equal to pi (RK4) = %8.8f \n', h2)

%% Ex 5 
%(Unfortunately, it takes 5 minutes) 
clearvars; close all; clc
A = @(alpha) [0,1;-1,2*cos(alpha)]; % A(alpha) operator
x0 = [1;1]; % Initial condition 
fevals = zeros(3,4); %Static allocation of the function evaluation matrix
t = [0 1]; % time range 
alpha = linspace(pi, 0,100); % vector containing alpha variaiton from pi to 0
x_an = @(a) expm(A(a))*x0;  % Analytic solution     
tol = [1e-3 1e-4 1e-5 1e-6]; % Tolerances

lambda = length(alpha); % Static allocation of the eigenvalues
for i=1:length(alpha)
    lambda(i) = max(eig(A(alpha(i))));  
end
guess = [5*tol;                  % Initial value of the fzero function 
    0.01 0.01 0.01 0.01; %0.01
    0.5 0.5 0.5 0.5]; % 0.5
tic
for i=1:3
    figure(i)
    hold on
    axis equal
    for j=1:length(tol)
        options = optimset('TolX',1e-2*tol(j));   % change the tolerance of fzero to decrease cpu time
        guessrk1 = tol(j);
        for k=1:length(alpha)
            a = alpha(k);
            if i==1
                h_m(k) = fzero(@(x) norm(switch_method(@(t,h) A(a)*h,t,x0,abs(x),i) ...
                    - x_an(a),Inf)-tol(j),guessrk1,options);
                if isnan(h_m(k))     % condition added to avoid singularities
                    guessrk1=tol(j);
                else
                    guessrk1=h_m(k);
                end
            else
                h_m(k) = fzero(@(x) norm(switch_method(@(t,h) A(a)*h,t,x0,abs(x),i) ...
                    - x_an(a),inf) - tol(j),guess(i,j));
            end
        end
        
        [~,fevals(i,j)] = switch_method(@(t,h) A(pi)*h,t,x0,abs(h_m(1)),i); %  number of function evaluation
        
        h = abs(h_m).*lambda;                
        plot([real(h) flip(real(h))],...
            [imag(h) -flip(imag(h))],'LineWidth',1.5)  % plot of the solution in complex plane
    end
    grid on
    legend('tol = $10^{-3}$','tol = $10^{-4}$','tol = $10^{-5}$','tol = $10^{-6}$','Interpreter','latex','Location','best')
    xlabel('$Re\{h\lambda\}$','Interpreter','latex')
    ylabel('$Im\{h\lambda\}$','Interpreter','latex')
    
end

figure(4)
loglog(tol,fevals,'LineWidth',1.5)
grid on
xlabel('tol','Interpreter','latex','FontSize',13)
ylabel('Function evaluations','Interpreter','latex','FontSize',13)
legend('RK1','RK2','RK4','Interpreter','latex','FontSize',13)
t_ex5=toc();
%% Ex 6 
clear; close all; clc;
alpha = linspace( pi,0,100);
A = @(alpha) [ 0 1; -1 2*cos(alpha)];
I = eye(2);
% BI2th Operator 
F = @(A,h,th)(I + A*th*h + 0.5*(A*th*h)^2)/(I - A*(1-th)*h + 0.5*(A*(1-th)*h)^2);
a = zeros(length(alpha),1);
h_m = zeros(length(alpha),1);
th = 0.4; % theta
for i = 1: length(alpha)
     a(i) = max(eig(A(alpha(i))));
     
     % try/catch statement added to avoid singularities
try 
    f = @(x) max(abs(eig(F(A(alpha(i)),x,th))))-1;
    h_m(i) = fzero(f,9);
  catch mess
            warning(mess.message) 
end   
end
% Plots
figure(1) 
hold on 
grid on

h1 = abs(h_m).*a;
X = [real(h1); flip(real(h1))];
Y = [imag(h1); -flip(imag(h1))];
plot(X,Y,'r','LineWidth',2)
xline(0)
yline(0)
axis equal
axis padded
legend('$BI2_{0.4}$','Interpreter','latex','Location','best')
xlabel('$Re\{h\lambda\}$','Interpreter','latex')
ylabel('$Im\{h\lambda\}$','Interpreter','latex')
% Second part 
theta = [0.1 0.3 0.7 0.9];
for k = 1:length(theta) 
    th = theta(k);
    for i = 1:length(alpha)
        try 
    f = @(x) max(abs(eig(F(A(alpha(i)),x,th))))-1;
    h_m(i) = fzero(f,4);
  catch mess
            warning(mess.message) 
        end
    end
figure(2) 
hold on 
grid on


h1 = abs(h_m).*a;
X = [real(h1); flip(real(h1))];
Y = [imag(h1); -flip(imag(h1))];
plot(X,Y,'LineWidth',2)

end
fig1 = xline(0);
fig2 = yline(0);
axis equal
axis padded
legend_entries = {'$BI2_{0.1}$','$BI2_{0.3}$','$BI2_{0.7}$','$BI2_{0.9}$'};
legend(legend_entries{:},'Interpreter','latex','Location','best')
set(get(get(fig1,'Annotation'),'LegendInformation'), ...
'IconDisplayStyle','off') % where H is the name of the plot
xlabel('$Re\{h\lambda\}$','Interpreter','latex')
ylabel('$Im\{h\lambda\}$','Interpreter','latex')
%% Ex 7
clear; clc; close all;
B = [ -180.5 219.5; 179.5 -220.5];
x0 = [1 1 ]';
h_ = 0.1;
I = eye(2);
% analytical solution 
x = @(t) expm(B*t)*x0;
t = [0 5];
t = t(1):h_:t(end);
% FR4 Operator 
Frk4 = @(A,h) I + h*A +1/2*h^2*A^2 + 1/6*h^3*A^3 + 1/24*h^4*A^4;
% BI2th Operator
BI2th = @(A,h,th)(I + A*th*h + 0.5*(A*th*h)^2)/(I - A*(1-th)*h + 0.5*(A*(1-th)*h)^2);
x_frk4 = zeros(2,51);
for i = 1:length(t)
    n = i-1;
x_frk4(:,i) = Frk4(B,h_)^n*x0;
x_an(:,i) = x(t(i));
end
% Plots
figure(1)
plot(t,x_an,'LineWidth',1)
hold on 
plot(t,x_frk4,'LineWidth',1)
grid on

legend('$x_1 - analytical$','$x_2 - analytical$',...
                    '$x_1$  $RK4$','$x_2$  $RK4$','Interpreter','latex','Location','best',...
    'FontSize',13)
xlabel('$time [s]$','Interpreter','latex','FontSize',13)
ylabel('$x_1$ and $x_2$','Interpreter','latex','FontSize',13)
% BI20.1
th = 0.1;
x_BI2 = zeros(2,51);
for i = 1:length(t)
    n = i-1;
x_BI2(:,i) = BI2th(B,h_,th)^n*x0;
x_an(:,i) = x(t(i));
end
figure(2)
plot(t,x_BI2,'LineWidth',1)
hold on 
plot(t,x_an,'LineWidth',1)
grid on
legend('$x_1$ $analytical$','$x_2$  $analytical$',...
                    '$x_1$  $BI2_{0.1}$','$x_2$  $BI2_{0.1}$','Interpreter','latex','Location','best',...
    'FontSize',13)
xlabel('$time [s]$','Interpreter','latex','FontSize',13)
ylabel('$x_1$ and $x_2$','Interpreter','latex','FontSize',13)
alpha = linspace( pi,0,100);
A = @(alpha) [ 0 1; -1 2*cos(alpha)];
I = eye(2);
F = @(A,h,th) (I + A*th*h + 0.5*(A*th*h)^2)/(I - A*(1-th)*h + 0.5*(A*(1-th)*h)^2);
Frk4 = @(A,h) I + h*A +1/2*h^2*A^2 + 1/6*h^3*A^3 + 1/24*h^4*A^4;
a = zeros(length(alpha),1);
h_m = zeros(length(alpha),1);
th = 0.1;
for i = 1: length(alpha)
    
    a(i) = max(eig(A(alpha(i))));
try 
    f = @(x) max(abs(eig(F(A(alpha(i)),x,th))))-1;
    h_m(i) = fzero(f,5);
  catch mess
            warning(mess.message) 
end

    
end
 ar = zeros(2,length(alpha));
 br = zeros(2,length(alpha));
 for i = 1 : length(alpha)
h0 = 3;
h_max_col(i) = fzero(@(h) max(abs(eig(Frk4(A(alpha(i)),h))))-1,h0);
ar(:,i) = eig(A(alpha(i)));
br(:,i) = ar(:,i)*h_max_col(i);

 end
figure(3)
a_  = h_*eig(B);
plot(real(a_),imag(a_),'*m','MarkerSize',10)
hold on 
h1 = abs(h_m).*a;
plot(real(br(1,:)),imag(br(1,:)),'r','LineWidth',2)
hold on
X = [real(h1); flip(real(h1))];
Y = [imag(h1); -flip(imag(h1))];
plot(X,Y,'b','LineWidth',2)
hold on 
plot(real(br(2,:)),imag(br(2,:)),'r','LineWidth',2)
grid on
axis equal
line(xlim(),[0,0], 'LineWidth',0.5,'Color',[0 0 0 ])
line([0,0],ylim(), 'LineWidth',0.5,'Color',[0 0 0 ])
legend('$h_{i}\lambda$ associated to the IVP','$RK4$','$BI2_{0.1}$','Interpreter','latex','Location','best')
xlabel('$Re\{h\lambda\}$','Interpreter','latex')
ylabel('$Im\{h\lambda\}$','Interpreter','latex')


%% Functions 

function [root,feval] = bisection(f,a,b,err)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             - Implementation of the bisection method knowing the 
%               initial value and the final value of the interval.
%
% INPUT
%     f:      function                                [function handle]
%     a:      initial value of the interval           [1x1]
%     b:      final value of the interval             [1x1]
%   err:      accuracy                                [1x1]
%
% OUTPUT
%     a:    zero of the function                      [1x1]
%     b:    zero of the function                      [1x1]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
feval = 0;
% Bisection Method
while abs(b-a)/2>err
 c = (a+b)/2;
 
if f(a)*f(c) < 0 
           b = c;
else
           a = c;
end 
feval = feval +2;
end
root = c;
end
function [root,feval] = secant(f,a,b,err)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             - Implementation of the secant method knowing the 
%               initial value and the final value of the interval.
%
% INPUT
%     f:      function                                [function handle]
%     a:      initial value of the interval           [1x1]
%     b:      final value of the interval             [1x1]
%   err:      accuracy                                [1x1]
%
% OUTPUT
%     a:    zero of the function                      [1x1]
%     b:    zero of the function                      [1x1]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = a; 
x1 = b;
feval = 0;

 while abs(x1-x0)/2 > err
 x2 = x1 - f(x1)*(x1 - x0)/( f(x1) - f(x0) );
 feval = feval +2;
 x0 = x1;
 x1 = x2;
 
 end
 root = x2;
end
function [root,feval] = regulafalsi(f,a,b,err)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             - Implementation of the regula-falsi method knowing the 
%               initial value and the final value of the interval.
%
% INPUT
%     f:      function                                [function handle]
%     a:      initial value of the interval           [1x1]
%     b:      final value of the interval             [1x1]
%   err:      accuracy                                [1x1]
%
% OUTPUT
%     a:    zero of the function                      [1x1]
%     b:    zero of the function                      [1x1]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

feval = 0;
while abs(b-a)/2>err
  c = b - (b-a)*f(b)/(f(b) - f(a)); 
% c = (a*f(b)-b*f(a))/(f(b)-f(a));
feval = feval + 2;
if f(c) > 0 
           b = c;
else
           a = c;
end
feval = feval +3;
end
root = c;
end
function x = realfun(x)

f_s =  [x(1).^2 - x(1) - x(2), x(1).^2/16 + x(2).^2 - 1];
x = f_s;
end

function [x,it] = newton_f_4 ( f ,x0, tol,y,z,nmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             - Implementation of the newton's method with analytical 
%               evaluation of the function derivative.
%
% INPUT
%     f:      function                                [1x2 sym]
%    x0:      initial value of the interval           [1x2]
%     y:      syms value                              [1x1 sym]
%     z:      syms value                              [1x1 sym]
%
% OUTPUT
%     x:    zeros of the function                     [1x2]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




k = 1;
t = zeros(1,3);
t(1) = tol + 1;

it = 0;
if abs(double(subs(f, [y, z], x0))) > 0 
while (it < nmax)


    

jnumeric = double(subs(jacobian( f, [y,z]),[y,z],x0));
fnumeric = double(subs(f,[y,z],x0));
deltax = - inv(jnumeric)*double(subs(f,[y,z],x0))';
x0 = deltax' + x0 ; 


t(k+1) = norm(fnumeric);
k = k+1;
it = it +1;
 end
else 
    error('using this x0: f = 0')
 end

x = x0;
end
function [xvect,it] = newton_finite(x0, nmax, tol, f, method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             - Implementation of the Newton's method finite difference method
%               knowing the initial value x0. 
%
% INPUT
%    x0:      initial value of the interval           [1x2]
%  nmax:      max steps                               [1x1]
%   tol:      tolerance                               [1x1]
%     f:      function                                [1x2 handle fun]
%
%
%
%
% OUTPUT
%     x:    zeros of the function                     [1x2]
%    it:    iterations                                [1x1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  err = tol + 1;
     it = 0;
     xvect = x0;
     
     if method == 0
             % Centered differences  
              
     else
            % Forward differences               
     end
     
     while (it < nmax && err >= tol)
         xv = xvect;
         
         h = max(sqrt(eps), sqrt(eps)*abs(xv));
         h1 = 2 * h(1);
         h2 = 2 * h(2);
         
         if method == 0
             % Centered differences approximated derivative 
               dfdx1 = @(x1, x2) ((feval(f, x1 + h(1), x2)) - ...
                   (feval(f, x1 - h(1), x2))) / h1;
               dfdx2 = @(x1, x2) ((feval(f, x1, x2 + h(2))) - ...
                   (feval(f, x1, x2 - h(2)))) / h2;
         else
             % Forward differences approximated derivatives
               dfdx1 = @(x1, x2) ((feval(f, x1 + h(1), x2)) - ...
                   (feval(f, x1, x2))) * 2 / h1;
               dfdx2 = @(x1, x2) ((feval(f, x1, x2 + h(2))) - ...
                   (feval(f, x1, x2))) * 2 / h2;
         end
         
       % Approximated Jacobian
         J = @(x1, x2) [ feval(dfdx1, x1, x2), feval(dfdx2, x1, x2)];
         
         xn = xv - feval(J, xv(1), xv(2))\feval(f, xv(1), xv(2));
         err = norm(xn - xv, Inf);
         xvect = xn;
         it = it + 1;
     end
     for i = 2 : length(xvect)
     end
   end
function [x, t,feval] = heun_m(x0,f,ti,tf,h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             - Implementation of the Heun's method
%
% INPUT
%    x0:      initial value of the interval           [1x2]
%     f:      function                                [1x2 handle fun]
%    ti:      initial time                            [1x1 sym]
%    tf:      final time                              [1x1 sym]
%
% OUTPUT
%     x:    zeros of the function                     [1x2]
%     t:    time                                      [1x1]
% feval:    function evaluations                      [1x1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = (tf-ti)/h;
t = linspace(ti, tf, n+1);
xp = zeros(n+1,1);
xc = zeros(n+1,1);
xc(1) = x0;
for k = 1:n
% Predictor 
xp(k+1) = x0(k) + h*f(x0(k),t(k));

% Corrector 
xc(k+1) = x0(k) + h/2 * (f(x0(k),t(k)) + f(xp(k+1),t(k+1)));
x0(k+1) = xc(k+1);
end
feval = 2*k;

x = xc;
end
function [x, t,feval] = rk4(x0,f,ti,tf,h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             - Implementation of the RK4 method
%
% INPUT
%
%    x0:      initial value of the interval           [1x2]
%     f:      function                                [1x2 handle fun]
%    ti:      initial time                            [1x1]
%    tf:      final time                              [1x1]
%     h:      step                                    [1x1]
%
% OUTPUT
%     x:    zeros of the function                     [1x2]
%     t:    time                                      [1x1]
% feval:    function evaluation                       [1x1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = (tf-ti)/h;
t = linspace(ti, tf, n+1);
xp1 = zeros(n+1,1);
xp2 = zeros(n+1,1);
xp3 = zeros(n+1,1);


xc = zeros(n+1,1);
xc(1) = x0;
for k = 1:n
% Predictor 
xp1(k+1) = x0(k) + h/2*f(x0(k),t(k));
xp2(k+1) = x0(k) + h/2*f(xp1(k+1),t(k)+1/2*h);
xp3(k+1) = x0(k) + h*f(xp2(k+1),t(k)+ 1/2 * h);
% Corrector 
xc(k+1) = x0(k) + h * (1/6*f(x0(k),t(k)) + 1/3*f(xp1(k+1),t(k)+1/2*h)...
    + 1/3*f(xp2(k+1),t(k)+1/2*h) + 1/6 *f(xp3(k+1) , t(k)+h));
x0(k+1) = xc(k+1);
end
feval = 4*k;
x = xc;
end

function [x,fevals,times] = switch_method(f,t,x0,h,method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             - Switch function among RK1, RK2, and RK3
%
% INPUT
%     f:      function                                [1x2 handle fun]
%     t:      initial time                            [1x1]
%    x0:      initial value of the interval           [1x2]
%     h:      step                                    [1x1]
%method:      the method                              [1x1]
%
% OUTPUT
%     x:    zeros of the function                     [1x2]
% feval:    function evaluations                      [1x1]
%     t:    time                                      [1x1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch method
    case 1
        [x,fevals,times] = RK1(f,t,x0,h);
    case 2
        [x,fevals,times] = RK2(f,t,x0,h);
    case 3
        [x,fevals,times] = RK4(f,t,x0,h);
end
x = x(:,end);
end
function [x,fevals,times] = RK1(f,tspan,x0,h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             - Implementation of RK1 method.
%
% INPUT
%     f:      function                                [1x2 handle fun]
% tspan:      time                                    [1x1]
%    x0:      initial value of the interval           [1x2]
%     h:      step                                    [1x1]
%
% OUTPUT
%     x:    zeros of the function                     [1x2]
% feval:    function evaluations                      [1x1]
%     t:    time                                      [1x1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = x0;
times = tspan(1):h:tspan(end);

t = tspan(1);

if mod(diff(tspan),h)~=0
    times=[times tspan(end)];
end

hs = diff(times);
for i=1:(length(times)-1)
    
    h = hs(i);
    fk = f(t,x);
   
    x = x + h*fk;
    t = t + h;
end
fevals = i;
end
function [x,fevals,times] = RK2(f,tspan,x0,h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             - Implementation of RK2 method.
%
% INPUT
%     f:      function                                [1x2 handle fun]
% tspan:      time                                    [1x1]
%    x0:      initial value of the interval           [1x2]
%     h:      step                                    [1x1]
%
% OUTPUT
%     x:    zeros of the function                     [1x2]
% feval:    function evaluations                      [1x1]
%     t:    time                                      [1x1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = x0;
times = tspan(1):h:tspan(end);
t = tspan(1);
if mod(diff(tspan),h)~=0
    times=[times tspan(end)];
end

hs = diff(times);

for i=1:(length(times)-1)
    xk = x(:,i);
    
    fk = f(t,xk);
    
    % Shrink the last timestep
    h = hs(i);
    
    % Predictor step
    xp = xk + h*fk;
    t = t + h;
    fp = f(t,xp);
     
    % Corrector step
    x(:,i+1) = xk + h*(0.5*fk + 0.5*fp);
    
end
fevals = 2*i;
end
function [x,fevals,times] = RK4(f,tspan,x0,h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             - Implementation of RK4 method.
%
% INPUT
%     f:      function                                [1x2 handle fun]
% tspan:      time                                    [1x1]
%    x0:      initial value of the interval           [1x2]
%     h:      step                                    [1x1]
%
% OUTPUT
%     x:    zeros of the function                     [1x2]
% feval:    function evaluations                      [1x1]
%     t:    time                                      [1x1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = x0;
times = tspan(1):h:tspan(end);
t = tspan(1);
if mod(diff(tspan),h)~=0
    times=[times tspan(end)];
end

hs = diff(times);

for k=1:(length(times)-1)
    h = hs(k);
    
    xk = x(:,k);
    
    % Predictor steps
    % #1
    k1 = f(t,xk);
    
    % #2
    t2 = t + h/2;
    x2= xk + k1*h/2;
    k2 = f(t2,x2);
    
    % #3
    t3 = t + h/2;
    x3 = xk + k2*h/2;
    k3 = f(t3,x3);
    
    % #4
    t4 = t + h;
    x4 = xk + k3*h;
    k4 = f(t4,x4);
    
    % Corrector step
    x(:,k+1) = xk + (k1 + 2*k2 + 2*k3 + k4)*h/6;
    
    t = t+h;
    
    
end
fevals = 4*k;
end