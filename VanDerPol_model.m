
%%  this code is modified from the koopman mpc £¨automatica 2018£©

clear all
close all
addpath('./Resources')
%rng(2141444)


%% *************************** Dynamics ***********************************

f_u =  @(t,x,u,w)(-[ -2*x( 2,:) ; 0.8*x(1,:) + 10*x(1,:).^2.*x(2,:) - 2*x(2,:) + u]+w );
%f_u =  @(t,x,u)(-[ -1*x( 2,:) ; 1*x(1,:) + 1*x(1,:).^2.*x(2,:) - 1*x(2,:) + u] );

n = 2;
m = 1; % number of control inputs


%% ************************** Discretization ******************************

deltaT = 0.01;
%Runge-Kutta 4
k1 = @(t,x,u,w) (  f_u(t,x,u,w) );
k2 = @(t,x,u,w) ( f_u(t,x + k1(t,x,u,w)*deltaT/2,u,w) );
k3 = @(t,x,u,w) ( f_u(t,x + k2(t,x,u,w)*deltaT/2,u,w) );
k4 = @(t,x,u,w) ( f_u(t,x + k1(t,x,u,w)*deltaT,u,w) );
f_ud = @(t,x,u,w) ( x + (deltaT/6) * ( k1(t,x,u,w) + 2*k2(t,x,u,w) + 2*k3(t,x,u,w) + k4(t,x,u,w)  )   );

%% ************************** Basis functions *****************************

basisFunction = 'rbf';
% RBF centers
Nrbf =2;
cent = 2*(rand(n,Nrbf)*1 - 0.5);
%cent=[0.381196454132004,0.267060960957992;-0.341397958310533,-0.889187981816313];
rbf_type = 'thinplate'; 

% Lifting mapping - RBFs + the state itself
liftFun = @(xx)( [xx;rbf(xx,cent,rbf_type)]-[[0;0];rbf([0;0],cent,rbf_type)]);
%liftFun = @(xx)( [xx;xx(1,:).*xx(2,:).^2;xx(1,:).^2.*xx(2,:)] );
% liftFun = @(xx)( [xx;xx(1,:).*xx(2,:).^2;xx(1,:).^2.*xx(2,:);...
%     xx(1,:).*xx(1,:);xx(1,:).*xx(1,:).^2;xx(1,:).^2.*xx(1,:);...
% xx(2,:).*xx(2,:);xx(2,:).*xx(2,:).^2;xx(2,:).^2.*xx(2,:);xx(1,:).*xx(2,:);xx(1,:).*xx(2,:).^3] );
% liftFun = @(xx)( [xx;xx(1,:).*xx(2,:).^2;xx(1,:).^2.*xx(2,:);...
%     xx(1,:).*xx(1,:);xx(1,:).*xx(1,:).^2;xx(1,:).^2.*xx(1,:);...
% xx(2,:).*xx(2,:);xx(2,:).*xx(2,:).^2;xx(2,:).^2.*xx(2,:);xx(1,:).*xx(2,:);xx(1,:).*xx(2,:).^3;...
% xx(1,:).*xx(2,:).^4;xx(1,:).^3.*xx(2,:);...
%     xx(1,:).*xx(1,:).^3;xx(1,:).^3.*xx(1,:).^2;xx(1,:).^4.*xx(1,:);...
% xx(1,:).*xx(2,:).^4;xx(2,:).*xx(2,:).^4;xx(2,:).^3.*xx(1,:);xx(1,:).^3.*xx(2,:);xx(1,:).^2.*xx(2,:).^3]);
Nlift = Nrbf + n;


%% ************************** Collect data ********************************
tic
disp('Starting data collection')
Nsim = 80;
Ntraj = 10000;   

% Random forcing
Ubig =2*rand([Nsim Ntraj]) - 1;
W_d=0.02*(rand([Nsim Ntraj])-0.5);
% Random initial conditions
Xcurrent = 2.5*(rand(n,Ntraj)*2 - 1);
Wd=[];
X = []; Y = []; U = [];
for i = 1:Nsim
    Xnext = f_ud(0,Xcurrent,Ubig(i,:),W_d(i,:));
    X = [X Xcurrent];
    Y = [Y Xnext];
    U = [U Ubig(i,:)];
    Wd = [Wd W_d(i,:)];
    Xcurrent = Xnext;
end
fprintf('Data collection DONE, time = %1.2f s \n', toc);


%% ******************************* Lift ***********************************

disp('Starting LIFTING')
tic
Xlift = liftFun(X);
Ylift = liftFun(Y);
fprintf('Lifting DONE, time = %1.2f s \n', toc);

%% ********************** Build predictor *********************************

disp('Starting REGRESSION')
tic
W = [Ylift ; X];
V = [Xlift; U;Wd];
VVt = V*V';
WVt = W*V';
alpha=diag([500,50,100,50,0,0]); %it is advisable to tune alpha such that A is Schur stable (if possible)
M = WVt * pinv(VVt+alpha'); % Matrix [A B; C 0]
Alift = M(1:Nlift,1:Nlift);
Blift = M(1:Nlift,Nlift+1:end-1);
Dlift = M(1:Nlift,Nlift+2:end);
Clift = M(Nlift+1:end,1:Nlift);

fprintf('Regression done, time = %1.2f s \n', toc);

%% *********************** Predictor comparison ***************************

Tmax = 0.5;
Nsim = Tmax/deltaT;
u_dt = @(i)((-1).^(round(i/5))); % control signal

% Initial condition
x0 = [0.5;0.5];
x_true = x0;

% Lifted initial condition
xlift = liftFun(x0);

% Local linearization predictor at x0
x = sym('x',[2;1]); u = sym('u',[1;1]);w = sym('w',[1;1]);
Ac_x0 = double(subs(jacobian(f_u(0,x,u,w),x),[x;u;w],[x0;0;0]));
Bc_x0 = double(subs(jacobian(f_u(0,x,u,w),u),[x;u;w],[x0;0;0]));
Dc_x0 = double(subs(jacobian(f_u(0,x,u,w),w),[x;u;w],[x0;0;0]));
c_x0 = double(subs(f_u(0,x,u,w),[x;u;w],[x0;0;0])) - Ac_x0*x0 - Bc_x0*0-Dc_x0*0;
ABc = expm([Ac_x0 [Bc_x0 Dc_x0 c_x0] ; zeros(3,5)]*deltaT); % discretize
Ad_x0 = ABc(1:2,1:2); Bd_x0 = ABc(1:2,3); Dc_x0 = ABc(1:2,4);cd_x0 = ABc(1:2,5);
X_loc_x0 = x0;

% Local linearization predictor at 0
x = sym('x',[2;1]); u = sym('u',[1;1]);w = sym('w',[1;1]);
Ac_0 = double(subs(jacobian(f_u(0,x,u,w),x),[x;u;w],[0;0;0;0]));
Bc_0 = double(subs(jacobian(f_u(0,x,u,w),u),[x;u;w],[0;0;0;0]));
c_0 = double(subs(f_u(0,x,u,w),[x;u;w],[0;0;0;0])) - Ac_0*[0;0] - Bc_0*0-Dc_x0*0;
ABc = expm([Ac_0 [Bc_0 Dc_x0 c_0] ; zeros(3,5)]*deltaT); % discretize
Ad_0 = ABc(1:2,1:2); Bd_0 = ABc(1:2,3);Dc_0 = ABc(1:2,4); cd_0 = ABc(1:2,4); 
X_loc_0 = x0;


% Simulate
for i = 0:Nsim-1
    % Koopman predictor
    xlift = [xlift, Alift*xlift(:,end) + Blift*u_dt(i)+Dlift*0.02*(rand(1)-0.5)]; % Lifted dynamics
    
    % True dynamics
    x_true = [x_true, f_ud(0,x_true(:,end),u_dt(i),0.02*(rand(1)-0.5)) ];
    
    % Local linearization predictor at x0
    X_loc_x0 = [X_loc_x0, Ad_x0*X_loc_x0(:,end) + Bd_x0*u_dt(i) + cd_x0+Dc_x0*0.02*(rand(1)-0.5)];
    
    % Local linearization predictor at 0
    X_loc_0 = [X_loc_0, Ad_0*X_loc_0(:,end) + Bd_0*u_dt(i) + c_0+Dc_0*0.02*(rand(1)-0.5)];
    
end
x_koop = Clift * xlift; % Koopman predictions



%% ****************************  Plots  ***********************************

lw = 4;

figure
plot([0:Nsim-1]*deltaT,u_dt(0:Nsim-1),'linewidth',lw); hold on
title('Control input $u$', 'interpreter','latex'); xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)

figure
plot([0:Nsim]*deltaT,x_true(2,:),'linewidth',lw); hold on
plot([0:Nsim]*deltaT,x_koop(2,:), '--r','linewidth',lw)
plot([0:Nsim]*deltaT,X_loc_x0(2,:), '--g','linewidth',lw-1)
plot([0:Nsim]*deltaT,X_loc_0(2,:), '--k','linewidth',lw-1)
axis([0 Tmax min(x_koop(2,:))-0.15 max(x_koop(2,:))+0.15])
title('Predictor comparison - $x_2$','interpreter','latex'); xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)
LEG = legend('True','Koopman','Local at $x_0$','Local at 0','location','southwest');
set(LEG,'interpreter','latex')

figure
plot([0:Nsim]*deltaT,x_true(1,:),'linewidth',lw); hold on
plot([0:Nsim]*deltaT,x_koop(1,:), '--r','linewidth',lw)
plot([0:Nsim]*deltaT,X_loc_x0(1,:), '--g','linewidth',lw-1)
plot([0:Nsim]*deltaT,X_loc_0(1,:), '--k','linewidth',lw-1)
axis([0 Tmax min(x_koop(1,:))-0.1 max(x_koop(1,:))+0.1])
title('Predictor comparison - $x_1$','interpreter','latex'); xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)
LEG = legend('True','Koopman','Local at $x_0$','Local at 0','location','southwest');
set(LEG,'interpreter','latex')



