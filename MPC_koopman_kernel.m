
%
% Formulation of KMPC
%% REQUIRE: mpt toolbox, yalmip toolbox
%close all
clear all
 load van_koopman_paper_test4-2.mat
 scalf=1.2;
 w=scalf*E_lift_max;
 v=scalf*E_x_max;
load rpiset_kernel_0223-1.mat;
load Xmpi_terminal_kernel_0223.mat;
nx=4;nu=1;
Ad=Alift;
Bd=Blift;
Cd=Clift;

tau=0.01;

alpha=1.1;
Av=[eye(2);-eye(2)];
bv=[v;v];
V= Polyhedron(Av,bv);
V=V.V;
V= Polyhedron(V)
Ax=[eye(2);-eye(2)];
bx=[2.5-v(1), 2.5-v(1),2.5+v(1), 2.5+v(1)]';
Xc = Polyhedron(Ax*Clift,bx);
Uc_vertex = [10; -10];
Uc = Polyhedron(Uc_vertex);
%% Simulate the closed loop
ydd0=[1.5;-1.5];
yd0=ydd0;
xd0 =liftFun(ydd0);
xd=xd0;
N_sim=400;
N=10;
T_CPU=[];
Q=diag([1,1,0.1*ones(1,nx-2)]);
R=0.1*eye(nu);
[K, P] = dlqr(Ad,Bd,diag([1,1,0.1*ones(1,nx-2)]),0.1);
K=-K;
Zc=Cd*Z;
Xc=Xc-alpha*Z;
Uc=Uc-alpha*K*Z;
Xc=Xc;
Uc=Uc;
uc = sdpvar(repmat(nu,1,N),repmat(1,1,N));
xc = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
x_ini=sdpvar(repmat(nx,1,1),repmat(1,1,1));
constraints=[];
objective = 0;
  for i = 1:N
    objective = objective + (xc{i})'*Q*(xc{i})+uc{i}'*R*uc{i};
    constraints = [constraints, xc{i+1} == Ad*xc{i}+Bd*uc{i}];
     constraints = [constraints, Uc.A*uc{i}<= Uc.b];%input constraint
    constraints = [constraints, Xc.A*xc{i}<=Xc.b];
  end
constraints = [constraints, Z.A*(x_ini-xc{1})<Z.b];
objective = objective + xc{N+1}'*P*xc{N+1};
constraints = [constraints,Xmpi_robust.A*xc{N+1}<=Xmpi_robust.b];
options =sdpsettings('solver','quadprog');%solver
controller = optimizer(constraints, objective,options,x_ini,{uc{1},xc{1}});
J=0;
Jx=0;
Ju=0;
trial=1;
U=zeros(trial,N_sim);
U_h=zeros(trial,N_sim);
X=zeros(4*trial,N_sim+1);
Y=zeros(2*trial,N_sim+1);
Y_h=zeros(2*trial,N_sim+1);
J_mpc_col=zeros(1,N_sim);
 for iter=1:trial

yd=yd0;
for k=1:N_sim  
ydd=yd;
 xc0=liftFun(ydd);
 xc0_hat=liftFun(Cd*xd);
tic;
UC= controller{xc0};%online computation
t1=toc;
T_CPU=[T_CPU,t1];
uc =real(UC{1});
xd=real(UC{2});
  ucc=uc+K*(xc0-xd);
 Y(iter*2-1:iter*2,k)=yd;
 X(iter*nx-(nx-1):iter*nx,k)=xd;
  Y_h(iter*2-1:iter*2,k)=Cd*xd;
    U(iter,k)=ucc;
    U_h(iter,k)=uc;
 xd =Ad*xd+Bd*uc;%system updating
  if rem(iter,10)-1==0
     ww=0.4*(rand(1)-0.5);
 end
 yd=f_ud(0,yd,ucc,0*ww+0*0.4*sin(10*k*0.01*pi));%0.4*sin(10*k*0.01*pi)0*0.4*sin(10*k*0.01*pi) 0.8*(rand(1)-0.5)

Jx=Jx+yd'*yd;
Ju=Ju+ucc'*ucc;
J=J+yd'*yd+0.1*ucc'*ucc;
J_mpc_col(k)=yd'*yd+0.1*ucc'*ucc;
    fprintf('Simulating... %d/%d\n',k,N_sim);
end

 end
 %%
 figure
subplot(2,1,1)
plot(1:N_sim,U,'LineWidth',1.5)
axis([0 N_sim -2 6]);
title('r-KMPC','Interpreter','latex')
ylabel('Real control $u$','Interpreter','latex');
 set(gca,'Fontname', 'Times New Roman','FontSize',14);
  grid
subplot(2,1,2)
plot(1:N_sim+1,Y(1,:),'k-',1:N_sim+1,Y(2,:),'b--','LineWidth',1.5)
xlabel('Simulation steps')
axis([0 N_sim -1.8 1.8]);
ylabel('Real state $x$','Interpreter','latex');
 set(gca,'Fontname', 'Times New Roman','FontSize',14);
legend('$x_1$','$x_2$','Interpreter','latex')
 set(gca,'Fontname', 'Times New Roman','FontSize',14);
 grid
%%
figure
subplot(2,2,1)
plot(1:N_sim,U_h,'LineWidth',1.5)
title('Results of r-KMPC using RBF under sinusoidal noise')
axis([0 N_sim -2 6]);
ylabel('Nominal control $\hat u$','Interpreter','latex');
 set(gca,'Fontname', 'Times New Roman','FontSize',14);
grid
subplot(2,2,2)
plot(1:N_sim,U,'LineWidth',1.5)
axis([0 N_sim -2 6]);
ylabel('Real control $u$','Interpreter','latex');
 set(gca,'Fontname', 'Times New Roman','FontSize',14);
  grid
subplot(2,2,3)
plot(1:N_sim+1,X(1,:),'k-',1:N_sim+1,X(2,:),'b--',1:N_sim+1,X(3,:),'r-.',1:N_sim+1,X(4,:),'c.','LineWidth',1.5)
xlabel('Simulation steps')
ylabel('Nominal state $\hat s$','Interpreter','latex');
 set(gca,'Fontname', 'Times New Roman','FontSize',14);
axis([0 N_sim -2.5 3.5]);
legend('$\hat s_1$','$\hat s_2$','$\hat s_3$','$\hat s_4$','Interpreter','latex')
 set(gca,'Fontname', 'Times New Roman','FontSize',14);
grid
subplot(2,2,4)
plot(1:N_sim+1,Y(1,:),'k-',1:N_sim+1,Y(2,:),'b--','LineWidth',1.5)
xlabel('Simulation steps')
axis([0 N_sim -1.8 1.8]);
ylabel('Real state $x$','Interpreter','latex');
 set(gca,'Fontname', 'Times New Roman','FontSize',14);
legend('$x_1$','$x_2$','Interpreter','latex')
 set(gca,'Fontname', 'Times New Roman','FontSize',14);
 grid
Kb=K*pinv(Cd)*Cd;
figure
subplot(2,2,2)
for i=1:2:N_sim
plot([i;U_h(:,i)]+[0,0,0,0;K]*alpha*Z,'color',[192,192,192]/255)%[192,192,192]/255
hold on
end
h1=plot(1:2:N_sim,U_h(1:2:N_sim),'k*',1:2:N_sim,U(1:2:N_sim),'b*','LineWidth',1.5)
legend(h1,'Nominal control','Real control','Interpreter','tex')
xlabel('Simulation steps')
title('Control $u/\hat u$ of r-KMPC','Interpreter','latex')
ylabel('$u/\hat u$','Interpreter','latex');
 set(gca,'Fontname', 'Times New Roman','FontSize',16);
subplot(2,2,1)
for i=1:2:N_sim
plot(Y_h(:,i)+alpha*Cd*Z+V,'color',[192,192,192]/255)
hold on
end
h1=plot(Y_h(1,1:2:N_sim),Y_h(2,1:2:N_sim),'k+',Y(1,1:2:N_sim),Y(2,1:2:N_sim),'b+','LineWidth',1.5)
legend(h1,'Nominal trajectory','Real trajectory','Interpreter','tex')
xlabel('$x_1/\hat x_1$','Interpreter','latex')
title('Trajectory of $x/\hat x$ of r-KMPC','Interpreter','latex')
ylabel('$x_2/\hat x_2$','Interpreter','latex');
 set(gca,'Fontname', 'Times New Roman','FontSize',16);