%% clear workspace

clear

%% load problem

##load('prob_quad_var2_ineq4.mat',...
##'n','P','xc','m','A','b','r','D','e','x0');
load('prob_quad_var100_ineq200.mat',...
'n','P','xc','m','A','b','r','D','e','x0');

f = @(x) sum((x-xc).'*P*(x-xc))/2;
df = @(x) P*(x-xc);

h = cell(1,m);
dh = cell(1,m);
for i=1:m
    h{i} = @(x) A(i,:)*x-b(i);
    dh{i} = @(x) A(i,:).';
end

g = cell(1,r);
dg = cell(1,r);
for i=1:r
    g{i} = @(x) D(i,:)*x-e(i);
    dg{i} = @(x) D(i,:).';
end

prob.f = f;
prob.df = df;
prob.h = h;
prob.dh = dh;
prob.g = g;
prob.dg = dg;
prob.n = n;

%% set parameters

param.c0 = 1e+0;
param.delta = 1+1e-0;
param.gamma = .25;

param.epsilon0 = 1e-0;
param.omega = .8;

param.s0 = 1;
param.alpha = .1;
param.beta = .6;

param.pth = 1e-3;
param.dth = 1e-3;

param.max_out_itr = 30;
param.max_in_itr = 10;
param.min_step = 1e-12;

param.x0 = x0;
param.lambda0 = zeros(numel(prob.h),1);
param.mu0 = zeros(numel(prob.g),1);

%% perform augmented Lagragian

res = augmented_Lagrangian(prob,param);

%% plot result

T =res.T;
x_hist = res.x_hist(:,1:T);
lambda_hist = res.lambda_hist(:,1:T);
mu_hist = res.mu_hist(:,1:T);
c_hist = res.c_hist(1:T);
PRx = cellfun(...
    res.PR,...
    num2cell(x_hist,1),...
    num2cell(mu_hist,1),...
    num2cell(c_hist,1)...
    );
DRx = cellfun(...
    res.DR,...
    num2cell(x_hist,1),...
    num2cell(lambda_hist,1),...
    num2cell(mu_hist,1),...
    num2cell(c_hist,1)...
    );

in_dur = res.in_dur(1:T);

epsilon_hist = res.epsilon_hist(1:T);

inv_c_hist = res.c_hist(1:T).^(-1);

figure
semilogy(in_dur,PRx,'-');
hold on;
semilogy(in_dur,DRx,'--');
plot(in_dur,epsilon_hist,':','LineWidth',0.1);
plot(in_dur,inv_c_hist,'-.','LineWidth',0.1);

legend({'PR','DR','epsilon','inv\_c'});
xlabel('time(s)');
set(gca,'FontSize',12);

%% save result

save result_AL.mat
