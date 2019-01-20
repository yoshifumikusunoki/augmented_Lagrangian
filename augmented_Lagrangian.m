function res = augmented_Lagrangian(prob,param)

f = prob.f;
df = prob.df;
h = prob.h;
dh = prob.dh;
g = prob.g;
dg = prob.dg;

n = prob.n;
m = numel(h);
r = numel(g);

c = param.c0;
gamma = param.gamma;
delta = param.delta;

epsilon = param.epsilon0;
omega = param.omega;

s0 = param.s0;
alpha = param.alpha;
beta = param.beta;

pth = param.pth;
dth = param.dth;

max_out_itr = param.max_out_itr;
max_in_itr = param.max_in_itr;
min_step = param.min_step;

AL = @(x,lambda,mu,c) f(x) ...
    + cellfun(@(F)F(x),h)*lambda ...
    + (c/2)*sum(cellfun(@(F)F(x),h).^2) ...
    + (1/(2*c))*sum(max(0,mu.'+c*cellfun(@(F)F(x),g)).^2-mu.'.^2);

dAL = @(x,lambda,mu,c) df(x)...
    + cell2mat(cellfun(@(F)F(x),dh,'UniformOutput',false))*(lambda+c*cellfun(@(F)F(x),h).') ...
    + cell2mat(cellfun(@(F)F(x),dg,'UniformOutput',false))*max(0,mu+c*cellfun(@(F)F(x),g).');

PR = @(x,mu,c) norm([cellfun(@(F)F(x),h),max(-mu.'/c,cellfun(@(F)F(x),g))]);
DR = @(x,lambda,mu,c) norm(dAL(x,lambda,mu,c));

x_hist = zeros(n,max_out_itr);
lambda_hist = zeros(m,max_out_itr);
mu_hist = zeros(r,max_out_itr);
c_hist = zeros(1,max_out_itr);
epsilon_hist = zeros(1,max_out_itr);
in_dur = zeros(1,max_out_itr);

x = param.x0;
lambda = param.lambda0;
mu = param.mu0;

k = 1;
PRxmc = PR(x,mu,c);
tic
while 1
    
    if k > max_out_itr
        fprintf('Terminated before both PR and DR reach thresholds.\n');
        k = max_out_itr;
        break;
    end
 
    lambda_hist(:,k) = lambda;
    mu_hist(:,k) = mu;
    c_hist(:,k) = c;
    epsilon_hist(:,k) = epsilon;
    
    l = 1;
    while 1
        
        if l > max_in_itr
            fprintf('Fail to minimize augmented Lagrangian function.\n');
            break;            
        end
        
        dALx = dAL(x,lambda,mu,c);
        
        if norm(dALx) <= epsilon
            break;
        end
        
        ALx = AL(x,lambda,mu,c);
        q = 0;
        while 1
            s = beta^q*s0;
            if s < min_step
                break;
            end
            ALxx = AL(x-s*dALx,lambda,mu,c);
            if ALx - ALxx >= alpha*s*(dALx.'*dALx)
                break
            end
            q = q + 1;
        end
        if s < min_step
            fprintf('Line search is failed.\n');
            break;
        end
        
        x = x - s*dALx;
        l = l+1;
    end
    in_dur(k) = toc;
    
    x_hist(:,k) = x;
    
    pre_PRxmc = PRxmc;
    PRxmc = PR(x,mu,c);
    DRxlm = DR(x,lambda,mu,c);
    
    fprintf('%3d, %4.2fs, (PR,DR) = (%e,%e)\n',k,in_dur(k),PRxmc,DRxlm);
    
    if PRxmc <= pth && DRxlm <= dth
        break;
    end
       
    lambda = lambda + c*cellfun(@(F)F(x),h).';
    mu = max(0,mu.'+c*cellfun(@(F)F(x),g)).';
    epsilon = max(.1*dth,omega*epsilon);
    
    if PRxmc > gamma*pre_PRxmc
        c = min(1e+9,delta*c);
    end
    
    k = k+1;
    
end

res.prob = prob;
res.param = param;
res.AL = AL;
res.dAL = dAL;
res.PR = PR;
res.DR = DR;
res.x_hist = x_hist;
res.lambda_hist = lambda_hist;
res.mu_hist = mu_hist;
res.c_hist = c_hist;
res.epsilon_hist = epsilon_hist;
res.in_dur = in_dur;
res.T = k;

end


