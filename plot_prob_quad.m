%% load problem and result

load('prob_quad_var2_ineq4.mat','D','e'); 
load('result_AL.mat','prob','x_hist','x0');

%% setup figure

x = linspace(-40,100);
y = linspace(-40,100);
[X,Y] = meshgrid(x,y);

figure
xlim([-40,100]);
ylim([-40,100]);
hold on;

%% calculate extreme points

expnt = [];
for i=1:size(D,1)
    for j=(i+1):size(D,1)
        if abs(det(D([i,j],:))) > 1e-3
            z = D([i,j],:)\e([i,j]);
            is_fes = true;
            for k=1:size(D,1)
                if D(k,:)*z-e(k) > 1e-3
                    is_fes = false;
                    break;
                end
            end
            if is_fes
                expnt = [expnt,z];
            end
        end
    end
end

%% plot feasible region

K = convhull(expnt(1,:),expnt(2,:));
plot(expnt(1,K),expnt(2,K),'k-');

%% plot object function

Z = arrayfun(@(x,y) prob.f([x;y]),X,Y);
contour(X,Y,Z,'LineColor','black');

%% plot trajectory

plot([x0(1),x_hist(1,:)],[x0(2),x_hist(2,:)],'k--+');
