
function [sr,ap] = consys_cent(k,sr,ap,gp)

mean_v  = sr.v(:,k); % measurement is the mean
mean_p  = sr.x(gp.Mp,k);
mean_f  = sr.x(gp.Mf,k);
mean_e  = sr.x(end,k);

yalmip('clear');

cons  = [];
cost  = 0;

% define decision variables for the windfarm
U  = sdpvar(gp.Nu*gp.Nh,1);

for kk=1:gp.Ns % loop over the samples
    
    sr.v(:,k)     = mean_v + ap.sigma_v*randn(gp.Na,1);
    sr.x(gp.Mp,k) = ap.cp(:,k).*sr.v(:,k).^3.*sr.u(:,k);
    sr.x(gp.Mf,k) = -ap.cf(:,k).*sr.v(:,k).^2.*sr.u(:,k);
    sr.x(end,k)   = gp.Pnref(k)-sum(sr.x(gp.Mp,k));
    xinit(:,kk)   = sr.x(:,k);
    
    
    % build wind farm model
    [gp,ap] = wfmodel(k,sr,ap,gp);
    
    % build contraints set
    gp = constset(k,ap,gp);
    
    % build matrices horizon
    gp = matsys(gp);
    
    
    %X    = gp.A*(gp.T(1:gp.n,:)*xinit) + gp.B*U + gp.Br*gp.Pnref(k:k+gp.Nh-1);%
    X    = gp.A*xinit(:,kk) + gp.B*U + gp.Br*gp.Pnref(k:k+gp.Nh-1);%
    Yobs = gp.C*X;
    
    Fobs       = Yobs(gp.MF);
    Pobs       = Yobs(gp.MP);
    Eobs(:,kk) = Yobs(gp.ME);
    
    cons = [cons, gp.ylb(gp.MP) <= Pobs <= gp.yub(gp.MP)];
    
end

q = 1/gp.Ns*ones(gp.Ns,1);
for kk=1:gp.Ns
    cost = cost + Eobs(:,kk)'*gp.Q*q(kk)*Eobs(:,kk);
end

sr.v(:,k)       = mean_v; % store the means as measurement again
sr.x(gp.Mp,k)   = mean_p;
sr.x(gp.Mf,k)   = mean_f;
sr.x(end,k)     = mean_e;

cons = [cons, gp.ulb <= U <= gp.uub];

if k == 1
    dUt = [ U(1:gp.Na)-sr.u(:,k) ; U(gp.Na+1:end)-U(1:end-gp.Na)];
    cons = [cons, -gp.duc <= dUt <= gp.duc];
else
    dUt = [ U(1:gp.Na)-sr.u(:,k-1) ; U(gp.Na+1:end)-U(1:end-gp.Na)];
    cons = [cons, -gp.duc <= dUt <= gp.duc];
end

cost = cost + dUt'*gp.R*dUt ;

%% finite horizon optimization problem

ops = sdpsettings('solver','cplex','verbose',0,'cachesolvers',1);

optimize(cons,cost,ops)


%% assign the decision variables
Yopt          = value(Yobs);
Uopt_hat      = Yopt(gp.MU);
temp          = reshape(Uopt_hat,[gp.Na,gp.Nh]);
for i=1:gp.Nh
    sr.U(i,k,:)   = temp(:,i); % full horizon optimal action
end

sr.u(:,k)   = sr.U(1,k,:); % 1st step optimal action

% check if no numerical error occured, otherwise take previous
if k>1
    temp         = sr.u(:,k)>ap.uM;
    sr.error     = sr.error+temp;
    sr.u(temp,k) = sr.u(temp,k-1);
end

% get unfiltered control signal
Uopt           = value(U);
temp           = reshape(Uopt,[gp.Na,gp.Nh]);
sr.u_uf(:,k)   = temp(:,1); % full horizon optimal action

end




