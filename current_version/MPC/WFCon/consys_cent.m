
function [sr,ap] = consys_cent(k,sr,ap,gp)

yalmip('clear');

cons  = [];

xinit = sr.x(:,k);

% define decision variables for the windfarm
U  = sdpvar(gp.Nu*gp.Nh,1);

% build wind farm model
[gp,ap] = wfmodel(k,sr,ap,gp);

% build contraints set
gp = constset(k,ap,gp);

% build matrices horizon
gp = matsys(gp);


X    = gp.A*xinit + gp.B*U + gp.Br*gp.Pnref(k:k+gp.Nh-1);%
Yobs = gp.C*X;

Fobs = Yobs(gp.MF);
Pobs = Yobs(gp.MP);
Eobs = Yobs(gp.ME);


cons = [cons, gp.ulb <= U <= gp.uub];

if k == 1
    dUt = [ U(1:gp.Na)-sr.u(:,k) ; U(gp.Na+1:end)-U(1:end-gp.Na)];
    cons = [cons, -gp.duc <= dUt <= gp.duc];
else
    dUt = [ U(1:gp.Na)-sr.u(:,k-1) ; U(gp.Na+1:end)-U(1:end-gp.Na)];
    cons = [cons, -gp.duc <= dUt <= gp.duc];
end

cons = [cons, gp.ylb(gp.MP) <= Pobs <= gp.yub(gp.MP)];
cons = [cons, gp.ylb(gp.MF) <= Fobs <= gp.yub(gp.MF)];

if k == 1
    dFt  = [ Fobs(1:gp.Na)-sr.y(ap.Mf,k) ; Fobs(gp.Na+1:end)-Fobs(1:end-gp.Na)];
else
    dFt = [ Fobs(1:gp.Na)-sr.y(ap.Mf,k-1) ; Fobs(gp.Na+1:end)-Fobs(1:end-gp.Na)];
end

cost = Eobs'*gp.Q*Eobs + U'*gp.R*U + dFt'*gp.S*dFt ;


%% finite horizon optimization problem

ops = sdpsettings('solver','cplex','verbose',0,'cachesolvers',1);

optimize(cons,cost,ops)


%% assign the decision variables
Xopt          = value(X);
Uopt_hat      = Xopt(gp.MU);
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


end




