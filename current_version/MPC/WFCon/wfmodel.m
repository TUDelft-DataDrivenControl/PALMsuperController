function [gp,ap] = wfmodel(k,sr,ap,gp)


% build B matrix of the turbine models for current times step
for i = 1:gp.Na
    ap.SwMes{i} = sr.v(i,k);
    v           = ap.SwMes{i};
    cf          = ap.cf(i,k);
    cp          = ap.cp(i,k);
    
    ap.b{i} = [ap.bcoef{i}(1)*(v*cosd(ap.yaw(i)))^2*cf ap.bcoef{i}(2)*(v*cosd(ap.yaw(i)))^3*cp ap.bcoef{i}(3)]';
end


% build wind farm model
gp.a = blkdiag(ap.a{:});
gp.b = blkdiag(ap.b{:});
gp.c = blkdiag(ap.c{:});
% add error channel
gp.a  = [gp.a zeros(gp.Nx-1,1)];
gp.a  = [gp.a;-repmat([0 1 0],1,gp.Na) 0];
gp.b  = [gp.b; zeros(1,gp.Na)];
gp.c  = blkdiag(gp.c, 1);
% build B matrix for wind farm reference
gp.br = [zeros(gp.Nx-1,1); 1];


end

