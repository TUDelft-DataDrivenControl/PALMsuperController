function [gp,ap] = wfmodel(k,sr,ap,gp)


% build B matrix of the turbine models for current times step
for i = 1:gp.Na
    ap.SwMes{i} = sr.v(i,k);
    v           = repmat(ap.SwMes{i},gp.Nh,1);
    cf          = ap.cf(i,k);
    cp          = ap.cp(i,k);
    
    ap.b{i} = [ap.bcoef{i}(1)*(v(1)*cosd(ap.yaw(i)))^2*cf ap.bcoef{i}(2)*(v(1)*cosd(ap.yaw(i)))^3*cp ap.bcoef{i}(3)]';


    for j = 1:gp.Nh
    ap.bt(:,i,j) = [ap.bcoef{i}(1)*(v(j)*cosd(ap.yaw(i)))^2*cf ap.bcoef{i}(2)*(v(j)*cosd(ap.yaw(i)))^3*cp ap.bcoef{i}(3)]';
    end
    
end

% build wind farm model
gp.a = blkdiag(ap.a{:});
gp.b = blkdiag(ap.b{:});
for j = 1:gp.Nh
    for i = 1:gp.Na 
    b{i}         = ap.bt(:,i,j);
    end
    gp.bt(:,:,j) = [blkdiag(b{:}); zeros(1,gp.Na)]; % add directly error channel
end
gp.c = blkdiag(ap.c{:});
% add error channel
gp.a  = [gp.a zeros(gp.Nx-1,1)];
gp.a  = [gp.a;-repmat([0 1 0],1,gp.Na) 0];
gp.b  = [gp.b; zeros(1,gp.Na)];
gp.c  = blkdiag(gp.c, 1);
% build B matrix for wind farm reference
gp.br = [zeros(gp.Nx-1,1); 1];

% get minimal realisation
sys                   = ss(gp.a,[gp.b gp.br],gp.c,0,1); % assumes one second sample period!
[sys_min,gp.T]        = minreal(sys,[],false);
[gp.am,gp.bm,gp.cm,~] = ssdata(sys_min); 
gp.n                  = size(gp.am,1);
gp.brm                = gp.bm(:,end);
gp.bm                 = gp.bm(:,1:end-1);

% uncomment following to get minimal
% gp.a  = gp.am; 
% gp.b  = gp.bm;
% gp.br = gp.brm;
% gp.c  = gp.cm;

end

