function [gp,ap,sr] =  InitializeController(Wp)


%% global parameters

gp = struct(...
    'Na','',...     % number of agents
    'Nx','',...     % the number of the states (imbalance_error) of each agent
    'Nu','',...     % the number of the real decision variables
    'Nh','',...     % prediction horizon
    'Nsim','');     % simulation steps

gp.Nsim       = Wp.N/(2*Wp.h);      % total simulation time (max=1200)
gp.Nh         = 10;                 % # samples in prediction horizon
gp.Ns         = 5;                 % # samples of uncertain variables
gp.Na         = Wp.Nt;              % # wind turbines/agents
gp.Nu         = gp.Na;              % #inputs
gp.Nx         = 3*gp.Na+1;          % #states

gp.MF         = logical(repmat([repmat([1 0 0]',gp.Na,1); 0],gp.Nh,1));
gp.Mf         = logical([repmat([1 0 0]',gp.Na,1); 0]);
gp.MP         = logical(repmat([repmat([0 1 0]',gp.Na,1); 0],gp.Nh,1));
gp.Mp         = logical([repmat([0 1 0]',gp.Na,1); 0]);
gp.MU         = logical(repmat([repmat([0 0 1]',gp.Na,1);0],gp.Nh,1));
gp.ME         = logical(repmat([zeros(gp.Nx-1,1); 1],gp.Nh,1));

gp.Q          = 1e4*eye(gp.Nh);                 % weigth on tracking
gp.S          = 0*eye(gp.Nh*gp.Na);             % weigth on variation of the force
gp.R          = 1e4*eye(gp.Nh*gp.Nu);         % weigth on control signal

gp.duc        = 2e-1;                           % limitation on du/dt
gp.dfc        = Inf;                            % limitation on dF/dt

% wind farm reference
gp.Pnref      = zeros(gp.Nsim+gp.Nh,1); 
load(strcat(Wp.controller,'/libraries/AGC_PJM_RegD_Norm2s'));
gp.AGCdata = AGCdata(:,2); 

gp.Pgreedy            =  1.125081579727473e+07; 
gp.Pnref(1:Wp.N0)     = .7*gp.Pgreedy; 
gp.Pnref(Wp.N0+1:end) = .7*gp.Pgreedy+.3*gp.Pgreedy*gp.AGCdata(1:gp.Nsim+gp.Nh-Wp.N0);

%% agents parameters

ap              = struct;

ap.Rr           = 60;                                            % rotor radius
ap.tau          = 5;                                             % time constant filter CT'
ap.Ts           = 1;
ap.cf           = 0.5*pi*ap.Rr^2*ones(gp.Na,gp.Nsim);
ap.cp           = 0.5*pi*ap.Rr^2*ones(gp.Na,gp.Nsim);
ap.uM           = 2;                                             % maximum CT'
ap.yaw          = zeros(1,gp.Na);
ap.sigma_v      = .3;
rng(1)

[num,den]       = tfdata(c2d(tf(1,[ap.tau 1]),ap.Ts),'v');       % filter on force
ap.sys          = c2d(tf(1,[ap.tau 1]),1);

for kk = 1:gp.Na    
    ap.a{kk}     = kron(eye(3),-den(2));
    ap.bcoef{kk} = num(2)*[-1;1;1];
    ap.c{kk}     = eye(3);                                       
end                                                                                                                              


% constraints definitions
ap.duc        = 5e-2;                                               % limitation on du/dt
ap.dfc        = Inf*ones(gp.Na,gp.Nsim);                            % limitation on dF/dt
ap.MF         = logical(repmat([1 0 0]',gp.Nh,1));
ap.Mf         = logical([repmat([1 0 0]',gp.Na,1); 0]);
ap.MP         = logical(repmat([0 1 0]',gp.Nh,1));
ap.Mp         = logical([repmat([0 1 0]',gp.Na,1); 0]);
ap.MU         = logical(repmat([0 0 1]',gp.Nh,1));
ap.Mu         = logical([repmat([0 0 1]',gp.Na,1); 0]);


%% simulation results

sr                = struct;

sr.x              = zeros(gp.Nx,gp.Nsim);  % state x=[F P]
sr.U              = zeros(gp.Nh,gp.Nsim,gp.Na);
sr.u              = zeros(gp.Na,gp.Nsim);  % control signal CT'
sr.CT             = zeros(gp.Na,gp.Nsim);  % control signal CT
sr.v              = 5*ones(gp.Na,gp.Nsim); % wind velocity at rotor (needs initial guess)
sr.y              = zeros(gp.Nx,gp.Nsim);  % output y=[F P]
sr.e              = zeros(1,gp.Nsim);      % wind farm error
sr.error          = zeros(gp.Na,1);        % #of wrong optimization steps

end
