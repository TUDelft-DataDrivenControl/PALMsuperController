function [yWake,UWake,fWake] = estimateWakeCenter(y,U,nW,y0,plotFig,y_ub,y_lb)
% y  = lateral distance measure, vector
% U  = wind speed along y (at hub height), vector
% nW = number of gaussians/number of wakes to separate, integer
% y0 = initial guess for wake locations, vector with nW entries
% plotFig = boolean to plot fit results (default: false)
% y_ub = upper bound on wake location   (default: max(y)+100)
% y_lb = lower bound on wake location   (default: min(y)-100)

if nargin <= 4
    plotFig = false;
end
if nargin <= 5
    ycenter_bnds = [min(y)-100 max(y)+100]; % Lower and upper bound of y_center
end

% Format U to a vector format for fit function
U_tmp(:,1) = U;
U          = U_tmp;

% Specify initial guesses
U_inf0  = median(U); % Initial guess for freestream wind speed (m/s)
U_min0  = min(U);    % Initial guess for minimum wind speed (m/s)
D_wake0 = 100;       % Initial guess for wake diameter (m)

% Setup (multi-)gaussian function
gaussianEqn = '';
for j = 1:nW
    % Add a gaussian term with corresponding ubs and lbs to the function
    gaussianEqn = [gaussianEqn '+ a' num2str(j) ' *exp(-(x-b' ...
                   num2str(j) ')^2/(2*c' num2str(j) '^2))'];
end
gaussianEqn = [gaussianEqn '+k'];

% Setup initial conditions, upper bounds, and lower bounds for fit
x0 = [ones(1,nW)*(U_min0-U_inf0), y0,                         D_wake0*ones(1,nW), U_inf0];
ub = [zeros(1,nW),                ones(1,nW)*ycenter_bnds(2), ones(1,nW)*Inf,     max(U)];  
lb = [ones(1,nW)*-max(U),         ones(1,nW)*ycenter_bnds(1), zeros(1,nW),        min(U)];

% Fit and report scores
fWake = fit(y,U,gaussianEqn,'Start', x0,'Upper',ub,'Lower',lb);
          
% plot
if plotFig
    figure;
    plot(y,U); hold on;
    plot(y,fWake(y),'--');
    legend('LES','Fit');
    xlabel('y (m)');
    ylabel('U (m/s)');
    title('(Multi-)gaussian fit of for wake deficit');
    grid minor;
end

% collect outputs
yWake = [];
UWake = [];
for j = 1:nW
    yWake = [yWake fWake.(['b' num2str(j)])]; % Wake center location (m)
    UWake = [UWake fWake.k+fWake.(['a' num2str(j)])]; % Wake deficit (m/s)
end
end
