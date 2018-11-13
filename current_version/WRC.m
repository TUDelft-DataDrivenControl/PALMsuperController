clear;clc;

load(strcat('WRC/libraries/mesh.mat'))

Wp.name          = '2turb';                  % farm
Wp.controller    = 'WRC';

addpath(genpath(strcat(Wp.controller,'/WFCon')));
addpath(genpath(strcat(Wp.controller,'/libraries')));
addpath(genpath(strcat('USER_CODE/6turb/matlab/mexcdf')));

Wp.filename      = strcat(Wp.name,'_m01.nc');
Wp.Nt            = 2;                        % #turbines
Wp.nz            = 5;                        % z-index of hubs (this should coincide with z-mask defined in _p3d and in the turbine model)

Wp.N             = 2*floor(2000/2);          % simulation time ( real simulation time is N-h )
Wp.h             = 1.0;                      % PALM sample period (measure and apply control signals every 2h)
Wp.k             = 0;                        % time index
Wp.N0            = 25;                       % first N0 samples of the simulation constant reference %round(.05*Wp.N/(2*Wp.h));
WrongMeasurement = 0;                        % amount of failed measurements
MaxTime          = 90*ones(Wp.N/(2*Wp.h),1); % maximum waiting time for measurement
MaxTime(1)       = 6*3600;

Wp.Ny    = 2;                                % #measurements [u,v]
Wp.ym    = zeros(Wp.N/(2*Wp.h),Wp.Ny,size(Wp.x,1),size(Wp.y,1)); % measurements [u,v]
Wp.time  = zeros(Wp.N/(2*Wp.h),1);           % time at which measurement is taken
Wp.phi   = zeros(Wp.N/(2*Wp.h),Wp.Nt);       % control signal

% read and write correct information from PALM input files
% set end_time in _p3df files equal to N1
unix(char(strcat('sed -i "/&d3par/c\&d3par end_time =',{' '},num2str(Wp.N,'%.1f'),'" JOBS/',Wp.name,'/INPUT/',Wp.name,'_p3df')));
% set initializing_actions = 'cyclic_fill' in _pd3f file
unix(char(strcat('sed -i "/initializing_actions/c\initializing_actions =','''','cyclic_fill','''',',','" JOBS/',Wp.name,'/INPUT/',Wp.name,'_p3df')));

% set sample period _p3d file
unix(char(strcat('sed -i "/dt =/c\         dt =',{' '},num2str(Wp.h,'%.1f'),',','" JOBS/',Wp.name,'/INPUT/',Wp.name,'_p3df')));
unix(char(strcat('sed -i "/dt_run_control =/c\         dt_run_control =',{' '},num2str(Wp.h,'%.1f'),',','" JOBS/',Wp.name,'/INPUT/',Wp.name,'_p3df')));
unix(char(strcat('sed -i "/dt_domask =/c\         dt_domask =',{' '},num2str(Wp.h,'%.1f'),',/','" JOBS/',Wp.name,'/INPUT/',Wp.name,'_p3df')));



% get turbine information from _p3d file
fid   = fopen(char(strcat( 'JOBS/',Wp.name,'/INPUT/',Wp.name,'_p3df' )));
tline = fgetl(fid);
while ischar(tline)
    tline = fgetl(fid);
    if strncmpi('rcx',tline,3)                   % the number 3 could be generalised (these digits should coincide with _p3df file)
        Wp.Crx       = str2num(tline(9:end));    % the number 9 could be generalised
    elseif strncmpi('rcy',tline,3)
        Wp.Cry       = str2num(tline(9:end));
    elseif strncmpi('rcz',tline,3)
        Wp.hubHeight = str2num(tline(9:end));
    elseif strncmpi('rr',tline,2)                % the number 2 could be generalised
        Wp.Dr        = str2num(tline(9:end));
    end
end
fclose(fid);
Wp.hubHeight = Wp.hubHeight(1);
Wp.Dr        = 2*Wp.Dr(1);
Wp.Crx       = Wp.Crx-Wp.xu(1);
Wp.Cry       = Wp.Cry-Wp.yv(1);
Wp.x         = Wp.x-Wp.x(1);
Wp.xu        = Wp.xu-Wp.xu(1);
Wp.y         = Wp.y-Wp.y(1);
Wp.yv        = Wp.yv-Wp.yv(1);
Wp.dx        = Wp.x(2)-Wp.x(1);

for ll = 1:length(Wp.Crx)
    [~,Wp.xline(ll,1)]  = min(abs(Wp.x-Wp.Crx(ll)));
    [~, L_prim ]        = min(abs(Wp.y- (Wp.Cry(ll)-Wp.Dr/2)));
    [~, R_prim ]        = min(abs(Wp.y- (Wp.Cry(ll)+Wp.Dr/2)));
    Wp.yline{ll}        = L_prim:1: R_prim;
end
clear L_prim R_prim

% lidar
Wp.nDr   = 3;                 % #rotor diameters downwind of lidar scan                    
xi       = Wp.nDr*Wp.Dr;      % lidar coordinate system
yi       = Wp.y'-Wp.Cry(1);
f        = sqrt(xi.^2+yi.^2); % points at which lidar measurement is taken

Parameter.Lidar.YawAngle    = 0;
Parameter.Lidar.RollAngle   = 0;
Parameter.Lidar.PitchAngle  = 0;
Parameter.Lidar.PositionLinI = [Wp.Crx,Wp.Cry,Wp.hubHeight];
Parameter.Wind.Elevation    = 0;
Parameter.Wind.Azimuth      = 0;
Parameter.Wind.PositionIinW = [0,0,0];
% no weighting function
Parameter.Lidar.a           = [0];
Parameter.Lidar.f_L_d       = [1];
Trajectory                  = load(strcat(Wp.controller,'/libraries/PalmWakeTracking10x10_1_0d5_3.mat'));

% design controller
Parameter.Time.dt                                       = Wp.h;
Parameter.Wind.URef                                     = 8; 
Parameter.Turbine.RotorDiameter                         = Wp.Dr;
Parameter.LidarWakeCenterEstimationPosition             = Wp.nDr*Wp.Dr;
Parameter.GainScheduling                                = 0;
Parameter.ControllerName                                = 'HinfController';
Parameter.WakeActuator.MaxYawRate                       = deg2rad(1);
Parameter.PadeOrder                                     = 10;
Parameter.PlotHinfController                            = false;
load(strcat(Wp.controller,'/libraries/PALMIdentificationWFSim_URef_8_StepID_7_MI_1_results.mat'),'logsout')
Parameter.Plant                                         = tf(logsout.ModelIdentification.sys)*180/pi;
Parameter.LoadExternalPlants                            = true;

HinfController = HinfThesisController('PALMHinfSecondTry',Parameter);

% clean up
unix(char(strcat('rm -f ',{' '},Wp.controller,'/Data/*')));
unix(strcat('rm -f JOBS/',Wp.name,'/MONITORING/*'));
unix(strcat('rm -f USER_CODE/',Wp.name,'/matlab/* 2> /dev/null'));
unix(strcat('rm -f USER_CODE/',Wp.name,'/matlab/.nfs*'));
unix(strcat('rm -f USER_CODE/',Wp.name,'/matlab/previous_y/*'));
unix(strcat('rm -f JOBS/',Wp.name,'/OUTPUT/*'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start PALM
[status,cmdout] = ...
    unix(char(strcat('mrun -d',{' '},Wp.name,' -h lchpc06 -K parallel -X 16 -T 8 -t 360000 -m 64000 -v -b -q dcsc -r "d3f maf"')));

%% Start the loop
while Wp.k<Wp.N/(2*Wp.h)
    
    tic
    Wp.k = Wp.k + 1;
    
    % 1) run controller
    if Wp.k>=1 && Wp.k<=Wp.N0
        
        deltay              = -20;                  % static wake deflection at 3D
        Wp.yWakeRef(Wp.k,:) = Wp.Cry(1)+deltay;
        Wp.phi(Wp.k,:)      = 0;
        
    else
        
        Wp.yWakeRef(Wp.k,:)  = Wp.Cry(1)+deltay;
        HinfController.DesiredWakeCenterPosition = deltay;
        Wp.phi(Wp.k,:)       = GetDemandedYawAngle(HinfController,squeeze(Wp.yWake(Wp.k-1,:,:))-Wp.yWakeRef(Wp.k-1,:));
        Wp.phi(Wp.k,2)       = 0;
    end
    
    Wp.phi(Wp.k,:) = min(max(Wp.phi(Wp.k,:),-35*pi/180),35*pi/180);  % saturate
   
    
    % 2) write control signals to a text file ( are applied at k=0,2h,...,(N-1)*h and hold in between )   
    for ll=1:Wp.Nt
        dlmwrite(strcat('USER_CODE/',Wp.name,'/matlab/u0',num2str(ll),'.txt'),Wp.phi(Wp.k,ll),'precision',8)
    end
    dlmwrite(strcat('USER_CODE/',Wp.name,'/matlab/uflag.txt'),1) % PALM waits for this flag
    fclose('all');
    
    % 3) wait until measurement from PALM is written to file      
    temp = zeros(1,2);
    while prod(temp)==0
        if toc > MaxTime(Wp.k)
            for ll=1:Wp.Nt
                unix(char(strcat('mv USER_CODE/',Wp.name,'/matlab/previous_y/uk.txt USER_CODE/',Wp.name,'/matlab/y0/uk.txt')));
                unix(char(strcat('mv USER_CODE/',Wp.name,'/matlab/previous_y/vk.txt USER_CODE/',Wp.name,'/matlab/y0/vk.txt')));
            end
            WrongMeasurement = WrongMeasurement+1;
            pause(.3)
            break;
        end
        pause(.3)
        temp(1) = exist(strcat('USER_CODE/',Wp.name,'/matlab/uk.txt'),'file');
        temp(2) = exist(strcat('USER_CODE/',Wp.name,'/matlab/vk.txt'),'file');
    end;   

    % 4) read measurement from PALM ( at k=h,3h,...,N*h )    
    Wp.ym(Wp.k,1,:,:) = dlmread( strcat('USER_CODE/',Wp.name,'/matlab/uk.txt') );    % u
    Wp.ym(Wp.k,2,:,:) = dlmread( strcat('USER_CODE/',Wp.name,'/matlab/vk.txt') );    % v
    unix(strcat('sed -i "s/\*/0/g" USER_CODE/',Wp.name,'/matlab/time.txt'));
    Wp.time(Wp.k)     = dlmread( strcat('USER_CODE/',Wp.name,'/matlab/time.txt') );  % t
    
    % remove measurement files
    unix(char(strcat('rm -f USER_CODE/',Wp.name,'/matlab/uk.txt')));
    unix(char(strcat('rm -f USER_CODE/',Wp.name,'/matlab/vk.txt')));
    
    
    % 5) estimate wake position
    nW              = 1;    % #wakes on the x-coord at which lidar measurement is taken
        
    uk              = squeeze(Wp.ym(Wp.k,1,:,:));
    vk              = squeeze(Wp.ym(Wp.k,2,:,:));
    
    Wp.vlos(Wp.k,:) = 1./f .* (xi.*uk(round( (Wp.Crx(1)+Wp.nDr*Wp.Dr) /Wp.dx ) , :) ...
                                + yi.*vk(round( (Wp.Crx(1)+Wp.nDr*Wp.Dr) /Wp.dx ) , :));   
    
    for ll=1:1
        [Wp.yWake(Wp.k,ll,:),Wp.UWake(Wp.k,ll,:),Wp.fWake] = estimateWakeCenter(Wp.y,Wp.vlos(Wp.k,:),nW,Wp.Cry(1),0);
        Wp.fitWake(Wp.k,ll,:)                              = Wp.fWake(Wp.y)';
    end
     
    disp(['Simulated t(' num2str(Wp.k) ')  = ' num2str(Wp.time(Wp.k),3) ...
        ' s.       CPU: ' num2str(toc,3) ' s.']);
    
    % backup after every 50 samples
    if mod(Wp.k,50)==0
        save(char(strcat(Wp.controller,'/Data/series_wt_',num2str(Wp.k))));
    end
end

tic
% fetch measurements from PALM simulation
while (exist(strcat('JOBS/',Wp.name,'/OUTPUT/',Wp.name,'_m01.nc'),'file')...
        *exist(strcat('JOBS/',Wp.name,'/MONITORING/',Wp.name,'_turbine_parameters0',num2str(Wp.Nt)),'file'))==0
    pause(.3)
    if toc > MaxTime(Wp.k)
        dlmwrite(strcat('USER_CODE/',Wp.name,'/matlab/uflag.txt'),1) % PALM waits for this flag
        fclose('all');
        tic
    end
end;
unix(strcat('bash USER_CODE/',Wp.name,'/matlab/mexcdf/convert_measurement.bash'));

%M = [Time   RSpeed  GSpeed  GenTorque  AeroTorque  Pitch Power(Gen)  Power(Rot)  RotThrust  WDirection  YawOrient]
for ll=1:Wp.Nt
    
    unix(strcat('sed -i "s/\*/0/g" JOBS/',Wp.name,'/MONITORING/',Wp.name,'_turbine_parameters0',num2str(ll),'.txt'));
    
    M = dlmread(strcat('JOBS/',Wp.name,'/MONITORING/',Wp.name,'_turbine_parameters0',num2str(ll),'.txt'),'',1,0);
    
    sr.Power(ll,:) = M(:,7)';
    sr.phi(ll,:)   = M(:,11)';
    sr.time(ll,:)  = M(:,1)';
end

sr.tflow = double(nc_varget(strcat('JOBS/',Wp.name,'/OUTPUT/',Wp.name,'/OUTPUT/',Wp.filename),'time'));
temp     = double(nc_varget(strcat('JOBS/',Wp.name,'/OUTPUT/',Wp.name,'/OUTPUT/',Wp.filename),'u'));
sr.u     = squeeze(temp(:,Wp.nz,:,:));
temp     = double(nc_varget(strcat('JOBS/',Wp.name,'/OUTPUT/',Wp.name,'/OUTPUT/',Wp.filename),'v'));
sr.v     = squeeze(temp(:,Wp.nz,:,:));
temp     = double(nc_varget(strcat('JOBS/',Wp.name,'/OUTPUT/',Wp.name,'/OUTPUT/',Wp.filename),'w'));
sr.w     = squeeze(temp(:,Wp.nz,:,:));
sr.x     = double(nc_varget(strcat('JOBS/',Wp.name,'/OUTPUT/',Wp.name,'/OUTPUT/',Wp.filename),'x'));
sr.xu    = double(nc_varget(strcat('JOBS/',Wp.name,'/OUTPUT/',Wp.name,'/OUTPUT/',Wp.filename),'xu'));
sr.y     = double(nc_varget(strcat('JOBS/',Wp.name,'/OUTPUT/',Wp.name,'/OUTPUT/',Wp.filename),'y'));
sr.yv    = double(nc_varget(strcat('JOBS/',Wp.name,'/OUTPUT/',Wp.name,'/OUTPUT/',Wp.filename),'yv'));
sr.z     = double(nc_varget(strcat('JOBS/',Wp.name,'/OUTPUT/',Wp.name,'/OUTPUT/',Wp.filename),'zu_3d'));

clear ll M tline fid MaxTime nw uk vk xi yi f

 