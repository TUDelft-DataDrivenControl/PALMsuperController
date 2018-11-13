clear;clc;

Wp.name          = '9turb';        % farm
Wp.controller    = 'SMPC';

addpath(genpath(strcat(Wp.controller,'/WFCon')));
addpath(genpath(strcat(Wp.controller,'/libraries')));
addpath(genpath(strcat('USER_CODE/6turb/matlab/mexcdf')));

filename         = strcat(Wp.name,'_m01.nc');
Wp.Nt            = 9;                        % #turbines
nz               = 1;                        % z-index of hubs (this should coincide with z-mask defined in _p3d)

Wp.N             = 2*floor(800/2);           % simulation time ( real simulation time is N-h )
Wp.h             = 0.5;                      % PALM sample period (measure and apply control signals every 2h)
Wp.k             = 0;                        % time index
Wp.N0            = round(.1*Wp.N/(2*Wp.h));  % First n% of the simulation constant reference
WrongMeasurement = 0;                        % Amount of failed measurements
MaxTime          = [30000;...
                  60*ones(Wp.N/(2*Wp.h),1)]; % Maximum waiting time for measurement
              
Wp.Ny    = 5;                                % #measurements [time,turbine_number,F,P,Ur]
Wp.y     = zeros(Wp.N/(2*Wp.h),Wp.Nt,Wp.Ny);

% initialize controller parameters/variables
[gp,ap,sr] =  InitializeController(Wp);

% clean up
unix(char(strcat('rm -f ',{' '},Wp.controller,'/Data/*')));
unix(strcat('rm -f JOBS/',Wp.name,'/MONITORING/*'));
unix(strcat('rm -f USER_CODE/',Wp.name,'/matlab/* 2> /dev/null'));
unix(strcat('rm -f USER_CODE/',Wp.name,'/matlab/.nfs*'));
unix(strcat('rm -f USER_CODE/',Wp.name,'/matlab/previous_y/*'));
unix(strcat('rm -f JOBS/',Wp.name,'/OUTPUT/*'));
%unix(strcat('rm -f ../../job_queue/*'));

% set simulation time in _p3d file
unix(char(strcat('sed -i "/&d3par/c\&d3par end_time =',{' '},num2str(Wp.N,'%.1f'),'" JOBS/',Wp.name,'/INPUT/',Wp.name,'_p3df')));

% set sample period _p3d file
unix(char(strcat('sed -i "/dt =/c\         dt =',{' '},num2str(Wp.h,'%.1f'),'" JOBS/',Wp.name,'/INPUT/',Wp.name,'_p3df')));
unix(char(strcat('sed -i "/dt_run_control =/c\         dt_run_control =',{' '},num2str(2*Wp.h,'%.1f'),'" JOBS/',Wp.name,'/INPUT/',Wp.name,'_p3df')));
unix(char(strcat('sed -i "/dt_domask =/c\         dt_domask =',{' '},num2str(2*Wp.h,'%.1f'),'" JOBS/',Wp.name,'/INPUT/',Wp.name,'_p3df')));

% set correct yaw angles in _p3d file
str = [];
for kk=1:Wp.Nt
    str = strcat(str,num2str(ap.yaw(kk),'%.1f'),',');
end
unix(char(strcat('sed -i "/phi_yaw =/c\         phi_yaw =',{' '},str,'" JOBS/',Wp.name,'/INPUT/',Wp.name,'_p3df')));


% start simulation PALM
[status,cmdout] = ...
    unix(char(strcat('mrun -d',{' '},Wp.name,' -h lchpc06 -K parallel -X 16 -T 8 -t 360000 -m 64000 -v -b -q dcsc -r "d3f maf"')));

% start loop
while Wp.k<Wp.N/(2*Wp.h)
    
    tic
    Wp.k = Wp.k + 1;
    
    
    % 1) run controller
    
    [sr,ap]        = consys_cent(Wp.k,sr,ap,gp);

    if Wp.k==1;    sr.u(:,Wp.k,:) = 2.0;end % start with initial developed flow (CT'=2.0)
         
    uc             = sr.u(:,Wp.k);                % CT'
    Wp.CT(Wp.k,:)  = 1 - 4*(uc./(4 + uc)-1/2).^2; % CT
    
    
    % 2) write control signals to a text file ( are applied at k=0,2h,...,(N-1)*h and hold until new ones are send )
    
    for ll=1:Wp.Nt
        dlmwrite(strcat('USER_CODE/',Wp.name,'/matlab/u0',num2str(ll),'.txt'),Wp.CT(Wp.k,ll),'precision',8)
    end
    dlmwrite(strcat('USER_CODE/',Wp.name,'/matlab/uflag.txt'),1) % PALM waits for this flag
    fclose('all');
    
    % 3) wait until measurement from PALM is written to file
    temp = zeros(1,Wp.Nt);
    while prod(temp)==0
        if toc > MaxTime(Wp.k)
            for ll=1:Wp.Nt
                unix(char(strcat('mv USER_CODE/',Wp.name,'/matlab/previous_y/y0',num2str(ll),'.txt USER_CODE/',Wp.name,'/matlab/y0',num2str(ll),'.txt')));
            end
            WrongMeasurement = WrongMeasurement+1;
            pause(.3)
            break;
        end
        pause(.3)
        for ll=1:Wp.Nt
            temp(ll) = exist(strcat('USER_CODE/',Wp.name,'/matlab/y0',num2str(ll),'.txt'),'file');
        end
    end;

    
    % 4) read measurement from PALM ( at k=0,2h,...,(N-1)*h )
    
    for ll=1:Wp.Nt
        % [time, turbine_number, force, power, Ur]
        fid   = fopen(strcat('USER_CODE/',Wp.name,'/matlab/y0',num2str(ll),'.txt')); % y(time,turbine,measurement)
        tline = fgetl(fid);
        timer = tic;
        while tline==-1
            tline = fgetl(fid);
            if toc(timer) > MaxTime(Wp.k)
                tline            = num2str(Wp.y(Wp.k-1,ll,:));
                WrongMeasurement = WrongMeasurement+1;
                break
            end
        end
        Wp.y(Wp.k,ll,:)  = str2num(tline); % y(time,turbine,measurement)
        fclose(fid);    
    end    

    % backup measurement
    for ll=1:Wp.Nt
        unix(char(strcat('mv USER_CODE/',Wp.name,'/matlab/y0',num2str(ll),'.txt USER_CODE/',Wp.name,'/matlab/previous_y/y0',num2str(ll),'.txt')));
    end 
    
    % store measurements in controller variables
    Wp.y(Wp.k,:,3)       = -Wp.y(Wp.k,:,3);
    xmes                 = [squeeze(Wp.y(Wp.k,:,3:4)) uc];  %=[F P uf]
    vmes                 = squeeze(Wp.y(Wp.k,:,5));         %=Ur
    umes                 = uc;   
    ymes                 = xmes;
    
    sr.e(Wp.k)         = gp.Pnref(Wp.k+1) - sum(ymes(:,2));   % wind farm error      
    sr.x(:,Wp.k+1)     = [reshape(xmes',size(xmes,2)*size(xmes,1),1); sr.e(Wp.k)]; % wind farm state
    sr.v(:,Wp.k+1)     = vmes'; % averaged rotor flow velocities
    sr.u(:,Wp.k+1)     = umes'; % CT'
    sr.y(:,Wp.k)       = sr.x(:,Wp.k+1); % wind farm state  

    
    
    save(strcat(Wp.controller,'/Data/',num2str(Wp.k)),'sr','gp','ap','Wp');
    
    disp(['Simulated t(' num2str(Wp.k) ')  = ' num2str(Wp.y(Wp.k,1,1),3) ...
        ' s.       CPU: ' num2str(toc,3) ' s.' ...
        '          error: ' num2str( (gp.Pnref(Wp.k)-sum(ymes(:,2)))/1e6,3 ) ' MW.']);
    
end
tic
% fetch measurements from .nc file
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

% Time   UR  Uinf  Ct_adm  a Yaw Thrust Power
for ll=1:Wp.Nt
    
    unix(strcat('sed -i "s/\*/0/g" JOBS/',Wp.name,'/MONITORING/',Wp.name,'_turbine_parameters0',num2str(ll),'.txt'));
    
    M = dlmread(strcat('JOBS/',Wp.name,'/MONITORING/',Wp.name,'_turbine_parameters0',num2str(ll),'.txt'),'',1,0);
    
    PALM.Power(ll,:) = M(:,8)';
    PALM.phi(ll,:)   = M(:,6)';
    PALM.time(ll,:)  = M(:,1)';
    PALM.CT(ll,:)    = M(:,4)';
    PALM.Ur(ll,:)    = M(:,2)';
    PALM.aif(ll,:)   = M(:,5)';
end

PALM.t       = double(nc_varget(strcat('JOBS/',Wp.name,'/OUTPUT/',filename),'time'));
temp         = double(nc_varget(strcat('JOBS/',Wp.name,'/OUTPUT/',filename),'u'));
PALM.u       = squeeze(temp(:,nz,:,:));
temp         = double(nc_varget(strcat('JOBS/',Wp.name,'/OUTPUT/',filename),'v'));
PALM.v       = squeeze(temp(:,nz,:,:));

% clear up before saving
clear temp M ll filename cmdout status xmes umes vmes ymes uc kk str MaxTime timer tline fid PALM

pause(10)
save('series_smpc')