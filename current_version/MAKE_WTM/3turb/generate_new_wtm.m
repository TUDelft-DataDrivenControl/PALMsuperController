clear;clc;


N = 3;              % number of turbines

CT = [.89*ones(1,250);
    .89*ones(1,250); 
    .89*ones(1,250)];


%% Create turbine model for PALM
CTT  = [];
for kk=1:size(CT,2)
    CTT = [CTT CT(:,kk)'];
end
[m,n] = size(CT);
% Make a new turbine model
fileID = fopen('wind_turbine_model_mod.f90','wb');
% Write first part wt model in new turbine model
copyfile wind_turbine_model_mod1.f90 wind_turbine_model_mod.f90
% Append with new input signal
fileID = fopen('wind_turbine_model_mod.f90','ab');
fprintf(fileID,'REAL(wp), DIMENSION(%d,%d) :: A\n\n',m,n);
fprintf(fileID,'       A = RESHAPE((/');
for kk=1:size(CTT,2)-1
    fprintf(fileID,'%6.4f,',CTT(kk));
end
fprintf(fileID,'%6.4f',CTT(end));
fprintf(fileID,'/),SHAPE(A))\n');
fclose(fileID);
% Read final part of turbine model
fileID = fopen('wind_turbine_model_mod2.f90','r');
% Get number of lines final part turbine model
numLines = 0;
while true
    t = fgetl(fileID);
    if ~ischar(t)
        break;
    else
        numLines = numLines + 1;
    end
end
fclose(fileID); fileID = fopen('wind_turbine_model_mod2.f90','r');
% Store final part of turbine model in array
mydata = cell(1, numLines);
for k = 1:numLines
    mydata{k} = fgetl(fileID);
end
fclose(fileID);
% Append the final part to the new turbine model
fileID = fopen('wind_turbine_model_mod.f90','a');
fprintf(fileID,'%s\n',mydata{:});
% Close the file
fclose('all');
% Move new turbine model to right folder
movefile('wind_turbine_model_mod.f90','../../USER_CODE/3turb')
% Set correct simulation time
unix(char(strcat('sed -i "/&d3par/c\&d3par end_time =',{' '},num2str(n,'%.1f'),'" ../../JOBS/3turb/INPUT/3turb_p3df')));


%% Start PALM simulation
% restart run
%[status,cmdout] = unix('(cd ../../;mrun -d 3turb -h lchpc06 -K parallel -X 16 -T 8 -t 3600000 -m 64000 -v -b -q dcsc -r "d3# ma# restart")');
% restart job
[status,cmdout] = unix('(cd ../../;mrun -d 3turb -h lchpc06 -K parallel -X 16 -T 8 -t 3600000 -m 64000 -v -b -q dcsc -r "d3f maf")');
