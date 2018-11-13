%% L2H
% Function: transforms from the Lidar to hub coordinates
% 
%% Usage:
% 
% [x_H,y_H,z_H] = L2H(x_L,y_L,z_L,Parameter)
%
% Parameter.Lidar.Position
%
% Parameter.Lidar.Elevation
%
% Parameter.Lidar.Azimuth
% 
%% Input:
% 
% * x_L,y_L,z_L      - Lidar coordinates
% * Parameter        - struct
% 
%% Output:
% 
% * x_H,y_H,z_H      - hub coordinates
% 
%
%% Modified:
%
% * David Schlipf on      01-Jul-2011
% - new nomenclature
% * David Schlipf on      14-Apr-2011
% - add Elevation and Azimuth
%
%% ToDo:
%
% 
% 
%% Created:
% David Schlipf on      22-Sep-2009
%
% (c) Universitaet Stuttgart
% 

function [x_H,y_H,z_H] = L2H(x_L,y_L,z_L,Parameter)


EL=Parameter.Lidar.Elevation;
AZ=Parameter.Lidar.Azimuth;

% Elevation is a rotation around y-axis
T_EL=[  cosd(EL)    0           sind(EL);
        0           1           0;
        -sind(EL)   0           cosd(EL)];

% Azimuth is a rotation around z-axis    
T_AZ=[  cosd(AZ)    -sind(AZ)   0;
        sind(AZ)    cosd(AZ)    0;
        0           0       	1];

T=T_AZ*T_EL;

x_R=T(1,1)*x_L+T(1,2)*y_L+T(1,3)*z_L;
y_R=T(2,1)*x_L+T(2,2)*y_L+T(2,3)*z_L;
z_R=T(3,1)*x_L+T(3,2)*y_L+T(3,3)*z_L;

% Translation 
x_H=x_R+Parameter.Lidar.Position(1);
y_H=y_R+Parameter.Lidar.Position(2);
z_H=z_R+Parameter.Lidar.Position(3);
