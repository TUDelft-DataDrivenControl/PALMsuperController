%% H2L
% Function: transforms from the hub to Lidar coordinates
% 
% 
%% Usage:
% 
% [x_L,y_L,z_L] = H2L(x_H,y_H,z_H,Parameter)
% 
%% Input:
% 
%  * x_H,y_H,z_H      - hub coordinates
%  * Parameter        -   struct
% 
%% Output:
% 
%  x_L,y_L,z_L     	- Lidar coordinates
%
%% Modified:
%
% * David Schlipf on      01-Jul-2011
% - new nomenclature
% * David Schlipf on      14-Apr-2011
% - add elevation and azimuth
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

function [x_L,y_L,z_L] = H2L(x_H,y_H,z_H,Parameter)

% Translation 
x_T=x_H-Parameter.Lidar.Position(1);
y_T=y_H-Parameter.Lidar.Position(2);
z_T=z_H-Parameter.Lidar.Position(3);


EL=Parameter.Lidar.Elevation;
AZ=Parameter.Lidar.Azimuth;

% elevation is a rotation around y-axis
T_EL=[  cosd(EL)    0           sind(EL);
        0           1           0;
        -sind(EL)   0           cosd(EL)];

% azimuth is a rotation around z-axis    
T_AZ=[  cosd(AZ)    -sind(AZ)   0;
        sind(AZ)    cosd(AZ)    0;
        0           0       	1];

T=inv(T_EL)*inv(T_AZ);

x_L=T(1,1)*x_T+T(1,2)*y_T+T(1,3)*z_T;
y_L=T(2,1)*x_T+T(2,2)*y_T+T(2,3)*z_T;
z_L=T(3,1)*x_T+T(3,2)*y_T+T(3,3)*z_T;
                