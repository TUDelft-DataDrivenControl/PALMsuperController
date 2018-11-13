%% I2W
% Function: transforms from the intertial/inflow to wind coordinates
% 
% 
%% Usage:
%  
% [x_W,y_W,z_W] = I2W(x_I,y_I,z_I,Parameter)
% The rotation order is Azimuth-Elevation with angles [rad]:
%
% * Parameter.Wind.Elevation
% * Parameter.Wind.Azimuth
%
% The translation [m] is equal to the position [x;y;z] of the inflow/
% inertial system in the wind system
% Parameter.Wind.PositionIinW
%  
%% Input:
%  
%  * x_I,y_I,z_I        - inflow/inertial coordinates
%  * Parameter          - struct
%  
%% Output:
%  
%  x_W,y_W,z_W          - wind coordinates
%
% 
%%  Modified:
% * David Schlipf on 	21-Dec-2013
% - Translation & new rotation order
% * David Schlipf on	01-Jul-2011
% - new nomenclature
% * David Schlipf on	14-Apr-2011
% - add elevation and azimuth
%
%% ToDo:
% 
%
%  
%% Created:
% David Schlipf on      14-Apr-2011
%
% (c) Universitaet Stuttgart
% 

function [x_W,y_W,z_W] = I2W(x_I,y_I,z_I,Parameter)

EL      = Parameter.Wind.Elevation;
AZ      = Parameter.Wind.Azimuth;

% elevation is a rotation around y-axis
T_EL    = [	cos(EL)     0           sin(EL);
            0         	1           0;
           -sin(EL)     0           cos(EL)];

% azimuth is a rotation around z-axis    
T_AZ    = [	cos(AZ)    -sin(AZ)     0;
            sin(AZ)     cos(AZ) 	0;
            0           0       	1];

T       = inv(T_AZ*T_EL); 

x_R     = T(1,1)*x_I+T(1,2)*y_I+T(1,3)*z_I;
y_R     = T(2,1)*x_I+T(2,2)*y_I+T(2,3)*z_I;
z_R     = T(3,1)*x_I+T(3,2)*y_I+T(3,3)*z_I;

if isfield(Parameter.Wind,'PositionIinW')
    % Translation 
    x_W     = x_R + Parameter.Wind.PositionIinW(1);
    y_W     = y_R + Parameter.Wind.PositionIinW(2);
    z_W     = z_R + Parameter.Wind.PositionIinW(3);
else
	% No Translation 
    x_W     = x_R;
    y_W     = y_R;
    z_W     = z_R;
end
