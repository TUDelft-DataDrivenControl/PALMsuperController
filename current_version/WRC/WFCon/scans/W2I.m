%% W2I
% Function: transforms from the wind to intertial/inflow coordinates
% 
% 
%% Usage:
%  
% [x_I,y_I,z_I] = W2I(x_W,y_W,z_W,Parameter)
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
%  * x_W,y_W,z_W        - wind coordinates
%  * Parameter          - struct
%  
%% Output:
%  
%  x_I,y_I,z_I          - inflow/inertial coordinates
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

function [x_I,y_I,z_I] = W2I(x_W,y_W,z_W,Parameter)
if isfield(Parameter.Wind,'PositionIinW')
    % Translation 
    x_T     = x_W - Parameter.Wind.PositionIinW(1);
    y_T     = y_W - Parameter.Wind.PositionIinW(2);
    z_T     = z_W - Parameter.Wind.PositionIinW(3);
else
	% No Translation 
    x_T     = x_W;
    y_T     = y_W;
    z_T     = z_W;
end

% Translation 


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

T       = T_AZ*T_EL;

x_I     = T(1,1)*x_T+T(1,2)*y_T+T(1,3)*z_T;
y_I     = T(2,1)*x_T+T(2,2)*y_T+T(2,3)*z_T;
z_I     = T(3,1)*x_T+T(3,2)*y_T+T(3,3)*z_T;
                