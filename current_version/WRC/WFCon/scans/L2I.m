%% L2I
% Function: transforms from the Lidar to inflow/inertial coordinates
% 
%% Usage:
% 
% [x_I,y_I,z_I] = L2I(x_L,y_L,z_L,Parameter)
%
% According to DIN 9300 the rotation order is Yaw-Pitch-Roll with the Euler
% angles [rad]:
%
% * Parameter.Lidar.YawAngle
% * Parameter.Lidar.PitchAngle
% * Parameter.Lidar.RollAngle
% 
% The translation [m] is equal to the position [x;y;z] of the Lidar in the 
% inflow/inertial system:
% Parameter.Lidar.PositionLinI
% 
%% Input:
% 
% * x_L,y_L,z_L         - Lidar coordinates
% * Parameter           - struct
%   .Lidar.YawAngle      - yaw angle of lidar coordinates [rad]
%   .Lidar.PitchAngle    - pitch angle of lidar coordinates [rad]
%   .Lidar.RollAngle     - roll angle of lidar coordinates [rad]
%   .Lidar.PositionLinI  - position of the lidar in the inertial system [rad]
% 
%% Output:
% 
% x_I,y_I,z_I           - inflow/inertial coordinates
% 
%
%% Modified:
%
% 
%
%% ToDo:
%
%
% 
%% Created:
% David Schlipf on 27-Oct-2012
%
% (c) Universitaet Stuttgart
% 

%% Code:
function [x_I,y_I,z_I] = L2I(x_L,y_L,z_L,Parameter)

Yaw     = Parameter.Lidar.YawAngle;
Pitch   = Parameter.Lidar.PitchAngle;
Roll    = Parameter.Lidar.RollAngle;


% Yaw is a rotation around z-axis    
T_Yaw 	= [ cos(Yaw)   -sin(Yaw)	0;
            sin(Yaw) 	cos(Yaw)    0;
            0           0           1];

% Pitch is a rotation around y-axis
T_Pitch = [	cos(Pitch)	0           sin(Pitch);
            0         	1        	0;
           -sin(Pitch)	0        	cos(Pitch)];

% Roll is a rotation around x-axis
T_Roll  = [	1           0       	0;
            0       	cos(Roll)  -sin(Roll);
            0          	sin(Roll)	cos(Roll)];

        
T       = T_Yaw*T_Pitch*T_Roll;


x_R     = T(1,1)*x_L + T(1,2)*y_L + T(1,3)*z_L;
y_R     = T(2,1)*x_L + T(2,2)*y_L + T(2,3)*z_L;
z_R     = T(3,1)*x_L + T(3,2)*y_L + T(3,3)*z_L;

% Translation 
x_I     = x_R + Parameter.Lidar.PositionLinI(1);
y_I     = y_R + Parameter.Lidar.PositionLinI(2);
z_I     = z_R + Parameter.Lidar.PositionLinI(3);
