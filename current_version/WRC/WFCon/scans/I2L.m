%% I2L
% Function: transforms from the inflow/inertial to lidar coordinates
% 
%% Usage:
% 
% [x_L,y_L,z_L] = I2L(x_I,y_I,z_I,Parameter);
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
% * x_I,y_I,z_I         - inflow/inertial coordinates
% * Parameter           - struct
% 
%% Output:
% 
%   x_L,y_L,z_L         - Lidar coordinates
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
% David Schlipf on 03-Aug-2014
%
% (c) Universitaet Stuttgart
% 

function [x_L,y_L,z_L] = I2L(x_I,y_I,z_I,Parameter)

Yaw     = Parameter.Lidar.YawAngle;
Pitch   = Parameter.Lidar.PitchAngle;
Roll    = Parameter.Lidar.RollAngle;

% Translation 
x_R     = x_I - Parameter.Lidar.PositionLinI(1);
y_R     = y_I - Parameter.Lidar.PositionLinI(2);
z_R     = z_I - Parameter.Lidar.PositionLinI(3);

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

        
T       = inv(T_Yaw*T_Pitch*T_Roll);


x_L     = T(1,1)*x_R + T(1,2)*y_R + T(1,3)*z_R;
y_L     = T(2,1)*x_R + T(2,2)*y_R + T(2,3)*z_R;
z_L     = T(3,1)*x_R + T(3,2)*y_R + T(3,3)*z_R;
