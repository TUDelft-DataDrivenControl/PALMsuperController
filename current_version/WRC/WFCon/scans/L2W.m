%% L2W
% Function: transforms from the Lidar to wind coordinates
% 
%% Usage:
% 
% [x_W,y_W,z_W] = L2W(x_L,y_L,z_L,Parameter)
%
% 
%% Input:
% 
% * x_L,y_L,z_L      - Lidar coordinates
% * Parameter        - struct
% - Parameter.Lidar.YawAngle
% - Parameter.Lidar.PitchAngle
% - Parameter.Lidar.RollAngle
% - Parameter.Lidar.PositionLinI
% - Parameter.Wind.Elevation
% - Parameter.Wind.Azimuth
% - Parameter.Wind.PositionIinW
% 
%% Output:
% 
%  x_W,y_W,z_W      - wind coordinates
% 
%% Modified:
% * David Schlipf on 17-May-2014
% - Update Parameter
%
%% ToDo:
%
% 
%% Created:
% David Schlipf on      27-Oct-2012
%
% (c) Universitaet Stuttgart
% 

function [x_W,y_W,z_W] = L2W(x_L,y_L,z_L,Parameter)
[x_I,y_I,z_I] = L2I(x_L,y_L,z_L,Parameter);
[x_W,y_W,z_W] = I2W(x_I,y_I,z_I,Parameter);