%% getVlos_UFWND
% Script: Gets the v_los velocety of an unfrozen windfield (class).
%
%% Usage:
% 
%  
% 
%  
%% Input:
% 
% 
% 
%% Output:
% 
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
% Steffen Raach on 16-Jun-2015
%
% (c) Universitaet Stuttgart
% 
%% Code:

function v_los = getVlos_UFWND(LidarData,UFwindfield,Parameter)

if ~isa(UFwindfield,'UnfrozenWindfield')
    error('windfield has wrong class type...')
end

% 1. internal variables
x_W             = LidarData.x_W(:);
y_W             = LidarData.y_W(:);
z_W             = LidarData.z_W(:);

nFocusDistances = LidarData.nFocusDistances;
nDataPoints   	= length(x_W);
nWeights        = length(Parameter.Lidar.a);

t               = repmat(LidarData.time.relSec, nFocusDistances,1);

if isfield(LidarData,'x_LW')
    x_LW        = repmat(LidarData.x_LW,        nFocusDistances,1); 
    y_LW        = repmat(LidarData.y_LW,        nFocusDistances,1); 
    z_LW        = repmat(LidarData.z_LW,        nFocusDistances,1);        
    x_LWdot     = repmat(LidarData.x_LWdot,     nFocusDistances,1); 
    y_LWdot     = repmat(LidarData.y_LWdot,     nFocusDistances,1); 
    z_LWdot     = repmat(LidarData.z_LWdot,     nFocusDistances,1);   
else
    [x_LW,y_LW,z_LW]= L2W(0,0,0,Parameter);
    x_LWdot     = zeros(nDataPoints,1);
    y_LWdot     = zeros(nDataPoints,1);
    z_LWdot     = zeros(nDataPoints,1);               
end        


% 2. normalized laser vector
LaserVector_W                       = [ x_W'-x_LW';
                                        y_W'-y_LW';
                                        z_W'-z_LW'];
NormLaserVector_W                   = sqrt(sum(LaserVector_W.^2));
NormedLaserVector_W                 = LaserVector_W./NormLaserVector_W(ones(3,1),:);
BackscatteredNormedLaserVector_W    = -NormedLaserVector_W;


% 3. considered points
Points_WF           = [ reshape(x_W(:,ones(nWeights,1)),1,[])+reshape(Parameter.Lidar.a(ones(1,nDataPoints),:),1,[]).*reshape(BackscatteredNormedLaserVector_W(1*ones(1,nWeights),:)',1,[]);  
                        reshape( y_W(:,ones(nWeights,1)),1,[])+reshape(Parameter.Lidar.a(ones(1,nDataPoints),:),1,[]).*reshape(BackscatteredNormedLaserVector_W(2*ones(1,nWeights),:)',1,[]);  
                        reshape( z_W(:,ones(nWeights,1)),1,[])+reshape(Parameter.Lidar.a(ones(1,nDataPoints),:),1,[]).*reshape(BackscatteredNormedLaserVector_W(3*ones(1,nWeights),:)',1,[]);
                        reshape(t(:,ones(nWeights,1)),1,[])];

% 4. Interpolation
u_W                 = UFwindfield.U(Points_WF(1,:),Points_WF(2,:),Points_WF(3,:),Points_WF(4,:));
v_W                 = UFwindfield.V(Points_WF(1,:),Points_WF(2,:),Points_WF(3,:),Points_WF(4,:));
w_W                 = UFwindfield.W(Points_WF(1,:),Points_WF(2,:),Points_WF(3,:),Points_WF(4,:));

if sum(isnan(u_W))+ sum(isnan(v_W)) + sum(isnan(w_W))>0
    x_error = Points_WF(1,~(Points_WF(1,:)>=min(UFwindfield.grid.x) & Points_WF(1,:)<=max(UFwindfield.grid.x)));
    y_error = Points_WF(2,~(Points_WF(2,:)>=min(UFwindfield.grid.y) & Points_WF(2,:)<=max(UFwindfield.grid.y)));
    z_error = Points_WF(3,~(Points_WF(3,:)>=min(UFwindfield.grid.z) & Points_WF(3,:)<=max(UFwindfield.grid.z)));
    error('hallo, stop, hier ging was falsch! out of windfield space...'); 
end
% 5. Calculation of line-of-sight velocity and weighting
RelativeWindVector_W    = [u_W;v_W;w_W]-repmat([x_LWdot';y_LWdot';z_LWdot'],1,nWeights); % wind vectors in all discrete points
v_los_d                 = dot(RelativeWindVector_W,BackscatteredNormedLaserVector_W(:, [1:nDataPoints]' * ones(1,nWeights)));
v_los                   = reshape((reshape(v_los_d,nDataPoints,nWeights)* Parameter.Lidar.f_L_d'),nDataPoints/nFocusDistances,nFocusDistances);
