function LidarData = LidarDataSimulationPALM(x,y,z,u,v,w,Trajectory,Parameter)

%% define wind field
PALMData.inputgrid.x = x;
PALMData.inputgrid.y = y;
PALMData.inputgrid.z = z;
PALMData.inputgrid.t = 0;
PALMData.inputgrid.ProvidedGrid = true;
PALMData.U           = permute(u,[4,3,2,1]);
PALMData.V           = permute(v,[4,3,2,1]);
PALMData.W           = permute(w,[4,3,2,1]);

PALMWindfield = UnfrozenWindfield('PALMWindfield',PALMData.inputgrid);
setU(PALMWindfield,PALMData.U);
setV(PALMWindfield,PALMData.V);
setW(PALMWindfield,PALMData.W);

%% lidar simulation
nTrajectories               = 1;
LidarData.time.relSec       = Trajectory.t*0;
LidarData.nFocusDistances   = Trajectory.nFocusDistances;
LidarData.nDataPoints       = length(LidarData.time.relSec);
LidarData.TrajectoryName    = Parameter.Scan.TrajectoryName;

LidarData.x_L               = repmat(Trajectory.x_L_AllDistances,nTrajectories,1);
LidarData.y_L               = repmat(Trajectory.y_L_AllDistances,nTrajectories,1);
LidarData.z_L               = repmat(Trajectory.z_L_AllDistances,nTrajectories,1);
LidarData.f                 = sqrt(LidarData.x_L.^2+LidarData.y_L.^2+LidarData.z_L.^2);

LidarData.IdxPointInTrajectory              = repmat([1:nPointsInTrajectory]',nTrajectories,1);

[LidarData.x_W,LidarData.y_W,LidarData.z_W] = L2W(LidarData.x_L,LidarData.y_L,LidarData.z_L,Parameter);

LidarData.v_los             =  getVlos_UFWND(LidarData,PALMWindfield,Parameter);