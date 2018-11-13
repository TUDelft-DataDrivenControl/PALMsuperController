classdef UnfrozenWindfield < handle
    properties (SetAccess = protected)
        Name
        grid
    end
    properties
        U
        V
        W
    end
    methods
        function WF = UnfrozenWindfield(Name,inputgrid)
            WF.Name = Name;
            setgrid(WF,inputgrid);
        end
        function setgrid(WF,inputgrid)
            if isfield(inputgrid,'ProvidedGrid') && inputgrid.ProvidedGrid
                % grid is externally provided
                WF.grid = inputgrid;
            else
            FieldNames = {'nt','nx','ny','nz','dt','dx','dy','dz','t0','x0','y0','z0'};
            for iFieldName = 1:length(FieldNames)
                if ~isfield(inputgrid,FieldNames{iFieldName});
                    error(['grid is not fully defined, ' FieldNames{iFieldName} ' is missing...']);
                end
            end
            WF.grid     = inputgrid;
            % calculate grid
            if ~isfield(WF.grid,'t')
                WF.grid.t   = WF.grid.t0:WF.grid.dt:(WF.grid.nt-1)*WF.grid.dt+WF.grid.t0;
                WF.grid.x   = WF.grid.x0:WF.grid.dx:(WF.grid.nx-1)*WF.grid.dx+WF.grid.x0;
                WF.grid.y   = linspace(-WF.grid.ny*WF.grid.dy/2,WF.grid.ny*WF.grid.dy/2,WF.grid.ny); % has perhaps to be adjusted
                WF.grid.z   = WF.grid.z0:WF.grid.dz:(WF.grid.nz-1)*WF.grid.dz+WF.grid.z0;
            end
            end
            [WF.grid.X,WF.grid.Y,WF.grid.Z,WF.grid.T] = ndgrid(WF.grid.x,WF.grid.y,WF.grid.z,WF.grid.t);            
        end
        
        % set velocities
        function setU(WF,U)
            WF.U = setVelocity(WF,U);
        end
        function setV(WF,V)
            WF.V = setVelocity(WF,V);
        end
        function setW(WF,W)
            WF.W = setVelocity(WF,W);
        end
        function Interpolant = setVelocity(WF,VEL)
            Interpolant = griddedInterpolant(WF.grid.X,WF.grid.Y,WF.grid.Z,WF.grid.T,VEL);
            Interpolant.ExtrapolationMethod = 'none';
        end
    end
end
        
        
        
        
        
%         windfield.grid.nt   = nTimeSteps;
%         windfield.grid.ny   = ny;
%         windfield.grid.nz   = nz;
% 
%         windfield.ny        = ny;
%         windfield.nz        = nz;
%         windfield.dt        = SliceDT;
% 
%         windfield.grid.dt   = SliceDT;
%         windfield.grid.dy   = dy; 
%         windfield.grid.dz   = dz;
% 
%         windfield.grid.t    = [0:windfield.dt:windfield.dt*(windfield.grid.nt-1)]';
% 
%         windfield.grid.z    = (0:dz:(nz-1)*dz) - HubHeight;
%         windfield.grid.y    = -dy*(ny-1)/2:dy:dy*(ny-1)/2;
%         
%         [windfield.grid.Y, windfield.grid.Z] = meshgrid(windfield.grid.y,windfield.grid.z);
%         
%         windfield.u         = zeros(windfield.ny,nTimeSteps,windfield.nz);
%         windfield.v         = zeros(windfield.ny,nTimeSteps,windfield.nz);
%         windfield.w         = zeros(windfield.ny,nTimeSteps,windfield.nz);
%         
%         for iTimeStep=1:nTimeSteps
%             interpolantU = scatteredInterpolant(y_L{iDistance,iTurbine},z_L{iDistance,iTurbine},u_L{iTimeStep,iDistance,iTurbine});
%             interpolantV = scatteredInterpolant(y_L{iDistance,iTurbine},z_L{iDistance,iTurbine},v_L{iTimeStep,iDistance,iTurbine});
%             interpolantW = scatteredInterpolant(y_L{iDistance,iTurbine},z_L{iDistance,iTurbine},w_L{iTimeStep,iDistance,iTurbine});
%             windfield.u(:,iTimeStep,:)         = interpolantU(windfield.grid.Y,windfield.grid.Z)';
%             windfield.v(:,iTimeStep,:)         = interpolantV(windfield.grid.Y,windfield.grid.Z)';
%             windfield.w(:,iTimeStep,:)         = interpolantW(windfield.grid.Y,windfield.grid.Z)';
%         end