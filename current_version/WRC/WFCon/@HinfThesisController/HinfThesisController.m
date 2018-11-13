classdef HinfThesisController < handle
    properties (SetAccess = protected)
        Name    % name of the instance
        URef    % reference wind speed
        dt      % dt for the controller
        RotorDiameter   % rotor diameter
        Controller      % Controller struct
        Filter      % filter struct
        LidarWakeCenterEstimationPosition % position, at which the wake center line is estimated (delta x behinde turbine)
        InternalModel   % internal model for wake deflection
        Plants
        Model
        ControllerName
    end
    properties (SetAccess = private)
        IntegratorState
        InternalModelOutput
        InternalModelDelayedOutput
        InternalModelBuffer
        InternalModelYaw
    end
    properties %public
        DesiredWakeCenterPosition % delta y with respect to the center line
    end
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        %% Initialization of the class
        function CLWSC = HinfThesisController(Name,Parameter)
            
            % basic parameter
            CLWSC.Name                                  = Name;
            CLWSC.dt                                    = Parameter.Time.dt;
            CLWSC.URef                                  = Parameter.Wind.URef;
            CLWSC.RotorDiameter                         = Parameter.Turbine.RotorDiameter;
            CLWSC.LidarWakeCenterEstimationPosition     = Parameter.LidarWakeCenterEstimationPosition;
            
            CLWSC.ControllerName                        = Parameter.ControllerName;
            CLWSC.Controller.GainScheduling             = Parameter.GainScheduling;
            CLWSC.Controller.MaxYawRate                 = Parameter.WakeActuator.MaxYawRate;
            
            CLWSC.InternalModel.PadeOrder               = Parameter.PadeOrder;
            
            % internal model            
            CLWSC.InternalModel.xi_init                 = @(Yaw,CT)1/2*cos(Yaw)^2*sin(Yaw)*CT;
            CLWSC.InternalModel.WakeCenterPositionFct   = @(Yaw,delta_dist,k_d,CT,RotorDiameter)CLWSC.InternalModel.xi_init(Yaw,CT)*(15*((2*k_d*delta_dist)/RotorDiameter+1)^4+CLWSC.InternalModel.xi_init(Yaw,CT)^2)/(30*k_d/RotorDiameter*(2*k_d*delta_dist/RotorDiameter+1)^5) - CLWSC.InternalModel.xi_init(Yaw,CT)*RotorDiameter*(15+CLWSC.InternalModel.xi_init(Yaw,CT)^2)/(30*k_d);
                        
            % important quantities
            CLWSC.InternalModel.k_d = 0.15;
            CLWSC.InternalModel.CT  = 0.7;
            CLWSC.InternalModel.Induction = .25;
            DelayTime                       = GetDelayTime(CLWSC);
            StaticGain                      = -GetInternalModelGain(CLWSC,deg2rad(20));
            CLWSC.InternalModel.DelayTime   = DelayTime;
            
            [Delay.num,Delay.den]           = pade(DelayTime,CLWSC.InternalModel.PadeOrder);
            Delay.tf = tf(Delay.num,Delay.den);
            
            if isfield(Parameter,'LoadExternalPlants') && Parameter.LoadExternalPlants
                CLWSC.Plants = Parameter.Plant;
            else
                CLWSC.Plants = StaticGain*Delay.tf;
            end
            Plant = CLWSC.Plants;
                        
            switch CLWSC.ControllerName
                case 'HinfController'
                    try
                        Parameter.PlotHinfController;
                    catch
                        Parameter.PlotHinfController = false;
                    end
                    %% design weights
                    s = tf('s');
                    display(['Delay time = ',num2str(DelayTime)])
                    wB1=1/DelayTime/2;     % desired closed-loop bandwidth
                    A=1/1000;     % desired disturbance attenuation inside bandwidth
                    M=2 ;       % desired bound on hinfnorm(S) & hinfnorm(T)
                    Wp=[(s/M+wB1)/(s+wB1*A) ]; % Sensitivity weight
                    
                    %Wu=[0.5*tf([1],[1]) ];
                    %omega=.01;
                    %Wu=[20*tf([1 2*(omega)*0.7 (omega)^2],[1 2*(10*omega)*0.7 (10*omega)^2])*100 ]; % Control weight
                    omega=0.1; omega2=0.01;
                    %Wu=[20*tf([1 2*(omega)*0.7 (omega)^2],[1 2*(10*omega)*0.7 (10*omega)^2])*100*tf([1 2*(omega2)*0.7 (omega2)^2],[1 2*(10*omega2)*0.7 (10*omega2)^2])*100 ]; % Control weight
                    %Wu=0*[0.5*150*tf([1 2*(omega)*0.7 (omega)^2],[1 2*(1000*omega)*0.7 (1000*omega)^2])*1000^2*tf([1 2*(omega2)*0.7 (omega2)^2],[1 2*(1000*omega2)*0.7 (1000*omega2)^2])*1000^2 ]; % Control weight
                    Wu=[0.05*40*tf([1 2*(omega)*0.7 (omega)^2],[1 2*(10*omega)*0.7 (10*omega)^2])*100*tf([1 2*(omega2)*0.7 (omega2)^2],[1 2*(10*omega2)*0.7 (10*omega2)^2])*100 ]; % Control weight

                      Wt=0;                  % Empty weight
                    %% Uses the robust control toolbox
                    T = iconnect;
                    u = icsignal(1);
                    w = icsignal(1);
                    output = icsignal(4);
                    T.Input = [w;u];
                    T.Equation{1} = equate(output,[Wp*(w+Plant*u);Wt*Plant*u;Wu*u;Plant*u+w]);
                    T.Output = [output];
                    
                    P = T.System;
                    
                    [K2,CL2,GAM33,INFO2] = hinfsyn(P,1,1);
                    
                    % plotting for analysis
                    if Parameter.PlotHinfController
                        % plots
                        figure('Name','open-loop');
                        margin(-K2*Plant); grid on;
                        % plots
                        figure
                        subplot(2,2,1)
                        bodemag(feedback(eye(1),series(K2,Plant),1,1,1)); hold on;grid on;
                        bodemag(inv(Wp))
                        title('sensitivity')
                        
                        subplot(2,2,2)
                        bodemag(K2*feedback(eye(1),series(K2,Plant),1,1,1)); hold on; grid on;
                        bodemag(inv(Wu))
                        title('performance')
                        
                        subplot(2,2,3)
                        bode(K2);grid on;
                        title('controller')
                        
                        subplot(2,2,4)
                        bodemag(feedback(series(K2,Plant),1,1,1,1)); hold on; grid on;
                        title('inverse sensitivity')
                        
                    end
                    % assignment
                    CLWSC.Controller.continous.tf   = tf(minreal(K2));
                    CLWSC.Controller.discrete.tf    = c2d(K2,CLWSC.dt);
                    CLWSC.Controller.discrete.ss    = ss(CLWSC.Controller.discrete.tf);
                    
                    CLWSC.IntegratorState               = zeros(length(CLWSC.Controller.discrete.ss.a),1);
                    CLWSC.DesiredWakeCenterPosition     = 0;
                    
            end
            % init outputs            
            CLWSC.InternalModelOutput           = 0;
            CLWSC.InternalModelDelayedOutput    = 0;
            CLWSC.InternalModelYaw              = 0;
            % output
            display(['CLWakeSteeringController: ' CLWSC.Name ' has been successfully created.']);
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        %% Helpers for initialization
        
        function [InternalModelGain] = GetInternalModelGain(CLWSC,Yaw_0)
            
            k_d = CLWSC.InternalModel.k_d;
            CT  = CLWSC.InternalModel.CT;
            LidarWakeCenterEstimationPosition = CLWSC.LidarWakeCenterEstimationPosition;
            RotorDiameter = CLWSC.RotorDiameter;
            
            InternalModelGain = (CLWSC.InternalModel.WakeCenterPositionFct(Yaw_0,LidarWakeCenterEstimationPosition,k_d,CT,RotorDiameter) - CLWSC.InternalModel.WakeCenterPositionFct(0,LidarWakeCenterEstimationPosition,k_d,CT,RotorDiameter))/(Yaw_0 - 0);
%             syms Yaw
%             Yaw_subs = (Yaw_0);
            
%             xi_init     = 1/2*cos(Yaw)^2*sin(Yaw)*CT;
%             xi_init*(15*((2*k_d*LidarWakeCenterEstimationPosition)/RotorDiameter+1)^4+xi_init^2)/(30*k_d/RotorDiameter*(2*k_d*LidarWakeCenterEstimationPosition/RotorDiameter+1)^5) - xi_init*RotorDiameter*(15+xi_init^2)/(30*k_d)
%             gain = subs(diff(xi_init*(15*((2*k_d*LidarWakeCenterEstimationPosition)/RotorDiameter+1)^4+xi_init^2)/(30*k_d/RotorDiameter*(2*k_d*LidarWakeCenterEstimationPosition/RotorDiameter+1)^5) - xi_init*RotorDiameter*(15+xi_init^2)/(30*k_d),Yaw),Yaw,Yaw_subs);
%             InternalModelGain   = eval(gain);
        end
        
        function [DelayTime] = GetDelayTime(CLWSC)
            DT = 1e-3;
            x = 0;
            iTime = 0;
            InductionFct = @(x,RotorDiameter,Induction)2/pi*Induction*(pi/2+atan(2*x/RotorDiameter)); 
            while x < CLWSC.LidarWakeCenterEstimationPosition
                x = x + DT * ((1-InductionFct(x,CLWSC.RotorDiameter,CLWSC.InternalModel.Induction))*CLWSC.URef);
                iTime = iTime + 1;
            end
            DelayTime = iTime * DT;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% main functions for CLWakeSteeringController
        
        function DemandedYawAngle = GetDemandedYawAngle(CLWSC,WakeCenterPosition,varargin)
            if nargin==3
                Yaw = varargin{1};
            else
                Yaw = NaN;    
            end
            
            switch CLWSC.ControllerName          
                case 'HinfController'
                    Error                   =  WakeCenterPosition - CLWSC.DesiredWakeCenterPosition;
                    CLWSC.IntegratorState   = CLWSC.Controller.discrete.ss.a * CLWSC.IntegratorState + CLWSC.Controller.discrete.ss.b * Error;
                    DemandedYawAngle        = CLWSC.Controller.discrete.ss.c * CLWSC.IntegratorState;
            end
            
        end
        
        function YawAngle = CLWakeSteeringStep(CLWSC,WakeCenterPosition)
            % call controller
            DemandedYawAngle    = GetDemandedYawAngle(CLWSC,WakeCenterPosition);
            % yaw actuator
            YawAngle            = CalculateFilterOutput(CLWSC.YawActuator.Handle,CLWSC.YawActuator.num,CLWSC.YawActuator.den,DemandedYawAngle);
        end
        
        
    end
end