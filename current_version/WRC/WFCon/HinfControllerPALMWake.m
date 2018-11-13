clear all
close all
clc

addpath(genpath('classes'));


%% parameter

Parameter.Time.dt                                       = 5;
Parameter.Wind.URef                                     = 8;
Parameter.Turbine.RotorDiameter                         = 120;
Parameter.LidarWakeCenterEstimationPosition             = 3*120;
Parameter.GainScheduling                                = 0;
Parameter.ControllerName                                = 'HinfController';
Parameter.WakeActuator.MaxYawRate                       = deg2rad(1);
Parameter.PadeOrder                                     = 10;
Parameter.PlotHinfController                            = true;
load('D:\05_code\Thesis\Results\PALMIdentification\PALMIdentificationWFSim_URef_8_StepID_7_MI_1_results.mat','logsout')
Parameter.Plant                                         = tf(logsout.ModelIdentification.sys)*180/pi;
Parameter.LoadExternalPlants                            = true;


HinfController = HinfThesisController('PALMHinfFirstTry',Parameter);

%%
Sensitivity = 1/(1-Parameter.Plant*tf(HinfController.Controller.continous.tf));

figure
step(1-Sensitivity)


figure
y = lsim(1-Sensitivity,0*(rand(1000,1)-0.5)/.5*20+5,0:999);
plot(y)