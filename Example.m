%% This file is an example to show how to use the functions, including:
%   DEP_evaluation.m: Function to evaluate deployable ERES power based on the rail operational data
%   Aggregation_model.m: Function to generate aggregation model from the ERES capacity
%   PlotResult.m: A simple plot function to show the result
%   TrainConfiguration_example.mat: The example data for functions, including rail operational data of 50 trains

%% 1, Deployable ERES power evaluation
% Example: Evaluating 5-minute upward power at 10:00 from example trains

% Input:
%   TrainConfiguration:1*N structure:
%       TrainPara:1*10 double
%           1 Fm        N/kg
%           2 Pm        W/kg
%           3 gamma    
%           4 a0        under m/s and unit mass
%           5 a1        under m/s and unit mass
%           6 a2        under m/s and unit mass
%           7 Bm        N/kg
%           8 m     kg
%           9 etaT
%           10 etaB
%       RoadPara:m*2 double
%           1 SpeedLimit        m/s
%           2 Length        m
%       DepartureTime       s
%       GuidancetRef 1-row list     s
%       GuidanceVRef 1-row list     m/s
%       GuidanceXRef 1-row list     m
%       GuidanceuRef 1-row list 
%       Addmissiblet:1*2 double
%           1 Start time        s
%           2 End time      s
%       Addmissiblex:1*2 double
%           1 Start position        m
%           2 End position      m
load('TrainConfiguration_example.mat');
TrainConfiguration = TrainConfiguration_example;
%   RegSet:1*5 cell:
%       list_Trainindex: Index of trains participating in the ERES
%       delta_t: Discrete time      s
%       T_reg: Regulation deploying period     min
%       Strategy: 1,Upward; 2,Downward
%       t_begin: Regulation start time     s
RegSet{1} = 1:50;
RegSet{2} = 60;
RegSet{3} = [0 5];
RegSet{4} = 1;
RegSet{5} = 36000;

% Output:
%   Power: Deployable ERES power capacity  MW
%   Energy: Deployable ERES energy  capacity  MWh
%   ChangeofPower: Power adjustment of every train  MW
[Power,Energy,ChangeofPower] = DEP_evaluation(TrainConfiguration,RegSet);

%% 2, Virtual energy storage aggregation model
% Example: Aggregation model of example trains at 10:00 for upward regulation

%AGGREGATION_MODEL Aggregation model of ERES
% Input:
%   ChangeofPower_MultiplePeriod:1*n cell
%           1 ChangeofPower when regulation deploying period = 0
%           2~n ChangeofPower when regulation deploying period = 1~n-1
%   RegDeployingPeriod_MultiplePeriod:1*n double
%           1 0
%           2~n RegDeployingPeriod     min
load('TrainConfiguration_example.mat');
TrainConfiguration = TrainConfiguration_example;
RegSet{1} = 1:50;
RegSet{2} = 60;
RegSet{4} = 1;
RegSet{5} = 36000;

RegDeployingPeriod_MultiplePeriod = 0:15;
for RegDeployingPeriod = RegDeployingPeriod_MultiplePeriod
    RegSet{3} = [0 RegDeployingPeriod];
    [Power,Energy,ChangeofPower] = DTP_evaluation(TrainConfiguration,RegSet);
    Power_MultiplePeriod(RegDeployingPeriod+1) = Power;
    Energy_MultiplePeriod(RegDeployingPeriod+1) = Energy;
    ChangeofPower_MultiplePeriod{RegDeployingPeriod+1} = ChangeofPower;
end

% Output:
%   EPC: ERES power capacity    MW
%   ESC: ERES storage capacity     MWh
%   RDR: Regulation deployable region:2*1500 double
%           1 Deploying period time point 1~15 min     
%           2 Deploy regulation power     MW
[EPC,ESC,RDR] = Aggregation_model(ChangeofPower_MultiplePeriod,RegDeployingPeriod_MultiplePeriod);

%% 3, Plot the result of 1 and 2
PlotResult(RDR,Power_MultiplePeriod,RegDeployingPeriod_MultiplePeriod);