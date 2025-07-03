function [EPC,ESC,RDR] = Aggregation_model(ChangeofPower_MultiplePeriod,RegDeployingPeriod_MultiplePeriod)
%AGGREGATION_MODEL Aggregation model of ERES
% Input:
%   ChangeofPower_MultiplePeriod:1*n cell
%           1 ChangeofPower when regulation deploying period = 0
%           2~n ChangeofPower when regulation deploying period = 1~n-1
%   RegDeployingPeriod_MultiplePeriod:1*n double
%           1 0
%           2~n RegDeployingPeriod     min

% Output:
%   EPC: ERES power capacity    MW
%   ESC: ERES storage capacity     MWh
%   RDR: Regulation deployable region:2*1500 double
%           1 Deploying period time point 1~15 min     
%           2 Deploy regulation power     MW

tauNum = size(ChangeofPower_MultiplePeriod,2);
RegulationPower = [];
RegulationDuration = [];
for tauG_count = 1:tauNum
    tauG = RegDeployingPeriod_MultiplePeriod(tauG_count);
    ChangeofPower = ChangeofPower_MultiplePeriod{tauG_count};
    Duration = RegDeployingPeriod_MultiplePeriod(tauG_count)/60;
    if tauG == 0
        RegulationPower = [RegulationPower abs(ChangeofPower(:,1))];
        RegulationDuration = [RegulationDuration zeros(size(ChangeofPower(:,1)))];
    else
        RegulationPower = [RegulationPower abs(mean(ChangeofPower(:,1:tauG),2))];
        RegulationDuration = [RegulationDuration Duration*ones(size(RegulationPower(:,1)))];
    end
end
[Batterylist, ~] = BatteryEqu(RegulationDuration,abs(RegulationPower));
EPC = sum(Batterylist(1,:));
ESC = sum(Batterylist(2,:));
Tlist = 0:0.01:15;
for t = 1:length(Tlist)
    Powerlist(t) = TotalBatteryPowerFunc(Tlist(t),Batterylist);
end
RDR = [Tlist; Powerlist];
end

function [Batterylist, Errorlist] = BatteryEqu(RegulationDuration,RegulationPower)
if size(RegulationDuration,1) ~= size(RegulationPower,1) || size(RegulationDuration,2) ~= size(RegulationPower,2)
    error('Invalid input!')
end
BatteryNum = size(RegulationDuration,1);
Batterylist = zeros(2,BatteryNum);
Errorlist = zeros(1,BatteryNum);
for b = 1:BatteryNum
    if ~isempty(find(RegulationPower(b,1:2) < 0.00001, 1))
        continue
    end
    [P,E,~,~,error] = PieceFuncFit(RegulationDuration(b,:),RegulationPower(b,:),1);
    Batterylist(1,b) = P;
    Batterylist(2,b) = E;
    Errorlist(1,b) = error;
end
end

function [Coef_constant, Coef_inverse_a, Coef_inverse_b, BestSplitPoint, MSE_relative] = PieceFuncFit(t_data, y_data, MethodFlag)
Excrange = [-1, 1]/1e+04;
[t_data, y_data] = ExcFilter(t_data, y_data, Excrange);

switch MethodFlag
    case 1
        C = y_data(1);
        piecewise_error = @(tc) calculate_error_M1(tc, t_data, y_data, C);
        if isempty(find(y_data <= 0.6 * C, 1))
            t_c_max = t_data(end);
        else
            t_c_max = t_data(find(y_data <= 0.6 * C, 1));
        end
        tc_opt = fminbnd(piecewise_error, t_data(1), t_c_max);
        a_opt = C * tc_opt;
        fit_values = (t_data(2:end) < tc_opt) * C + (t_data(2:end) >= tc_opt) .* (a_opt ./ t_data(2:end)); % 计算拟合值
        fit_values = [C fit_values];
        residuals = y_data - fit_values; 
        MSE = mean(residuals.^2); 
        MSE_relative = MSE / mean(y_data); 
        Coef_constant = C;
        Coef_inverse_a = a_opt;
        Coef_inverse_b = 0;
        BestSplitPoint = tc_opt;
    case 2

        C = y_data(1);
        piecewise_error = @(params) calculate_error_M2(params, t_data, y_data, C);
        initial_params = [1, 1];
        opt_params = fminsearch(piecewise_error, initial_params);
        a_opt = opt_params(1);
        b_opt = opt_params(2);
        tc_opt = (a_opt / C) - b_opt;
        fit_values = (t_data < tc_opt) * C + (t_data >= tc_opt) .* (a_opt ./ (b_opt + t_data));
        residuals = y_data - fit_values;
        MSE = mean(residuals.^2);
        MSE_relative = MSE / mean(y_data);
        Coef_constant = C;
        Coef_inverse_a = a_opt;
        Coef_inverse_b = b_opt;
        BestSplitPoint = tc_opt;
end
end

function error = calculate_error_M1(tc, t_data, y_data, C)
a = C * tc;
fit_values = (t_data < tc) * C + (t_data >= tc) .* (a ./ t_data);
error = sum((fit_values(t_data >= tc) - y_data(t_data >= tc)).^2);
end

function error = calculate_error_M2(params, t_data, y_data, C)
a = params(1);
b = params(2);
tc = (a / C) - b;
fit_values = (t_data < tc) * C + (t_data >= tc) .* (a ./ (b + t_data));
error = sum((fit_values - y_data).^2);
end

function [t_data_fted, y_data_fted] = ExcFilter(t_data, y_data, Excrange)
Excidx = find(y_data >= Excrange(1) & y_data <= Excrange(2));
y_data(Excidx) = [];
t_data(Excidx) = [];
y_data_fted = y_data;
t_data_fted = t_data;
end

function [Power] = TotalBatteryPowerFunc(t,Batterylist)
Power = 0;
if isempty(Batterylist)
    return
end

Powerlist = Batterylist(1,:);
Energylist = Batterylist(2,:);
Durationlist = Energylist./Powerlist*60;

BatteryNum = size(Batterylist,2);
for bt = 1:BatteryNum
    if t < Durationlist(bt) || t == 0
        Power = Power + Powerlist(bt);
    else
        Power = Power + Energylist(bt)/t*60;
    end
end
end