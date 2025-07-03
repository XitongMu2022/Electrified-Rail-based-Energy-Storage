function [Power,Energy,ChangeofPower] = DEP_evaluation(TrainConfiguration,RegSet)
%DTP_EVALUATION Deployable ERES power evaluation
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
%           1 Start position        m
%           2 End position      m
%       Addmissiblex:1*2 double
%           1 Start time        sdelta_t
%           2 End time      s
%   RegSet:1*5 cell:
%       list_Trainindex: Index of trains participating in the ERES
%       delta_t: Discrete time      s
%       T_reg: Regulation deploying period     min
%       Strategy: 1,Upward; 2,Downward
%       t_begin: Regulation start time     s

% Output:
%   Power: Deployable ERES power capacity  MW
%   Energy: Deployable ERES energy  capacity  MWh
%   ChangeofPower: Power adjustment of every train  MW

%% LP problem data preparation

list_Trainindex = RegSet{1};
delta_t = RegSet{2};
T_reg = RegSet{3}; T_regb = T_reg(1); T_rege = T_reg(2); 
Strategy = RegSet{4};
t_begin = RegSet{5};


NumofTrain = length(list_Trainindex);
for i=1:NumofTrain
    ERESConfiguration_p(i).tau=[];
end

for pcount=1:NumofTrain
    i = list_Trainindex(pcount);
    TrainConfiguration_p(pcount).TrainPara = TrainConfiguration(i).TrainPara;
    TrainConfiguration_p(pcount).RoadPara = TrainConfiguration(i).RoadPara;
    TrainConfiguration_p(pcount).DepartureTime = TrainConfiguration(i).DepartureTime - t_begin;
    TrainConfiguration_p(pcount).GuidancetRef = TrainConfiguration(i).GuidancetRef - t_begin;
    TrainConfiguration_p(pcount).GuidanceVRef = TrainConfiguration(i).GuidanceVRef;
    TrainConfiguration_p(pcount).GuidanceXRef = TrainConfiguration(i).GuidanceXRef;
    TrainConfiguration_p(pcount).GuidanceuRef = TrainConfiguration(i).GuidanceuRef;
    TrainConfiguration_p(pcount).Addmissiblet = TrainConfiguration(i).Addmissiblet - t_begin;
    TrainConfiguration_p(pcount).Addmissiblex = TrainConfiguration(i).Addmissiblex;
end

for i=1:NumofTrain
    Addmissiblete=TrainConfiguration_p(i).Addmissiblet(2);
    T_mode = T_rege;
    if Addmissiblete >= T_mode*delta_t
        T_mode = floor(Addmissiblete/delta_t) + 1;
    end
end

for i=1:NumofTrain

    if isempty(TrainConfiguration_p(i).TrainPara)
        GSCConfiguration_p(i).tau=[];
        continue
    end

    veoutwithfullbraking=0;
    veoutwith0=0;
    vcru = [];

    TrainPara=TrainConfiguration_p(i).TrainPara;
    RoadPara=TrainConfiguration_p(i).RoadPara;
    DepartureTime=TrainConfiguration_p(i).DepartureTime;
    
    % parameters extraction
    Fm=TrainPara(1);
    Pm=TrainPara(2);
    gamma=TrainPara(3);
    a0=TrainPara(4);
    a1=TrainPara(5);
    a2=TrainPara(6);
    Bm=TrainPara(7);
    m=TrainPara(8);
    etaT=TrainPara(9);
    etaB=TrainPara(10);
    vturn=Pm/Fm;
    
    tRefi=TrainConfiguration_p(i).GuidancetRef;
    VRefi=TrainConfiguration_p(i).GuidanceVRef;
    XRefi=TrainConfiguration_p(i).GuidanceXRef;
    uRefi=TrainConfiguration_p(i).GuidanceuRef;
    
    if length(TrainConfiguration_p(i).Addmissiblet)==0
        GSCConfiguration_p(i).tau=[];
    else
        Addmissibletb=TrainConfiguration_p(i).Addmissiblet(1);
        Addmissiblete=TrainConfiguration_p(i).Addmissiblet(2);
        vcru = VRefi(find(XRefi==TrainConfiguration_p(i).Addmissiblex(1),1)+1);
        if isempty(vcru)
            continue
        end
        
        for btemp=0:T_mode
            if Addmissibletb<btemp*delta_t
                break;
            end
        end
        for etemp=T_mode:-1:0
            if Addmissiblete>etemp*delta_t
                break;
            end
        end
        if etemp>btemp
            GSCConfiguration_p(i).tau=[btemp,etemp];
        else
            GSCConfiguration_p(i).tau=[];
        end 
    end
    
    if isempty(GSCConfiguration_p(i).tau)
        GSCConfiguration_p(i).MechEnergy=[]; % [Mec traction, Mec braking ]J/kg
    else
        taub=GSCConfiguration_p(i).tau(1);
        taue=GSCConfiguration_p(i).tau(2);
        OriginalEnergyM=[];        
        for j=taub:taue-1
            ubegin=hInterpolation(tRefi,uRefi,j*delta_t);
            vbegin=hInterpolation(tRefi,VRefi,j*delta_t);
            PMTbegin=max(ubegin,0)*Pm;
            PMBbegin=-min(ubegin,0)*Bm*vbegin;
            ttemp=[j*delta_t];
            PMTtemp=[PMTbegin];
            PMBtemp=[PMBbegin];
            for itempbegin=1:length(tRefi)
                if tRefi(itempbegin)>j*delta_t
                    break;
                end
            end
            for itempend=1:length(tRefi)
                if tRefi(itempend)>(j+1)*delta_t
                    break;
                end
            end
            for itemp=itempbegin:itempend-1
                ttemp=[ttemp,tRefi(itemp)];
                PMTtemp=[PMTtemp,max(uRefi(itemp),0)*Pm];
                PMBtemp=[PMBtemp,-min(uRefi(itemp),0)*VRefi(itemp)*Bm];
            end
            uend=hInterpolation(tRefi,uRefi,(j+1)*delta_t);
            vend=hInterpolation(tRefi,VRefi,(j+1)*delta_t);
            PMTend=max(uend,0)*Pm;
            PMBend=-min(uend,0)*Bm*vend;
            ttemp=[ttemp,(j+1)*delta_t];
            PMTtemp=[PMTtemp,PMTend];
            PMBtemp=[PMBtemp,PMBend];
            OriginalEnergyM=[OriginalEnergyM;trapz(ttemp,PMTtemp),trapz(ttemp,PMBtemp)]; % [Mec traction, Mec braking(>=0) ]J/kg
        end
        GSCConfiguration_p(i).MechEnergy=OriginalEnergyM; % [Mec traction, Mec braking ]J/kg
    end
    
    % Constraints
    if isempty(GSCConfiguration_p(i).tau)
        GSCConfiguration_p(i).MecPowerBound=[];  % W/kg
    else
        taub=GSCConfiguration_p(i).tau(1);
        taue=GSCConfiguration_p(i).tau(2);
        vbegin=hInterpolation(tRefi,VRefi,taub*delta_t);
        xbegin=hInterpolation(tRefi,XRefi,taub*delta_t);
        MecEnergyBound=[]; % J/kg
        for j=taub:taue-1
            xend=hInterpolation(tRefi,XRefi,(j+1)*delta_t);
            [SecNum,~]=size(RoadPara);
            for itempend=1:SecNum
                if xend<=sum(RoadPara(1:itempend,2))
                    break;
                end
            end
            vlim=RoadPara(itempend,1);
            MecEnergyBound2 = 0.5 * (1+gamma) * (vlim^2 - vcru^2);
            MecEnergyBound3 = 0.5 * (1+gamma) * (0 - vcru^2);
            MecEnergyBound=[MecEnergyBound;MecEnergyBound2,MecEnergyBound3];

        end     
        GSCConfiguration_p(i).MecEnergyBound=MecEnergyBound; % >0 for traction; <0 for braking; W/kg
    end

end

if T_rege == 0

jmax = 1;
yalmip('clear');
z=sdpvar(1);
DPM=sdpvar(NumofTrain,T_mode,'full');

Constraints=[];

coefficientDPE=zeros(NumofTrain,T_mode);
for i=1:NumofTrain
    if isempty(GSCConfiguration_p(i).tau)
        Constraints=[Constraints,DPM(i,:)==0];
    else
        
        taub=GSCConfiguration_p(i).tau(1);
        taue=GSCConfiguration_p(i).tau(2);
        MechEnergy=GSCConfiguration_p(i).MechEnergy;
        MecEnergyBound=GSCConfiguration_p(i).MecEnergyBound;

        if MecEnergyBound(1,1)<=0 || MecEnergyBound(1,2)>=0 
            Constraints=[Constraints,DPM(i,:)==0];
            continue
        end
        
        TrainPara=TrainConfiguration_p(i).TrainPara;
        m=TrainPara(8);
        etaT=TrainPara(9);
        etaB=TrainPara(10);
        Pm=TrainPara(2);       
        
        for j=1:T_mode
            if j<=taub
                Constraints=[Constraints,DPM(i,j)==0];
            elseif j>taue
                Constraints=[Constraints,DPM(i,j)==0];
            else
                TractionEnergy=MechEnergy(j-taub,1);
                BrakingEnergy=MechEnergy(j-taub,2);

                if BrakingEnergy ~= 0 
                    Constraints=[Constraints,DPM(i,j)==0];
                    continue
                end

                Constraints=[Constraints,delta_t*sum(DPM(i,taub+1:j))<=MecEnergyBound(j-taub,1),delta_t*sum(DPM(i,taub+1:j))>=MecEnergyBound(j-taub,2)];
                
                PMij=(TractionEnergy-BrakingEnergy)/delta_t;

                if PMij<=-Pm || PMij>=Pm 
                    Constraints=[Constraints,DPM(i,:)==0];
                    continue
                end

                if TractionEnergy>5*BrakingEnergy || BrakingEnergy<0.01 % traction dominates
                    if Strategy == 1 || Strategy == 2
                        Constraints=[Constraints,-Pm<=DPM(i,j)+PMij,DPM(i,j)+PMij<=Pm]; 
                    else
                        error('Invalid RegSet!')
                    end
                    coefficientDPE(i,j)=m/etaT;
                elseif 5*TractionEnergy<BrakingEnergy % braking dominates
                    Constraints=[Constraints,-Pm<=DPM(i,j)+PMij,DPM(i,j)+PMij<=Pm]; 
                    coefficientDPE(i,j)=m*etaB;
                else
                    Constraints=[Constraints,DPM(i,j)==0];
                end
            end
        end
        TotalElecEnergy=sum(sum(MechEnergy));
        Constraints=[Constraints,sum(DPM(i,:))==0.00*TotalElecEnergy/delta_t];
    end
end
coefficientDPE=coefficientDPE/10^6;
AuxiliaryConstraints=[];
options = sdpsettings('verbose',0,'solver','gurobi');
AuxiliaryConstraints=[];
for j = 1:jmax
    if Strategy ==1
        AuxiliaryConstraints=[AuxiliaryConstraints,z>=coefficientDPE(:,j)'*DPM(:,j)];
    else
        AuxiliaryConstraints=[AuxiliaryConstraints,z<=coefficientDPE(:,j)'*DPM(:,j)];
    end
end

% Solve the problem
if Strategy ==1
    sol = optimize([Constraints,AuxiliaryConstraints],z,options);
else
    sol = optimize([Constraints,AuxiliaryConstraints],-z,options);
end
if sol.problem == 0
    value_DPM=value(DPM);
    ChangeofPower_Ref=coefficientDPE.*value_DPM;
else
    disp('Hmm, something went wrong!');
    sol.info
end

N_index = find(ChangeofPower_Ref(:,1));
coefficientDPE=zeros(NumofTrain,1);
DPM = zeros(NumofTrain,1);
for i = 1:NumofTrain
    if isempty(find(i == N_index,1))
        continue
    else       
        taub=GSCConfiguration_p(i).tau(1);
        taue=GSCConfiguration_p(i).tau(2);
        MechEnergy=GSCConfiguration_p(i).MechEnergy;
        MecEnergyBound=GSCConfiguration_p(i).MecEnergyBound;
        if MecEnergyBound(1,1)<=0 || MecEnergyBound(1,2)>=0 
            continue
        end
        
        TrainPara=TrainConfiguration_p(i).TrainPara;
        m=TrainPara(8);
        etaT=TrainPara(9);
        etaB=TrainPara(10);
        Pm=TrainPara(2);

        TractionEnergy=MechEnergy(1,1);
        BrakingEnergy=MechEnergy(1,2);

        PMij=(TractionEnergy-BrakingEnergy)/delta_t;
        if PMij<=-Pm || PMij>=Pm
            continue
        end

        switch Strategy
            case 2
        DPM(i,:) = Pm - PMij;
            case 1
        DPM(i,:) = - Pm - PMij;
            otherwise
        error('Invalid RegSet!')
        end
        coefficientDPE(i,:)=m/etaT;
    end
    
end
    ChangeofPower=coefficientDPE.*DPM/1e+06;
    Power = abs(sum(ChangeofPower));
    Energy = 0;

elseif T_rege >= 1
jmax = T_rege;
yalmip('clear');
z=sdpvar(1);
DPM=sdpvar(NumofTrain,T_mode,'full');

Constraints=[];

coefficientDPE=zeros(NumofTrain,T_mode);
for i=1:NumofTrain
    if isempty(GSCConfiguration_p(i).tau)
        Constraints=[Constraints,DPM(i,:)==0];
    else
        
        taub=GSCConfiguration_p(i).tau(1);
        taue=GSCConfiguration_p(i).tau(2);
        MechEnergy=GSCConfiguration_p(i).MechEnergy;
        MecEnergyBound=GSCConfiguration_p(i).MecEnergyBound;

        if MecEnergyBound(1,1)<=0 || MecEnergyBound(1,2)>=0 
            Constraints=[Constraints,DPM(i,:)==0];
            continue
        end
        
        TrainPara=TrainConfiguration_p(i).TrainPara;
        m=TrainPara(8);
        etaT=TrainPara(9);
        etaB=TrainPara(10);
        Pm=TrainPara(2);       
        
        for j=1:T_mode
            if j<=taub
                Constraints=[Constraints,DPM(i,j)==0];
            elseif j>taue
                Constraints=[Constraints,DPM(i,j)==0];
            else
                TractionEnergy=MechEnergy(j-taub,1);
                BrakingEnergy=MechEnergy(j-taub,2);

                if BrakingEnergy ~= 0 
                    Constraints=[Constraints,DPM(i,j)==0];
                    continue
                end

                Constraints=[Constraints,delta_t*sum(DPM(i,taub+1:j))<=MecEnergyBound(j-taub,1),delta_t*sum(DPM(i,taub+1:j))>=MecEnergyBound(j-taub,2)];
                
                PMij=(TractionEnergy-BrakingEnergy)/delta_t;

                if PMij<=-Pm || PMij>=Pm 
                    Constraints=[Constraints,DPM(i,:)==0];
                    continue
                end

                if TractionEnergy>5*BrakingEnergy || BrakingEnergy<0.01 % traction dominates
                    if Strategy == 1 || Strategy == 2
                        Constraints=[Constraints,-Pm<=DPM(i,j)+PMij,DPM(i,j)+PMij<=Pm]; 
                    else
                        error('Invalid RegSet!')
                    end
                    coefficientDPE(i,j)=m/etaT;
                elseif 5*TractionEnergy<BrakingEnergy % braking dominates
                    Constraints=[Constraints,-Pm<=DPM(i,j)+PMij,DPM(i,j)+PMij<=Pm]; 
                    coefficientDPE(i,j)=m*etaB;
                else
                    Constraints=[Constraints,DPM(i,j)==0];
                end
            end
        end
        TotalElecEnergy=sum(sum(MechEnergy));
        Constraints=[Constraints,sum(DPM(i,:))==0.00*TotalElecEnergy/delta_t];
    end
end
coefficientDPE=coefficientDPE/10^6;
AuxiliaryConstraints=[];
options = sdpsettings('verbose',0,'solver','gurobi');
AuxiliaryConstraints=[];
for j = 1:jmax
    if Strategy ==1
        AuxiliaryConstraints=[AuxiliaryConstraints,z>=coefficientDPE(:,j)'*DPM(:,j)];
    else
        AuxiliaryConstraints=[AuxiliaryConstraints,z<=coefficientDPE(:,j)'*DPM(:,j)];
    end
end

% Solve the problem
if Strategy ==1
    sol = optimize([Constraints,AuxiliaryConstraints],z,options);
else
    sol = optimize([Constraints,AuxiliaryConstraints],-z,options);
end
if sol.problem == 0
    Power = abs(value(z));
    Energy = Power*jmax*delta_t/3600;
    value_DPM=value(DPM);
    ChangeofPower=coefficientDPE.*value_DPM;
else
    disp('Hmm, something went wrong!');
    sol.info
end

else
    error('Invalid RegSet!')
end

end

function b=hInterpolation(A,B,a)
if length(A)<2
    b=[];
    return;
end

for locindex=1:length(A)
    if A(locindex)>a
        break;
    end
end
if locindex==1
    b=B(1);
elseif locindex==length(A)
    b=B(locindex);
else
    b=interp1([A(locindex-1),A(locindex)]',[B(locindex-1),B(locindex)]',a);
end
end