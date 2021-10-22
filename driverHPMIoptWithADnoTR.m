% ShowMxyPub Shows the flip angle schedules and the resulting magnetization
% evolutions for the paper
% Author: Chris Walker
% Date: 8/6/2018

clear all
close all
clc
N_vect = round(linspace(5,60,15));
%% Variable Setup
Ntime = 5;
TR = 2;
TR_list = (0:(Ntime-1))*TR;
T1a = 43;
T1b = 33;
Kpl = 0.1;
alpha = 2.5;
beta = 4.5;
M0 = [0,0];
kve = 0.02;
ve = 0.95;
VIF_scale_fact = [1;0];
bb_flip_angle = 20;
opts = optimset('lsqcurvefit');
opts.TolFun = 1e-09;
opts.TolX = 1e-09;
opts.Display = 'off';
params = struct('t0',[0;0],'gammaPdfA',[alpha;1],'gammaPdfB',[beta;1],...
    'scaleFactor',VIF_scale_fact,'T1s',[T1a,T1b],'ExchangeTerms',[0,Kpl;0,0],...
    'TRList',TR_list,'PerfusionTerms',[kve,0],'volumeFractions',ve,...
    'fitOptions', opts);
model = HPKinetics.NewMultiPoolTofftsGammaVIF();

%% Tissue Parameters
T1pmean = [ 30 ]; % s
T1pstdd = [ 10 ]; % s
T1lmean = [ 25 ]; % s
T1lstdd = [ 10 ]; % s
kplmean = [ 5 ];       % s
kplstdd = [ 5 ];       % s
kvemean = [ 0.15 ];       % s
kvestdd = [ .05  ];       % s
tisinput=[T1pmean; T1pstdd; T1lmean; T1lstdd; kplmean; kplstdd; kvemean; kvestdd];

%% Get true Mz
%% Choose Excitation Angle
FAType = {'Const'};
%% HACK- @cmwalker code for initial conditions - https://github.com/fuentesdt/TumorHPMRI/blob/master/models/gPC/walker/ShowMxyPub.m
for i = 1:numel(FAType)
    switch (FAType{i})
        case('Const') % Nagashima for lactate const 10 pyruvate
            tic
            E1(1) = exp(-TR*(1/T1a+Kpl));
            E1(2) = exp(-TR/T1b);
            for n = 1:Ntime
                %flips(2,n) = acos(sqrt((E1(2)^2-E1(2)^(2*(N-n+1)))/(1-E1(2)^(2*(N-n+1)))));
                flips(2,n) = 15*pi/180;
                flips(1,n) = 20*pi/180;
            end
            params.FaList = flips;
    end
    tic
    %% Fitting
    [t_axis,Mxy,Mz] = model.compile(M0.',params);
    toc
    save_Mxy{i} = Mxy;
    save_Mz{i} = Mz;
    save_t_axis{i} = t_axis;
    save_flip_angles{i} = params.FaList;
end

%% Plot initial guess
plotinit = true;
if plotinit
    % plot initial guess
    figure(1)
    plot(TR_list,Mxy(1,:),'b',TR_list,Mxy(2,:),'k')
    ylabel('Const Mxy')
    xlabel('sec')
    figure(2)
    plot(TR_list,flips(1,:),'b',TR_list,flips(2,:),'k')
    ylabel('Const FA')
    xlabel('sec')
    % save('tmpShowMxyPub')
    figure(3)
    plot(TR_list,Mz(1,:),'b',TR_list,Mz(2,:),'k')
    ylabel('Const Mz')
    xlabel('sec')

    % plot gamma
    jmA0    = 10.
    jmalpha = 2.5
    jmbeta  = 4.5
    jmt0    = 0
    jmaif   = jmA0  * gampdf(TR_list - jmt0  , jmalpha , jmbeta);
    figure(20)
    plot(TR_list,jmaif ,'b')
    ylabel('aif')
    xlabel('sec')
end


%% compute TR's from TR_list
TRi = getTR(TR_list); % vector of 'size(TR_list) - 1'

%% optimize MI for TR and FA
optf = true;
if optf
    % Pulse Sequence Bounds
    pmin =  [TRi'-1.5; flips(:)*0];     % <-- constraints on TR's and not TR_list
    pmax =  [TRi'+1.5; flips(:)*0+pi/2];% <-- constraints on TR's and not TR_list
    findiffrelstep=1.e-6;
    tolx=1.e-9;%1.e-5;
    tolfun=1.e-9;%1.e-5;QALAS_synphan_MIcalc.m
    maxiter=500;

    tic;
    % Convert this function file to an optimization expression.
    
    %% 
    % Furthermore, you can also convert the |rosenbrock| function handle, which 
    % was defined at the beginning of the plotting routine, into an optimization expression.
    

    % setup optimization variables
    Nspecies = 2
    FaList = optimvar('FaList',2,Ntime);
    TRList = TR_list;
    diffTR = diff(TRList);
    NGauss  = 2
    [x,xn,xm,w,wn]=GaussHermiteNDGauss(NGauss,[tisinput(1:2:7)],[tisinput(2:2:8)]);
    lqp=length(xn{1}(:));
    statevariable  = optimexpr([Nspecies,Ntime,lqp]);

    signu = 10 ; % TODO - FIXME
    [x2,xn2,xm2,w2,wn2]=GaussHermiteNDGauss(NGauss,0,signu);
    lqp2=length(xn2{1}(:));


    disp('build state variable')
    for iqp = 1:lqp
      T1Pqp = xn{1}(iqp);
      T1Lqp = xn{2}(iqp);
      kplqp = xn{3}(iqp);
      klpqp =    0 ;     % @cmwalker where do I get this from ? 
      kveqp = xn{4}(iqp);
      currentTR = 2;
      % >> syms a  kpl d currentTR    T1P kveqp T1L 
      % >> expATR = expm([a,  0; kpl, d ] * currentTR )
      % 
      % expATR =
      % 
      % [                                     exp(a*currentTR),                0]
      % [(kpl*exp(a*currentTR) - kpl*exp(currentTR*d))/(a - d), exp(currentTR*d)]
      % 
      % >> a = -1/T1P - kpl - kveqp
      % >> d = -1/T1L
      % >> eval(expATR)
      % 
      % ans =
      % 
      % [                                                              exp(-currentTR*(kpl + kveqp + 1/T1P)),                   0]
      % [(kpl*exp(-currentTR/T1L) - kpl*exp(-currentTR*(kpl + kveqp + 1/T1P)))/(kpl + kveqp - 1/T1L + 1/T1P), exp(-currentTR/T1L)]
      %    
      %expATR = fcn2optimexpr(@expm,A*currentTR );
      % A = [-1/T1P - kpl - kveqp,  0; kpl, -1/T1L ];
      expATR = [ exp(-currentTR*(kplqp + kveqp + 1/T1Pqp)),                   0; (kplqp*exp(-currentTR/T1Lqp) - kplqp*exp(-currentTR*(kplqp + kveqp + 1/T1Pqp)))/(kplqp + kveqp - 1/T1Lqp + 1/T1Pqp), exp(-currentTR/T1Lqp)];
      %pyruvatematrix = diag(ones(Ntime,1)) + diag( - expATR(1,1)*  cos(FaList(1,1:(Ntime-1))),-1);
      pyruvateRHS = rand(Ntime,1) ;
      %statevariable(1,:,iqp) = pyruvatematrix \pyruvateRHS ;
      for iii = 2:Ntime
        statevariable(1,iii,iqp) = pyruvateRHS(iii)  + expATR(1,1)* cos(FaList(1,iii-1)) * statevariable(1,iii-1,iqp);
      end 
      %lactatematrix = diag(ones(Ntime,1)) + diag( - expATR(2,2)*  cos(FaList(2,1:(Ntime-1))),-1);
      lactateRHS = cos(FaList(2,:)).* statevariable(1,:,iqp);
      %statevariable(2,:,iqp) = lactatematrix \lactateRHS ;
       
      for iii = 2:Ntime
        statevariable(2,iii,iqp) = lactateRHS(iii)  + expATR(2,2)* cos(FaList(2,iii-1)) * statevariable(2,iii-1,iqp);
      end 
    end

    disp('build objective function')
    sumstatevariable = squeeze(sum(statevariable,2));
    %statematrix = optimexpr([lqp,lqp]);
    diffsummone = repmat(sumstatevariable(1,:)',1,lqp) - repmat(sumstatevariable(1,:) ,lqp,1);
    diffsummtwo = repmat(sumstatevariable(2,:)',1,lqp) - repmat(sumstatevariable(2,:) ,lqp,1);
    Hz = 0;
    for jjj=1:lqp2
      znu=xn2{1}(jjj) ;
      Hz = Hz + wn2(jjj) * (wn(:)' * log(exp(-(znu + diffsummone).^2/sqrt(2)/signu   - (znu + diffsummtwo).^2/sqrt(2)/signu  ) * wn(:)));
    end
    %% Hz = 0;
    %% for iii=1:lqp
    %%     for jjj=1:lqp2
    %%         znu=xn2{1}(jjj) ;
    %%         lntermtmp=0;
    %%         for kkk=1:lqp
    %%             lntermtmp=lntermtmp + wn(kkk) * exp(-(znu+sum(statevariable(1,:,iii))- sum(statevariable(1,:,kkk)))^2/sqrt(2)/signu);
    %%             lntermtmp=lntermtmp + wn(kkk) * exp(-(znu+sum(statevariable(2,:,iii))- sum(statevariable(2,:,kkk)))^2/sqrt(2)/signu);
    %%         end
    %%         % for function
    %%         lnterm = log(lntermtmp)+log(pi^(-1.5));
    %%         Hz = Hz + wn(iii) * wn2(jjj) * lnterm;
    %%     end
    %% end
    MIGaussObj = -pi^(-1.5-2.5)*Hz; 

    %% 
    % Create an optimization problem using these converted optimization expressions.
    
    disp('create optim prob')
    convprob = optimproblem('Objective',MIGaussObj );
    %% 
    % View the new problem.
    
    %show(convprob)
    problem = prob2struct(convprob,'ObjectiveFunctionName','generatedObjective');
    %% 
    % Solve the new problem. The solution is essentially the same as before.
    
    x0.FaList = params.FaList;
    myoptions = optimoptions(@fminunc,'Display','iter-detailed','SpecifyObjectiveGradient',true)
    [popt,fval,exitflag,output] = solve(convprob,x0,'Options',myoptions, 'ObjectiveDerivative', 'auto-reverse' )
    %[popt,fval,exitflag,output] = solve(convprob,x0 )


    toc;

    params.FaList = popt.FaList;
    [t_axisopt,Mxyopt,Mzopt] = model.compile(M0.',params);
    figure(4)
    plot(params.TRList,Mxyopt(1,:),'b',params.TRList,Mxyopt(2,:),'k')
    ylabel('MI Mxy')
    xlabel('sec')
    figure(5)
    plot(params.TRList,params.FaList(1,:),'b',params.TRList,params.FaList(2,:),'k')
    ylabel('MI FA')
    xlabel('sec')
end 


%% convert time sequence to TR and TR to time sequence
function TR = getTR(t)
% compute TR from time sequence
    N = size(t,2);
    TR = zeros(1,N-1);
    for i=1:(N-1)
        TR(i) = t(i+1) - t(i);
    end
end

