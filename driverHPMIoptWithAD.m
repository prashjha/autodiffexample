% ShowMxyPub Shows the flip angle schedules and the resulting magnetization
% evolutions for the paper
% Author: Chris Walker
% Date: 8/6/2018

clear all
close all
clc

%% Tissue Parameters
T1pmean = [ 30 ]; % s
T1pstdd = [ 10 ]; % s
T1lmean = [ 25 ]; % s
T1lstdd = [ 10 ]; % s
kplmean = [ 5 ];       % s
kplstdd = [ 5 ];       % s
kvemean = [ 0.15 ];       % s
kvestdd = [ .05  ];       % s
t0mean  = [ 4    ];       % s
t0sttd  = [ 1.3  ] ;       % s
alphamean  =  [2.5];
betamean  =  [4.5];
tisinput=[T1pmean; T1pstdd; T1lmean; T1lstdd; kplmean; kplstdd; kvemean; kvestdd;t0mean;t0sttd];

%% Variable Setup
Ntime = 23;
TR = 2;
TR_list = (0:(Ntime-1))*TR;
M0 = [0,0];
%ve = 0.95;
ve = 1.;
VIF_scale_fact = [1;0];
bb_flip_angle = 20;
opts = optimset('lsqcurvefit');
opts.TolFun = 1e-09;
opts.TolX = 1e-09;
opts.Display = 'off';
params = struct('t0',[t0mean(1);0],'gammaPdfA',[alphamean(1)  ;1],'gammaPdfB',[betamean(1);1],...
    'scaleFactor',VIF_scale_fact,'T1s',[T1pmean(1),T1lmean(1)],'ExchangeTerms',[0,kplmean(1) ;0,0],...
    'TRList',TR_list,'PerfusionTerms',[kvemean(1),0],'volumeFractions',ve,...
    'fitOptions', opts);
model = HPKinetics.NewMultiPoolTofftsGammaVIF();
  

%% Get true Mz
%% Choose Excitation Angle
FAType = {'Const'};
%% HACK- @cmwalker code for initial conditions - https://github.com/fuentesdt/TumorHPMRI/blob/master/models/gPC/walker/ShowMxyPub.m
for i = 1:numel(FAType)
    switch (FAType{i})
        case('Const') % Nagashima for lactate const 10 pyruvate
            tic
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
    jmA0    = VIF_scale_fact(1);
    jmalpha = alphamean(1);
    jmbeta  = betamean(1);
    jmt0    = t0mean(1);
    jmaif   = jmA0  * gampdf(TR_list - jmt0  , jmalpha , jmbeta);
    figure(4)
    plot(TR_list,jmaif ,'b')
    ylabel('aif')
    xlabel('sec')
end


%% compute TR's from TR_list
TRi = getTR(TR_list); % vector of 'size(TR_list) - 1'

%% optimize MI for TR and FA
optf = true;
if optf
    tic;
    % Convert this function file to an optimization expression.
    
    %% 
    % Furthermore, you can also convert the |rosenbrock| function handle, which 
    % was defined at the beginning of the plotting routine, into an optimization expression.
    

    % setup optimization variables
    Nspecies = 2
    FaList = optimvar('FaList',Nspecies,Ntime,'LowerBound',0, 'UpperBound',35*pi/180);
    TRList = TR_list;
    NGauss  = 3
    NumberUncertain=4;
    [x,xn,xm,w,wn]=GaussHermiteNDGauss(NGauss,[tisinput(1:2:7)],[tisinput(2:2:8)]);
    lqp=length(xn{1}(:));
    statevariable = optimvar('state',Nspecies,Ntime,lqp,'LowerBound',0);
    stateconstraint  = optimconstr(    [Nspecies,Ntime,lqp]);
    % scaling important for the optimization step length update
    scalestate = 1.e-2;

    modelSNR = 10 ; % TODO - FIXME
    signuImage = (max(Mxy(1,:))+max(Mxy(2,:)))/2/modelSNR;
    signu = Ntime * signuImage;
    [x2,xn2,xm2,w2,wn2]=GaussHermiteNDGauss(NGauss,0,signu);
    lqp2=length(xn2{1}(:));


    disp('build state variable')
    stateconstraint(:,1,:)  = statevariable(:,1,:) ==0;
    for iqp = 1:lqp
      for iii = 1:Ntime-1
        % disp([Nspecies*Ntime*(iqp-1) + Nspecies*(iii-1)+1 , Nspecies*Ntime*(iqp-1) + Nspecies*(iii-1)+2 ,Ntime* Nspecies* lqp])
        T1Pqp = xn{1}(iqp);
        T1Lqp = xn{2}(iqp);
        kplqp = xn{3}(iqp);
        klpqp =    0 ;     % @cmwalker where do I get this from ? 
        kveqp = xn{4}(iqp);
        %
        currentTR = TR ;
        nsubstep = 5;
        deltat = currentTR /nsubstep ;
        % setup AIF
        integratedt =TRList(iii)+ [1:2:2*nsubstep+1]*deltat/2;
        integrand = jmA0 * gampdf(integratedt(1:nsubstep )'-jmt0,jmalpha,jmbeta) ;
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
        expATR = [ exp(-currentTR*(kplqp + kveqp + 1/T1Pqp)),                   0; (kplqp*exp(-currentTR/T1Lqp) - kplqp*exp(-currentTR*(kplqp + kveqp + 1/T1Pqp)))/(kplqp + kveqp - 1/T1Lqp + 1/T1Pqp), exp(-currentTR/T1Lqp)];
        % mid-point rule integration
        aifterm = kveqp * deltat * [ exp((-1/T1Pqp - kplqp - kveqp)*deltat*[.5:1:nsubstep] );
    kplqp*(-exp((-1/T1Pqp - kplqp - kveqp)*deltat*[.5:1:nsubstep] ) + exp(-1/T1Lqp *deltat*[.5:1:nsubstep] ))/(1/T1Pqp + kplqp + kveqp - 1/T1Lqp )] * integrand ;
        stateconstraint(:,iii+1,iqp) = scalestate*statevariable(:,iii+1,iqp) ==  expATR *(scalestate*cos(FaList(:,iii)).*statevariable(:,iii,iqp ))   + aifterm ;
      end
    end

    disp('build objective function')
    sumstatevariable = optimexpr([Nspecies,lqp]);
    for jjj = 1:lqp
       sumstatevariable(:,jjj) =  sum(scalestate*sin(FaList).*statevariable(:,:,jjj),2);
    end 
    %statematrix = optimexpr([lqp,lqp]);
    expandvar  = ones(1,lqp);
    diffsummone = sumstatevariable(1,:)' * expandvar   - expandvar' * sumstatevariable(1,:);
    diffsummtwo = sumstatevariable(2,:)' * expandvar   - expandvar' * sumstatevariable(2,:);
    Hz = 0;
    for jjj=1:lqp2
      znu=xn2{1}(jjj) ;
      Hz = Hz + wn2(jjj) * (wn(:)' * log(exp(-(znu + diffsummone).^2/2/signu^2   - (znu + diffsummtwo).^2/2/signu^2  ) * wn(:)));
    end
    MIGaussObj = Hz/sqrt(pi)^(NumberUncertain+1); 

    %% 
    % Create an optimization problem using these converted optimization expressions.
    
    disp('create optim prob')
    convprob = optimproblem('Objective',MIGaussObj , "Constraints",stateconstraint);
    %% 
    % View the new problem.
    
    %show(convprob)
    %problem = prob2struct(convprob,'ObjectiveFunctionName','generatedObjective');
    %% 
    % Solve the new problem. The solution is essentially the same as before.
    
    x0.FaList = params.FaList;
    x0.state  = repmat( 1/scalestate * Mz./cos(params.FaList),1,1,lqp);
    % truthconstraint = infeasibility(stateconstraint,x0);
    myoptions = optimoptions(@fmincon,'Display','iter-detailed','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'MaxFunctionEvaluations',1e7,'ConstraintTolerance',2.e-9, 'OptimalityTolerance',2.5e-9,'Algorithm','interior-point','StepTolerance',1.000000e-12,'MaxIterations',1000,'PlotFcn',{'optimplotfvalconstr', 'optimplotconstrviolation', 'optimplotfirstorderopt' },'HonorBounds',true, 'HessianApproximation', 'lbfgs' ,'Diagnostic','on','FunValCheck','on' )
    [popt,fval,exitflag,output] = solve(convprob,x0,'Options',myoptions, 'ObjectiveDerivative', 'auto-reverse' , 'ConstraintDerivative', 'auto-reverse')
    %[popt,fval,exitflag,output] = solve(convprob,x0 )


    toc;
    handle = figure(5)

    optparams = params;
    optparams.FaList = popt.FaList;
    [t_axisopt,Mxyopt,Mzopt] = model.compile(M0.',params);
    figure(6)
    plot(optparams.TRList,Mxyopt(1,:),'b',optparams.TRList,Mxyopt(2,:),'k')
    ylabel('MI Mxy')
    xlabel('sec')
    figure(7)
    plot(optparams.TRList,optparams.FaList(1,:),'b',optparams.TRList,optparams.FaList(2,:),'k')
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

