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
kplmean = [ .15 ];       % s
kplstdd = [ .03 ];       % s
kvemean = [ 0.05 ];       % s
kvestdd = [ .01  ];       % s
t0mean  = [ 4    ];       % s
t0sttd  = [ 1.3  ] ;       % s
alphamean  =  [2.5];
alphasttd  =  [.3];
betamean  =  [4.5];
betasttd  =  [.3];
tisinput=[T1pmean; T1pstdd; T1lmean; T1lstdd; kplmean; kplstdd; kvemean; kvestdd;t0mean;t0sttd;alphamean; alphasttd; betamean ; betasttd ];

%% Variable Setup
Ntime = 30;
TR = 3;
TR_list = (0:(Ntime-1))*TR;
M0 = [0,0];
%ve = 0.95;
ve = 1.;
VIF_scale_fact = [100;0];
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
            for n = 1:Ntime
                %flips(2,n) = acos(sqrt((E1(2)^2-E1(2)^(2*(N-n+1)))/(1-E1(2)^(2*(N-n+1)))));
                flips(2,n) = 15*pi/180;
                flips(1,n) = 20*pi/180;
            end
            params.FaList = flips;
    end
    %% Fitting
    [t_axis,Mxy,Mz] = model.compile(M0.',params);
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
    NumberUncertain=3;
    switch (NumberUncertain)
       case(3)
         [x,xn,xm,w,wn]=GaussHermiteNDGauss(NGauss,[tisinput(5:2:9)],[tisinput(6:2:10)]);
       case(4)
         [x,xn,xm,w,wn]=GaussHermiteNDGauss(NGauss,[tisinput(1:2:7)],[tisinput(2:2:8)]);
    end 
    lqp=length(xn{1}(:));
    statevariable    = optimvar('state',Nspecies,Ntime,lqp,'LowerBound',0);
    auxvariable      = optimexpr(    [Nspecies,Ntime,lqp]);
    stateconstraint  = optimconstr(    [Nspecies,Ntime,lqp]);

    modelSNR = 10 ; % TODO - FIXME
    signuImage = (max(Mxy(1,:))+max(Mxy(2,:)))/2/modelSNR;
    % variance for Gauss RV is sum. sqrt for std
    signu = sqrt(2* Ntime) * signuImage;
    [x2,xn2,xm2,w2,wn2]=GaussHermiteNDGauss(NGauss,0,signu);
    lqp2=length(xn2{1}(:));


    disp('build state variable')
    stateconstraint(:,1,:)  = statevariable(:,1,:) ==0;
    auxvariable(:,1,:) =0;
    for iqp = 1:lqp
      for iii = 1:Ntime-1
        switch (NumberUncertain)
           case(3)
             T1Pqp   = T1pmean;
             T1Lqp   = T1lmean;
             kplqp   = xn{1}(iqp);
             klpqp   =    0 ;     % @cmwalker where do I get this from ? 
             kveqp   = xn{2}(iqp);
             t0qp    = xn{3}(iqp); 
           case(4)
             T1Pqp   = xn{1}(iqp);
             T1Lqp   = xn{2}(iqp);
             kplqp   = xn{3}(iqp);
             klpqp   =    0 ;     % @cmwalker where do I get this from ? 
             kveqp   = xn{4}(iqp);
             t0qp    = t0mean(1); 
        end 
        %
        currentTR = TR ;
        nsubstep = 5;
        deltat = currentTR /nsubstep ;
        % setup AIF
        integratedt = [TRList(iii):deltat:TRList(iii+1)] +deltat/2  ;
        integrand = jmA0 * gampdf(integratedt(1:nsubstep )'-t0qp,jmalpha,jmbeta) ;
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
        % linear approximation of aif    
        % >> syms tk tkm1 tau a gammak gammakm1
        % >> 
        % >> aifexpr = exp(a*(tk-tau)) *  (gammak * (tau - tkm1)/(tk-tkm1 ) + gammakm1 * (tk - tau )/(tk-tkm1 ) )  
        % >> aifint = int(aifexpr,tau,tkm1,tk)
        % aifint = -(gammak - gammakm1 - gammak*exp(a*tk - a*tkm1) + gammakm1*exp(a*tk - a*tkm1) + a*gammak*tk - a*gammak*tkm1 - a*gammakm1*tk*exp(a*tk - a*tkm1) + a*gammakm1*tkm1*exp(a*tk - a*tkm1))/(a^2*(tk - tkm1))

        expATR = [ exp(-currentTR*(kplqp + kveqp + 1/T1Pqp)),                   0; (kplqp*exp(-currentTR/T1Lqp) - kplqp*exp(-currentTR*(kplqp + kveqp + 1/T1Pqp)))/(kplqp + kveqp - 1/T1Lqp + 1/T1Pqp), exp(-currentTR/T1Lqp)];
        % mid-point rule integration
        aifterm = kveqp * deltat * [ exp((-1/T1Pqp - kplqp - kveqp)*deltat*[.5:1:nsubstep] );
    kplqp*(-exp((-1/T1Pqp - kplqp - kveqp)*deltat*[.5:1:nsubstep] ) + exp(-1/T1Lqp *deltat*[.5:1:nsubstep] ))/(1/T1Pqp + kplqp + kveqp - 1/T1Lqp )] * integrand ;
        auxvariable(:,iii+1,iqp) =  expATR *(cos(FaList(:,iii)).*auxvariable(:,iii,iqp ))   + aifterm ;
        stateconstraint(:,iii+1,iqp) = statevariable(:,iii+1,iqp) ==  expATR *(cos(FaList(:,iii)).*statevariable(:,iii,iqp ))   + aifterm ;
      end
    end

    disp('build objective function')
    sumstatevariable = optimexpr([Nspecies,lqp]);
    for jjj = 1:lqp
       sumstatevariable(:,jjj) =  sum(sin(FaList).*statevariable(:,:,jjj),2);
    end 
    %statematrix = optimexpr([lqp,lqp]);
    expandvar  = ones(1,lqp);
    diffsumm =(sumstatevariable(1,:)+sumstatevariable(2,:))' * expandvar   - expandvar' * (sumstatevariable(1,:)+sumstatevariable(2,:));
    Hz = 0;
    for jjj=1:lqp2
      znu=xn2{1}(jjj) ;
      Hz = Hz + wn2(jjj) * (wn(:)' * log(exp(-(znu + diffsumm).^2/2/signu^2 - log(signu) -log(2*pi)/2   ) * wn(:)));
      %Hz = Hz + wn2(jjj) * (wn(:)' * log(exp(-(znu + diffsumm).^2/2/signu^2                             ) * wn(:)));
    end
    %% MIGaussObj = Hz/sqrt(pi)^(NumberUncertain+1); 
    MIGaussObj = Hz;

    %% 
    % Create an optimization problem using these converted optimization expressions.
    
    disp('create optim prob')
    convprob = optimproblem('Objective',MIGaussObj , "Constraints",stateconstraint);
    myidx = varindex(convprob )
    %% 
    % View the new problem.
    
    %show(convprob)
    problem = prob2struct(convprob,'ObjectiveFunctionName','reducedObjective','ConstraintFunctionName','reducedConstraint');
    %% extraParamsobj = functions(problem.objective).workspace{1}.extraParams;
    %% extraParamscon = functions(problem.nonlcon).workspace{1}.extraParams;
    %% 
    % Solve the new problem. The solution is essentially the same as before.
    

    % truthconstraint = infeasibility(stateconstraint,x0);
    InitialGuess =  [flips(:)];   
    pmin =  [flips(:)*0];     
    pmax =  [flips(:)*0+35*pi/180];
    tolx=1.e-9;
    tolfun=5.e-4;
    maxiter=400;

    Fx = @(x) MIGHQuadHPTofts(x, problem, myidx,Nspecies,Ntime,auxvariable);
    %% debug info
    %% x0.FaList = params.FaList;
    %% x0.state  = evaluate(auxvariable ,x0);
    %% mystate = evaluate( sumstatevariable ,x0);
    %% Xfull = [ x0.FaList(:); x0.state(:)];
    %% [MIobjfun,initVals.g] = problem.objective(Xfull);
    %% [initConst.ineq,initConst.ceq, initConst.ineqGrad,initConst.ceqGrad] = problem.nonlcon(Xfull);
    %% [myobjfun, myobjfun_Der]= Fx(InitialGuess)
    [designopt,fval,exitflag,output,lambda,grad,hessian] ...
     =fmincon(Fx, InitialGuess ,[],[],[],[],pmin,pmax,[],...
        optimset('TolX',tolx,'TolFun',tolfun,'MaxIter', ...
        maxiter,'Display','iter-detailed',... 
        'GradObj','on','PlotFcn',{'optimplotfvalconstr', 'optimplotconstrviolation', 'optimplotfirstorderopt' }) ...
        );

    toc;
    handle = figure(5)
    optparams = params;
    optparams.FaList = reshape(designopt(:),size(params.FaList ));
    [t_axisopt,Mxyopt,Mzopt] = model.compile(M0.',params);
    figure(6)
    plot(optparams.TRList,Mxyopt(1,:),'b',optparams.TRList,Mxyopt(2,:),'k')
    ylabel('MI Mxy')
    xlabel('sec')
    figure(7)
    plot(optparams.TRList,optparams.FaList(1,:)*180/pi,'b',optparams.TRList,optparams.FaList(2,:)*180/pi,'k')
    ylabel('MI FA')
    xlabel('sec')
    handle = figure(8)
    plot(optparams.TRList,Mzopt(1,:),'b--',optparams.TRList,Mzopt(2,:),'k--')
    ylabel('MI Mz ')
    xlabel('sec'); legend('Pyr','Lac')
end 


function [MIobjfun, MIobjfun_Der]=MIGHQuadHPTofts(xopt,problem,myidx,Nspecies,Ntime,auxvariable)
    x0.FaList = reshape(xopt,Nspecies,Ntime);
    x0.state  = evaluate(auxvariable ,x0);
    Xfull = [ x0.FaList(:); x0.state(:)];
    [MIobjfun,initVals.g] = problem.objective(Xfull);
    [initConst.ineq,initConst.ceq, initConst.ineqGrad,initConst.ceqGrad] = problem.nonlcon(Xfull);
    objectiveGradFA    = initVals.g(myidx.FaList);
    objectiveGradState = initVals.g(myidx.state);
    jacobianFA    = initConst.ceqGrad(myidx.FaList,:);
    jacobianState = initConst.ceqGrad(myidx.state,:);
    adjointvar =-jacobianState \objectiveGradState ;
    MIobjfun_Der = objectiveGradFA +  jacobianFA *   adjointvar ;
end
