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
lnsr2pi = 0.9189385332046727; % log(sqrt(2*pi))
Ntime = 30;
TR = 3;
TR_list = (0:(Ntime-1))*TR;
nsubstep = 3;
SubTimeList = (0:(Ntime-1)*nsubstep )*TR/nsubstep;
M0 = [0,0];
ve = 0.95;
%ve = 1.;
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
    plot(TR_list,flips(1,:)*180/pi,'b',TR_list,flips(2,:)*180/pi,'k')
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
    jmaif   = jmA0  * gampdf(TR_list + jmt0  , jmalpha , jmbeta);
    jmaifad = jmA0 * 1/(betamean^alphamean* gamma(alphamean)) *(TR_list+jmt0).^(alphamean-1).* exp(-(TR_list+jmt0)/ betamean);
    figure(4)
    plot(TR_list,jmaif ,'b', TR_list,jmaifad ,'r')
    ylabel('aif')
    xlabel('sec')
end


%% optimize MI for TR and FA
optf = true;
if optf
    % Convert this function file to an optimization expression.
    
    % Furthermore, you can also convert the |rosenbrock| function handle, which 
    % was defined at the beginning of the plotting routine, into an optimization expression.

    % setup optimization variables
    Nspecies = 2
    FaList = optimvar('FaList',Nspecies,Ntime,'LowerBound',0, 'UpperBound',35*pi/180);
    subFaList = optimexpr(      [Nspecies,(Ntime-1)*nsubstep+1 ]);
    subFaList(:,1:nsubstep:(Ntime-1)*nsubstep+1) = FaList ;
    %TRList = optimvar('TR',Ntime-1,1,'LowerBound',0, 'UpperBound',5);%TR_list;
    TRList = optimvar('TR','LowerBound',0, 'UpperBound',5);%TR_list;
    %TRList = TR;
    % [0;cumsum( TR* ones(Ntime-1,1))]

    NGauss  = 3
    NumberUncertain=3;
    switch (NumberUncertain)
       case(3)
         [x,xn,xm,w,wn]=GaussHermiteNDGauss(NGauss,[tisinput(5:2:9)],[tisinput(6:2:10)]);
       case(4)
         [x,xn,xm,w,wn]=GaussHermiteNDGauss(NGauss,[tisinput(1:2:7)],[tisinput(2:2:8)]);
    end 
    lqp=length(xn{1}(:));
    statevariable    = optimvar('state',Nspecies,(Ntime-1)*nsubstep+1 ,'LowerBound',0);
    auxvariable      = optimexpr(      [Nspecies,(Ntime-1)*nsubstep+1 ]);
    stateconstraint  = optimconstr(    [Nspecies,(Ntime-1)*nsubstep+1 ]);

    modelSNR = 20 ; % TODO - FIXME
    signuImage = (max(Mxy(1,:))+max(Mxy(2,:)))/2/modelSNR;
    % walker paper is peak pyruvate only
    signuImage = max(Mxy(1,:))/modelSNR;
    % variance for Gauss RV is sum. sqrt for std
    signu = sqrt(2* Ntime) * signuImage;
    [x2,xn2,xm2,w2,wn2]=GaussHermiteNDGauss(NGauss,0,signu);
    lqp2=length(xn2{1}(:));

    disp('build state variable')
    stateconstraint(:,1,:)  = statevariable(:,1,:) ==0;
    auxvariable(:,1,:) =0;
    switch (NumberUncertain)
       case(3)
         T1Pqp   = T1pmean;
         T1Lqp   = T1lmean;
         kplqp   = optimvar('kpl');
         klpqp   =    0 ;     % @cmwalker where do I get this from ? 
         kveqp   = optimvar('kve');
         t0qp    = optimvar('t0');
       %case(4)
       %  T1Pqp   = xn{1}(iqp);
       %  T1Lqp   = xn{2}(iqp);
       %  kplqp   = xn{3}(iqp);
       %  klpqp   =    0 ;     % @cmwalker where do I get this from ? 
       %  kveqp   = xn{4}(iqp);
       %  t0qp    = t0mean(1); 
    end 

    % precompute
    A = [(-kplqp -kveqp/ve -1/T1Pqp), 0; kplqp, -1/T1Lqp ]
    A_inv = inv(A)
    A_inv_sq = A_inv^2
    subTR = TRList /nsubstep;

    % loop over time
    for iii = 1:(Ntime-1)*nsubstep
      % >> syms a  kpl d subTR    T1P kveqp T1L 
      % >> expATR = expm([a,  0; kpl, d ] * subTR )
      % 
      % expATR =
      % 
      % [                                     exp(a*subTR),                0]
      % [(kpl*exp(a*subTR) - kpl*exp(subTR*d))/(a - d), exp(subTR*d)]
      % 
      % >> a = -1/T1P - kpl - kveqp
      % >> d = -1/T1L
      % >> eval(expATR)
      % 
      % ans =
      % 
      % [                                                              exp(-subTR*(kpl + kveqp + 1/T1P)),                   0]
      % [(kpl*exp(-subTR/T1L) - kpl*exp(-subTR*(kpl + kveqp + 1/T1P)))/(kpl + kveqp - 1/T1L + 1/T1P), exp(-subTR/T1L)]
      %    

      expAsubTR = [ exp(-subTR*(kplqp + kveqp/ve + 1/T1Pqp)),                   0; (kplqp*exp(-subTR/T1Lqp) - kplqp*exp(-subTR*(kplqp + kveqp/ve + 1/T1Pqp)))/(kplqp + kveqp/ve - 1/T1Lqp + 1/T1Pqp), exp(-subTR/T1Lqp)];
      aifterm = - kveqp/ve*A_inv*(eye(2) - expAsubTR)*VIF_scale_fact.*[gampdf(SubTimeList(iii)-t0qp , alphamean, betamean);0] ...
                + A_inv_sq*(expAsubTR-(A*subTR)-eye(2))*kveqp/ve/subTR*VIF_scale_fact.* [gampdf(SubTimeList(iii+1)-t0qp, alphamean, betamean)-gampdf(SubTimeList(iii)-t0qp, alphamean, betamean);0];
      auxvariable(:,iii+1) = expAsubTR*(cos(subFaList(:,iii)).*(auxvariable(:,iii))) + aifterm ;
      stateconstraint(:,iii+1) = statevariable(:,iii+1) ==  expAsubTR *(cos(subFaList(:,iii)).*statevariable(:,iii))   + aifterm ;
    end

    disp('build objective function')
    sumstatevariable =  sum(sum(sin(subFaList).*(ve*statevariable(:,:)  + (1-ve) *jmA0  * [gampdf( SubTimeList - t0qp  , jmalpha , jmbeta);zeros(1,(Ntime-1)*nsubstep+1)]  ),2));

    %% 
    % Create an optimization problem using these converted optimization expressions.
    
    disp('create optim prob')
    convprob = optimproblem('Objective',sumstatevariable , "Constraints",stateconstraint);
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
    %InitialGuess =  [flips(:);TR* ones(Ntime-1,1) ];   
    %pmin =  [flips(:)*0;zeros(Ntime-1,1)   ];     
    %pmax =  [flips(:)*0+35*pi/180;5*ones(Ntime-1,1) ];
    InitialGuess =  [flips(:);TR ];   
    pmin =  [flips(:)*0;0 ];     
    pmax =  [flips(:)*0+35*pi/180;5 ];
    %InitialGuess =  [flips(:) ];   
    %pmin =  [flips(:)*0 ];     
    %pmax =  [flips(:)*0+35*pi/180 ];
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
        %%maxiter,'Display','iter-detailed','Hessian',{'lbfgs',1}, 'HessMult',@myHessMultFcn,...
    tic;
    [designopt,fval,exitflag,output,lambda,grad,hessian] ...
     =fmincon(Fx, InitialGuess ,[],[],[],[],pmin,pmax,[],...
        optimset('TolX',tolx,'TolFun',tolfun,'MaxIter', ...
        maxiter,'Display','iter-detailed','Hessian',{'lbfgs',1}, ...
        'GradObj','on','PlotFcn',{'optimplotfvalconstr', 'optimplotconstrviolation', 'optimplotfirstorderopt' }) ...
        );

    toc;
    handle = figure(5)
      switch (NumberUncertain)
         case(3)
           optparams = struct('t0',[t0qp(1);0],'gammaPdfA',[alphamean(1)  ;1],'gammaPdfB',[betamean(1);1],...
               'scaleFactor',VIF_scale_fact,'T1s',[T1pmean(1),T1lmean(1)],'ExchangeTerms',[0,kplqp(1) ;0,0],...
               'TRList',TR_list,'PerfusionTerms',[kveqp(1),0],'volumeFractions',ve,...
               'fitOptions', opts);
         case(4)
           error("WIP")
         case(5)
           error("WIP")
      end
    optparams.FaList = reshape(designopt(myidx.FaList),size(params.FaList ));
    [t_axisopt,Mxyopt,Mzopt] = model.compile(M0.',params);
    mystate  = evaluate(auxvariable ,optparams);
    figure(6)
    plot(optparams.TRList,Mxyopt(1,:),'b',optparams.TRList,Mxyopt(2,:),'k')
    ylabel('MI Mxy')
    xlabel('sec')
    figure(7)
    plot(optparams.TRList,optparams.FaList(1,:)*180/pi,'b',optparams.TRList,optparams.FaList(2,:)*180/pi,'k')
    ylabel('MI FA')
    xlabel('sec')
    handle = figure(8)
    plot(optparams.TRList,Mzopt(1,:),'b',optparams.TRList,Mzopt(2,:),'k',optparams.TRList,cos(optparams.FaList(1,:)).*mystate(1,1:nsubstep:(Ntime-1)*nsubstep+1,1),'b--',optparams.TRList,cos(optparams.FaList(2,:)).*mystate(2,1:nsubstep:(Ntime-1)*nsubstep+1,1),'k--')
    ylabel('MI Mz ')
    xlabel('sec'); legend('Pyr','Lac')
end 

% SD
function W = myHessMultFcn(x,lambda,v)
   W = v;
end


function [MIobjfun, MIobjfun_Der]=MIGHQuadHPTofts(xopt,problem,myidx,Nspecies,Ntime,auxvariable)
    x0.FaList = reshape(xopt(myidx.FaList),Nspecies,Ntime);
    x0.TR     = xopt(myidx.TR);
    x0.state  = evaluate(auxvariable ,x0);
    %Xfull = [ x0.FaList(:); x0.TR(:); x0.state(:)];
    Xfull = [ x0.FaList(:); x0.TR; x0.state(:)];
    %Xfull = [ x0.FaList(:); x0.state(:)];
    [MIobjfun,initVals.g] = problem.objective(Xfull);
    [initConst.ineq,initConst.ceq, initConst.ineqGrad,initConst.ceqGrad] = problem.nonlcon(Xfull);
    objectiveGradFA    = initVals.g(myidx.FaList);
    objectiveGradTR    = initVals.g(myidx.TR);
    objectiveGradState = initVals.g(myidx.state);
    jacobianFA    = initConst.ceqGrad(myidx.FaList,:);
    jacobianTR    = initConst.ceqGrad(myidx.TR,:);
    jacobianState = initConst.ceqGrad(myidx.state,:);
    adjointvar =-jacobianState \objectiveGradState ;
    %MIobjfun_Der = [objectiveGradFA] +  [jacobianFA] *   adjointvar ;
    MIobjfun_Der = [objectiveGradFA;objectiveGradTR] +  [jacobianFA;jacobianTR    ] *   adjointvar ;
    for iqp = 1:lqp
    end
    for jjj = 1:lqp
    end 
    expandvar  = ones(1,lqp);
    diffsumm =(sumstatevariable(1,:)+sumstatevariable(2,:))' * expandvar   - expandvar' * (sumstatevariable(1,:)+sumstatevariable(2,:));
    Hz = 0;
    for jjj=1:lqp2
      znu=xn2{1}(jjj) ;
      Hz = Hz + wn2(jjj) * (wn(:)' * log(exp(-(znu + diffsumm).^2/2/signu^2 - log(signu) -log(2*pi)/2   ) * wn(:)));
    end
    %% MIGaussObj = Hz/sqrt(pi)^(NumberUncertain+1); 
    MIGaussObj = Hz;
end
