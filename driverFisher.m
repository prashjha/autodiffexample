% Maximize signal

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
NumberUncertain = 7;

%% Variable Setup
lnsr2pi = 0.9189385332046727; % log(sqrt(2*pi))
Ntime = 30;
TR = 3;
TR_list = (0:(Ntime-1))*TR;
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
    % ai3 = jmalpha -1;
    % zzz = (TR_list+jmt0)/jmbeta;
    % jmaifad = zeros(size(TR_list));
    % nonzeroidx  = find(zzz >0);
    % jmaifad(nonzeroidx  ) = jmA0  * exp(-lnsr2pi -0.5*log(ai3) - stirlerr(ai3) - (ai3.*log(ai3./zzz(nonzeroidx  ))+zzz(nonzeroidx  )-ai3) ) ./ jmbeta;
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
    
    % Furthermore, you can also convert the |rosenbrock| function handle, which 
    % was defined at the beginning of the plotting routine, into an optimization expression.

    % setup optimization variables
    Nspecies = 2
    %FaList = optimvar('FaList',Nspecies,Ntime,'LowerBound',0, 'UpperBound',35*pi/180);
    FaList = flips;
    %TRList = optimvar('TR',Ntime-1,1,'LowerBound',0, 'UpperBound',5);%TR_list;
    %TRList = optimvar('TR','LowerBound',0, 'UpperBound',5);%TR_list;
    TRList = TR;
    % [0;cumsum( TR* ones(Ntime-1,1))]

    statevariable    = optimvar('state',Nspecies,Ntime,'LowerBound',0);
    auxvariable      = optimexpr(    [Nspecies,Ntime]);
    stateconstraint  = optimconstr(  [Nspecies,Ntime]);

    modelSNR = 20 ; % TODO - FIXME
    signuImage = (max(Mxy(1,:))+max(Mxy(2,:)))/2/modelSNR;
    % walker paper is peak pyruvate only
    signuImage = max(Mxy(1,:))/modelSNR;
    % variance for Gauss RV is sum. sqrt for std
    signu = sqrt(2* Ntime) * signuImage;

    disp('build state variable')
    stateconstraint(:,1)  = statevariable(:,1)==0;
    auxvariable(:,1) =0;
    TimeList = (0:(Ntime-1))*TRList ;
    %TimeList = [0;cumsum( TRList)]

      klpqp   =    0 ;  
    % sample point
      samplepoint = randn(NumberUncertain ,1).* tisinput(2:2:14) + tisinput(1:2:13);
      T1Pqp   = optimvar('T1P');
      T1Lqp   = optimvar('T1L');
      kplqp   = optimvar('kpl');
      kveqp   = optimvar('kve');
      t0qp    = t0mean;
      alphaqp = alphamean;
      betaqp  = betamean ;

      % loop over time
      for iii = 1:Ntime-1
        %currentTR = TRList(iii) ;
        currentTR = TRList ;
        nsubstep = 5;
        deltat = currentTR /nsubstep ;
        % setup AIF
        integratedt = TimeList(iii)+ [1:2:2*nsubstep]*deltat/2;
        integrand = jmA0 * gampdf(integratedt(1:nsubstep )'-t0qp,alphaqp,betaqp) ;

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
        % >> syms tk tkm1 tau a gammak gammakm1 kplqp d  T1Pqp kveqp T1Lqp 
        % >> aifexpr = exp([a,  0; kplqp, d ] *(tk-tau)) * [gammak * (tau - tkm1)/(tk-tkm1 ) + gammakm1 * (tk - tau )/(tk-tkm1 ) ;0]
        % 
        % aifexpr =
        % 
        %     -exp(-a*(tau - tk))*((gammakm1*(tau - tk))/(tk - tkm1) - (gammak*(tau - tkm1))/(tk - tkm1))
        % -exp(-kplqp*(tau - tk))*((gammakm1*(tau - tk))/(tk - tkm1) - (gammak*(tau - tkm1))/(tk - tkm1))
        % 
        % >> aifint = int(aifexpr,tau,tkm1,tk)
        % >> a = -1/T1Pqp - kplqp - kveqp
        % >> d = -1/T1Lqp
        % >> simplify(eval(aifint ))
        % 
        % ans =
        % 
        % -(gammak - gammakm1 - gammak*exp(tkm1*(kplqp + kveqp + 1/T1Pqp) - tk*(kplqp + kveqp + 1/T1Pqp)) + gammakm1*exp(tkm1*(kplqp + kveqp + 1/T1Pqp) - tk*(kplqp + kveqp + 1/T1Pqp)) - gammak*tk*(kplqp + kveqp + 1/T1Pqp) + gammak*tkm1*(kplqp + kveqp + 1/T1Pqp) + gammakm1*tk*exp(tkm1*(kplqp + kveqp + 1/T1Pqp) - tk*(kplqp + kveqp + 1/T1Pqp))*(kplqp + kveqp + 1/T1Pqp) - gammakm1*tkm1*exp(tkm1*(kplqp + kveqp + 1/T1Pqp) - tk*(kplqp + kveqp + 1/T1Pqp))*(kplqp + kveqp + 1/T1Pqp))/((tk - tkm1)*(kplqp + kveqp + 1/T1Pqp)^2)
        % -(gammak - gammakm1 - gammak*exp(kplqp*(tk - tkm1)) + gammakm1*exp(kplqp*(tk - tkm1)) + gammak*kplqp*tk - gammak*kplqp*tkm1 - gammakm1*kplqp*tk*exp(kplqp*(tk - tkm1)) + gammakm1*kplqp*tkm1*exp(kplqp*(tk - tkm1)))/(kplqp^2*(tk - tkm1))

        % mid-point rule integration gampdf
        aifterm = kveqp/ve * deltat * [ exp((-1/T1Pqp - kplqp - kveqp/ve)*(TimeList(iii+1)-deltat*[.5:1:nsubstep] -TimeList(iii))); kplqp*(-exp((-1/T1Pqp - kplqp - kveqp/ve)*(TimeList(iii+1)-deltat*[.5:1:nsubstep] -TimeList(iii)) ) + exp(-1/T1Lqp *(TimeList(iii+1)-deltat*[.5:1:nsubstep] -TimeList(iii)) ))/(1/T1Pqp + kplqp + kveqp/ve - 1/T1Lqp )] * integrand ;

        % evaluate signal model
        expATR = [ exp(-currentTR*(kplqp + kveqp/ve + 1/T1Pqp)),                   0; (kplqp*exp(-currentTR/T1Lqp) - kplqp*exp(-currentTR*(kplqp + kveqp/ve + 1/T1Pqp)))/(kplqp + kveqp/ve - 1/T1Lqp + 1/T1Pqp), exp(-currentTR/T1Lqp)];
        auxvariable(:,iii+1) =  expATR *(cos(FaList(:,iii)).*auxvariable(:,iii ))   + aifterm ;
        stateconstraint(:,iii+1) = statevariable(:,iii+1) ==  expATR *(cos(FaList(:,iii)).*statevariable(:,iii ))   + aifterm ;
      end

    disp('build objective function')
    sumstatevariable = optimexpr([Nspecies]);
    sumstatevariable =  sum(sin(FaList).*(ve*statevariable  + (1-ve) *jmA0  * [gampdf( TimeList - t0qp  , alphaqp , betaqp);zeros(1,Ntime)]  ),2);
    %statematrix = optimexpr([lqp,lqp]);
    signalmodel =  sum(sumstatevariable,1);

    %% 
    % Create an optimization problem using these converted optimization expressions.
    
    disp('create optim prob')
    convprob = optimproblem('Objective',signalmodel , "Constraints",stateconstraint);
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
    %InitialGuess =  [flips(:);TR ];   
    %pmin =  [flips(:)*0;0 ];     
    %pmax =  [flips(:)*0+35*pi/180;5 ];
    xinit.T1P =T1pmean;
    xinit.T1L =T1lmean;
    xinit.kpl =kplmean;
    xinit.kve =kvemean;
    %pmin =  [flips(:)*0 ];     
    %pmax =  [flips(:)*0+35*pi/180 ];
    tolx=1.e-9;
    tolfun=5.e-4;
    maxiter=400;

    Fx = @(x) hpConditionalProbability(x, problem, myidx,Nspecies,Ntime,auxvariable);
    %% debug info
    %% x0.FaList = params.FaList;
    %% x0.state  = evaluate(auxvariable ,x0);
    %% mystate = evaluate( sumstatevariable ,x0);
    %% Xfull = [ x0.FaList(:); x0.state(:)];
    %% [MIobjfun,initVals.g] = problem.objective(Xfull);
    %% [initConst.ineq,initConst.ceq, initConst.ineqGrad,initConst.ceqGrad] = problem.nonlcon(Xfull);
    [initobj,initgrad] = Fx(xinit);
    fishermatrix = 1/signu * (initgrad * initgrad')
    crlb = inv(fishermatrix ) 
    %% [designopt,fval,exitflag,output,lambda,grad,hessian] ...
    %%  =fmincon(Fx, InitialGuess ,[],[],[],[],pmin,pmax,[],...
    %%     optimset('TolX',tolx,'TolFun',tolfun,'MaxIter', ...
    %%     maxiter,'Display','iter-detailed','Hessian',{'lbfgs',1}, ...
    %%     'GradObj','on','PlotFcn',{'optimplotfvalconstr', 'optimplotconstrviolation', 'optimplotfirstorderopt' }) ...
    %%     );

    %% toc;
    %% handle = figure(5)
    %% optparams = struct('t0',[t0qp;0],'gammaPdfA',[alphaqp  ;1],'gammaPdfB',[betaqp;1],...
    %%     'scaleFactor',VIF_scale_fact,'T1s',[T1Pqp,T1Lqp ],'ExchangeTerms',[0,kplqp ;0,0],...
    %%     'TRList',TR_list,'PerfusionTerms',[kveqp,0],'volumeFractions',ve,...
    %%     'fitOptions', opts);
    %% optparams.FaList = reshape(designopt(myidx.FaList),size(params.FaList ));
    %% [t_axisopt,Mxyopt,Mzopt] = model.compile(M0.',optparams);
    %% mystate  = evaluate(auxvariable ,optparams);
    %% figure(6)
    %% plot(optparams.TRList,Mxyopt(1,:),'b',optparams.TRList,Mxyopt(2,:),'k')
    %% ylabel('MI Mxy')
    %% xlabel('sec')
    %% figure(7)
    %% plot(optparams.TRList,optparams.FaList(1,:)*180/pi,'b',optparams.TRList,optparams.FaList(2,:)*180/pi,'k')
    %% ylabel('MI FA')
    %% xlabel('sec')
    %% handle = figure(8)
    %% plot(optparams.TRList,Mzopt(1,:),'b',optparams.TRList,Mzopt(2,:),'k',optparams.TRList,cos(optparams.FaList(1,:)).*mystate(1,:),'b--',optparams.TRList,cos(optparams.FaList(2,:)).*mystate(2,:),'k--')
    %% ylabel('MI Mz ')
    %% xlabel('sec'); legend('Pyr','Lac')
end 

function [MIobjfun, MIobjfun_Der]=hpConditionalProbability(xopt,problem,myidx,Nspecies,Ntime,auxvariable)
    %x0.TR     = xopt(myidx.TR);
    x0.state  = evaluate(auxvariable ,xopt);
    %Xfull = [ x0.FaList(:); x0.TR(:); x0.state(:)];
    Xfull = [xopt.T1L; xopt.T1P; xopt.kpl; xopt.kve; x0.state(:)];
      
    [MIobjfun,initVals.g] = problem.objective(Xfull);
    [initConst.ineq,initConst.ceq, initConst.ineqGrad,initConst.ceqGrad] = problem.nonlcon(Xfull);
    objectiveGradFA    = initVals.g([myidx.T1L; myidx.T1P; myidx.kpl; myidx.kve]);
    
    %objectiveGradTR    = initVals.g(myidx.TR);
    objectiveGradState = initVals.g(myidx.state);
    jacobianFA    = initConst.ceqGrad([myidx.T1L; myidx.T1P; myidx.kpl; myidx.kve],:);
    %jacobianTR    = initConst.ceqGrad(myidx.TR,:);
    jacobianState = initConst.ceqGrad(myidx.state,:);
    adjointvar =-jacobianState \objectiveGradState ;
    MIobjfun_Der = [objectiveGradFA] +  [jacobianFA] *   adjointvar ;
    %MIobjfun_Der = [objectiveGradFA;objectiveGradTR] +  [jacobianFA;jacobianTR    ] *   adjointvar ;
end
