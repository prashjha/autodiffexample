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
    ai3 = jmalpha -1;
    zzz = (TR_list+jmt0)/jmbeta;
    jmaifad = zeros(size(TR_list));
    nonzeroidx  = find(zzz >0);
    jmaifad(nonzeroidx  ) = jmA0  * exp(-lnsr2pi -0.5*log(ai3) - stirlerr(ai3) - (ai3.*log(ai3./zzz(nonzeroidx  ))+zzz(nonzeroidx  )-ai3) ) ./ jmbeta;
    figure(4)
    plot(TR_list,jmaif ,'b', TR_list,jmaifad ,'r')
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
    syms FaList [Nspecies,Ntime] positive real
    %TRList = optimvar('TR',Ntime-1,1,'LowerBound',0, 'UpperBound',5);%TR_list;
    %TRList = optimvar('TR','LowerBound',0, 'UpperBound',5);%TR_list;
    TRList = TR;
    % [0;cumsum( TR* ones(Ntime-1,1))]

    NGauss  = 3
    NumberUncertain=3;
    switch (NumberUncertain)
       case(3)
         [x,xn,xm,w,wn]=GaussHermiteNDGauss(NGauss,[tisinput(5:2:9)],[tisinput(6:2:10)]);
         [x1,xn1,xm1,w1,wn1]=GaussHermiteNDGauss(NGauss,kplmean,kplstdd);
         [x2,xn2,xm2,w2,wn2]=GaussHermiteNDGauss(NGauss,kvemean,kvestdd);
         [x3,xn3,xm3,w3,wn3]=GaussHermiteNDGauss(NGauss,t0mean ,t0sttd );
       case(4)
         [x,xn,xm,w,wn]=GaussHermiteNDGauss(NGauss,[tisinput(1:2:7)],[tisinput(2:2:8)]);
    end 
    lqp=length(xn{1}(:));
    syms statevariable    [Nspecies,Ntime,NGauss,NGauss,NGauss] real
    syms auxstatevariable [Nspecies,Ntime,NGauss,NGauss,NGauss] real
    syms kplsymvar   [1 NGauss] positive real
    syms kvesymvar   [1 NGauss] positive real
    syms t0symvar    [1 NGauss] positive real

    modelSNR = 20 ; % TODO - FIXME
    signuImage = (max(Mxy(1,:))+max(Mxy(2,:)))/2/modelSNR;
    % walker paper is peak pyruvate only
    signuImage = max(Mxy(1,:))/modelSNR;
    % variance for Gauss RV is sum. sqrt for std
    signu = sqrt(2* Ntime) * signuImage;
    [x2,xn2,xm2,w2,wn2]=GaussHermiteNDGauss(NGauss,0,signu);
    lqp2=length(xn2{1}(:));

    disp('build state variable')
    statevariable(:,1,:,:,:) ==0;
    auxstatevariable(:,1,:,:,:) ==0;
    TimeList = (0:(Ntime-1))*TRList ;
    %TimeList = [0;cumsum( TRList)]
    for iqp = 1:NGauss
      for jqp = 1:NGauss
        for kqp = 1:NGauss
      iqp,jqp,kqp
      switch (NumberUncertain)
         case(3)
           T1Pqp   = T1pmean;
           T1Lqp   = T1lmean;
           kplqp   = kplsymvar(iqp);
           klpqp   =    0 ;     % @cmwalker where do I get this from ? 
           kveqp   = kvesymvar(jqp);
           t0qp    = t0symvar(kqp); 
         case(4)
           T1Pqp   = xn{1}(iqp);
           T1Lqp   = xn{2}(iqp);
           kplqp   = xn{3}(iqp);
           klpqp   =    0 ;     % @cmwalker where do I get this from ? 
           kveqp   = xn{4}(iqp);
           t0qp    = t0mean(1); 
      end 
      % setup AIF
      % ai3 = jmalpha -1;
      % zzz = (TimeList +t0qp)/jmbeta;
      % integrand = jmA0  * exp(-lnsr2pi -0.5*log(ai3) - stirlerr(ai3) - (ai3.*log(ai3./zzz)+zzz-ai3) ) ./ jmbeta;

      % loop over time
      for iii = 1:Ntime-1
        %currentTR = TRList(iii) ;
        currentTR = TRList ;
        nsubstep = 5;
        deltat = currentTR /nsubstep ;
        % setup AIF
        integratedt = TimeList(iii)+ [1:2:2*nsubstep]*deltat/2;
        ai3 = jmalpha -1;
        zzz = (integratedt +t0qp)/jmbeta;
        integrand = jmA0  * exp(-lnsr2pi -0.5*log(ai3) - stirlerr(ai3) - (ai3.*log(ai3./zzz)+zzz-ai3) ) ./ jmbeta;
        %integrand = jmA0 * gampdf(integratedt(1:nsubstep )'-t0qp,jmalpha,jmbeta) ;

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

        % 
        % mid-point rule integration gampdf
        %aifterm = kveqp/ve * deltat * [ exp((-1/T1Pqp - kplqp - kveqp/ve)*(TimeList(iii+1)-deltat*[.5:1:nsubstep] -TimeList(iii))); kplqp*(-exp((-1/T1Pqp - kplqp - kveqp/ve)*(TimeList(iii+1)-deltat*[.5:1:nsubstep] -TimeList(iii)) ) + exp(-1/T1Lqp *(TimeList(iii+1)-deltat*[.5:1:nsubstep] -TimeList(iii)) ))/(1/T1Pqp + kplqp + kveqp/ve - 1/T1Lqp )] * integrand ;
        % mid-point rule integration df
        aifterm = kveqp/ve * deltat * [ exp((-1/T1Pqp - kplqp - kveqp/ve)*deltat*[.5:1:nsubstep] ); kplqp*(-exp((-1/T1Pqp - kplqp - kveqp/ve)*deltat*[.5:1:nsubstep] ) + exp(-1/T1Lqp *deltat*[.5:1:nsubstep] ))/(1/T1Pqp + kplqp + kveqp/ve - 1/T1Lqp )] * integrand' ;
        % symbolic integration
        %aifterm =  [-(integrand(iii+1) - integrand(iii) - integrand(iii+1)*exp(TimeList(iii)*(kplqp + kveqp + 1/T1Pqp) - TimeList(iii+1)*(kplqp + kveqp + 1/T1Pqp)) + integrand(iii)*exp(TimeList(iii)*(kplqp + kveqp + 1/T1Pqp) - TimeList(iii+1)*(kplqp + kveqp + 1/T1Pqp)) - integrand(iii+1)*TimeList(iii+1)*(kplqp + kveqp + 1/T1Pqp) + integrand(iii+1)*TimeList(iii)*(kplqp + kveqp + 1/T1Pqp) + integrand(iii)*TimeList(iii+1)*exp(TimeList(iii)*(kplqp + kveqp + 1/T1Pqp) - TimeList(iii+1)*(kplqp + kveqp + 1/T1Pqp))*(kplqp + kveqp + 1/T1Pqp) - integrand(iii)*TimeList(iii)*exp(TimeList(iii)*(kplqp + kveqp + 1/T1Pqp) - TimeList(iii+1)*(kplqp + kveqp + 1/T1Pqp))*(kplqp + kveqp + 1/T1Pqp))/((TimeList(iii+1) - TimeList(iii))*(kplqp + kveqp + 1/T1Pqp)^2) ; -(integrand(iii+1) - integrand(iii) - integrand(iii+1)*exp(kplqp*TimeList(iii+1) - kplqp*TimeList(iii)) + integrand(iii)*exp(kplqp*TimeList(iii+1) - kplqp*TimeList(iii)) + integrand(iii+1)*kplqp*TimeList(iii+1) - integrand(iii+1)*kplqp*TimeList(iii) - integrand(iii)*kplqp*TimeList(iii+1)*exp(kplqp*TimeList(iii+1) - kplqp*TimeList(iii)) + integrand(iii)*kplqp*TimeList(iii)*exp(kplqp*TimeList(iii+1) - kplqp*TimeList(iii)))/(kplqp^2*(TimeList(iii+1) - TimeList(iii)))];

        expATR = [ exp(-currentTR*(kplqp + kveqp/ve + 1/T1Pqp)),                   0; (kplqp*exp(-currentTR/T1Lqp) - kplqp*exp(-currentTR*(kplqp + kveqp/ve + 1/T1Pqp)))/(kplqp + kveqp/ve - 1/T1Lqp + 1/T1Pqp), exp(-currentTR/T1Lqp)];
        statevariable(:,iii+1,iqp,jqp,kqp) =  expATR *(cos(FaList(:,iii+1)).* statevariable(:,iii,iqp,jqp,kqp)) +  aifterm ;
        auxstatevariable(:,iii+1,iqp,jqp,kqp) =  expATR *(cos(FaList(:,iii+1)).* auxstatevariable(:,iii,iqp,jqp,kqp)) +  aifterm ;
      end
    end
    end
    end
    %% % debug
    %% dbgparams = struct('t0',[xn{3}(1);0],'gammaPdfA',[alphamean(1)  ;1],'gammaPdfB',[betamean(1);1],...
    %%     'scaleFactor',VIF_scale_fact,'T1s',[T1pmean(1),T1lmean(1)],'ExchangeTerms',[0,xn{1}(1) ;0,0],...
    %%     'TRList',TR_list,'PerfusionTerms',[xn{2}(1),0],'volumeFractions',ve, 'fitOptions', opts);
    %% dbgparams.FaList = flips;
    %% dbgparams.TR = TR* ones(Ntime-1,1);
    %% [t_axis,Mxydbg,Mzdbg] = model.compile(M0.',dbgparams);
    %% dbgparams.FaList = flips;
    %% mystate  = evaluate(auxvariable ,dbgparams);
    %% figure(11)
    %% plot(dbgparams.TRList,Mzdbg(1,:),'b',dbgparams.TRList,Mzdbg(2,:),'k',dbgparams.TRList,mystate(1,:,1),'b--',dbgparams.TRList,mystate(2,:,1),'k--')

    disp('build objective function')
    syms sumstatevariable [NGauss,NGauss,NGauss] real;
    syms sumauxstatevariable [NGauss,NGauss,NGauss] real;
    for iqp = 1:NGauss
      for jqp = 1:NGauss
        for kqp = 1:NGauss
       sumstatevariable(iqp,jqp,kqp) =  sum(sum(sin(FaList) .*    statevariable(:,:,iqp,jqp,kqp),2),1);
       sumauxstatevariable(iqp,jqp,kqp) =  sum(sum(sin(FaList) .*   auxstatevariable(:,:,iqp,jqp,kqp),2),1);
       %sumstatevariable(:,jjj) =  sum(sin(FaList).*(ve*statevariable(:,:,jjj)  + (1-ve) *jmA0  * [gampdf( TimeList - t0qp  , jmalpha , jmbeta);zeros(1,Ntime)]  ),2);
    end 
    end 
    end 
    %statematrix = optimexpr([lqp,lqp]);
    expandvar  = ones(1,lqp);
    %diffsumm =(sumstatevariable(1,:)+sumstatevariable(2,:))' * expandvar   - expandvar' * (sumstatevariable(1,:)+sumstatevariable(2,:));
    Hz = 0;
    %for jjj=1:lqp2
    %  znu=xn2{1}(jjj) ;
    %  Hz = Hz + wn2(jjj) * (wn(:)' * log(exp(-(znu + diffsumm).^2/2/signu^2 - log(signu) -log(2*pi)/2   ) * wn(:)));
    %end
    %for iii=1:lqp
    myintegrand = exp(-(sumstatevariable- sumauxstatevariable).^2/sqrt(2)/signu);
    integrandsum  = sum(repmat(wn3,1,3,3).*subs(myintegrand , t0symvar',xn3{1}),3);
    integrandsum2 = sum(repmat(wn2,1,3).*  subs(integrandsum,kvesymvar',xn2{1}),2);
    integrandsum3 = sum(wn1.*              subs(integrandsum,kplsymvar',xn1{1})  );
    %% for iqp = 1:NGauss
    %%   for jqp = 1:NGauss
    %%     for kqp = 1:NGauss
    %%     for jjj=1:lqp2
    %%         znu=xn2{1}(jjj) ;
    %%         lntermtmp=0;
    %%         %for kkk=1:lqp
    %%         for lqp = 1:NGauss
    %%           for mqp = 1:NGauss
    %%           lntermtmp  = wn3' *subs(exp(-(sumstatevariable(1,1,1,1)- sumstatevariable(1,1,1,:)).^2/sqrt(2)/signu)  , t0symvar', xn3{1} )
    %%          % for nqp = 1:NGauss
    %%          %   lntermtmp=lntermtmp + wn1(lqp) *wn2(mqp) *wn3(nqp) * exp(-(znu+sumstatevariable(1,iqp,jqp,kqp)- sumstatevariable(1,lqp,mqp,nqp))^2/sqrt(2)/signu);
    %%          %   lntermtmp=lntermtmp + wn1(lqp) *wn2(mqp) *wn3(nqp) * exp(-(znu+sumstatevariable(2,iqp,jqp,kqp)- sumstatevariable(2,lqp,mqp,nqp))^2/sqrt(2)/signu);
    %%          % end
    %%         end
    %%         end
    %%         % for function
    %%         lnterm = log(lntermtmp)+log(pi^(-1.5));
    %%         Hz = Hz + wn1(iqp) *wn2(jqp) *wn3(kqp) * wn2(jjj) * lnterm;
    %%     end
    %% end
    %% end
    %% end
    %% MIGaussObj = Hz/sqrt(pi)^(NumberUncertain+1); 
    MIGaussObj = integrandsum3 ;

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
    InitialGuess =  [flips(:);TR* ones(Ntime-1,1) ];   
    pmin =  [flips(:)*0;zeros(Ntime-1,1)   ];     
    pmax =  [flips(:)*0+35*pi/180;5*ones(Ntime-1,1) ];
    %InitialGuess =  [flips(:);TR ];   
    %pmin =  [flips(:)*0;0 ];     
    %pmax =  [flips(:)*0+35*pi/180;5 ];
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
    figure(6)
    plot(optparams.TRList,Mxyopt(1,:),'b',optparams.TRList,Mxyopt(2,:),'k')
    ylabel('MI Mxy')
    xlabel('sec')
    figure(7)
    plot(optparams.TRList,optparams.FaList(1,:)*180/pi,'b',optparams.TRList,optparams.FaList(2,:)*180/pi,'k')
    ylabel('MI FA')
    xlabel('sec')
    handle = figure(8)
    plot(optparams.TRList,Mzopt(1,:),'b',optparams.TRList,Mzopt(2,:),'k')
    plot(optparams.TRList,Mzopt(1,:),'b',optparams.TRList,Mzopt(2,:),'k')
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
    Xfull = [ x0.FaList(:); x0.TR(:); x0.state(:)];
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
end
