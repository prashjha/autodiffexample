% ShowMxyPub Shows the flip angle schedules and the resulting magnetization
% evolutions for the paper
% Author: Chris Walker
% Date: 8/6/2018

clear all
clc

myoptions = optimoptions(@fmincon,'Display','iter-detailed','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'MaxFunctionEvaluations',1e7,'ConstraintTolerance',2.e-9, 'OptimalityTolerance',2.5e-9,'Algorithm','interior-point','StepTolerance',1.000000e-12,'MaxIterations',1000,'PlotFcn',{'optimplotfvalconstr', 'optimplotconstrviolation', 'optimplotfirstorderopt' },'HonorBounds',true, 'HessianApproximation', 'lbfgs' ,'Diagnostic','on','FunValCheck','on' )
driverHPMIopt(3,3, 2,myoptions)
driverHPMIopt(3,3, 5,myoptions)
driverHPMIopt(3,3, 8,myoptions)
driverHPMIopt(3,3,10,myoptions)
driverHPMIopt(3,3,12,myoptions)
driverHPMIopt(3,3,15,myoptions)
driverHPMIopt(3,3,20,myoptions)
driverHPMIopt(3,3,22,myoptions)
driverHPMIopt(3,3,25,myoptions)
driverHPMIopt(4,3, 2,myoptions)
driverHPMIopt(4,3, 5,myoptions)
driverHPMIopt(4,3, 8,myoptions)
driverHPMIopt(4,3,10,myoptions)
driverHPMIopt(4,3,12,myoptions)
driverHPMIopt(4,3,15,myoptions)
driverHPMIopt(4,3,20,myoptions)
driverHPMIopt(4,3,22,myoptions)
driverHPMIopt(4,3,25,myoptions)
%driverHPMIopt(5,3,10,myoptions)
      %myoptions = optimoptions(@fmincon,'Display','iter-detailed','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'MaxFunctionEvaluations',1e7,'ConstraintTolerance',2.e-6, 'OptimalityTolerance',2.5e-6,'Algorithm','interior-point','StepTolerance',1.000000e-12,'MaxIterations',1000,'PlotFcn',{'optimplotfvalconstr', 'optimplotconstrviolation', 'optimplotfirstorderopt' },'SubproblemAlgorithm','cg','HonorBounds',false, 'HessianApproximation', 'finite-difference' ,'Diagnostic','on','FunValCheck','on','BarrierParamUpdate','predictor-corrector' )
      %myoptions = optimoptions(@fmincon,'Display','iter-detailed','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'MaxFunctionEvaluations',1e7,'ConstraintTolerance',1.e-7, 'OptimalityTolerance',1.e-16,'Algorithm','active-set','StepTolerance',1.000000e-16)
myoptions = optimoptions(@fmincon,'Display','iter-detailed','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'MaxFunctionEvaluations',1e7,'ConstraintTolerance',1.e-14, 'OptimalityTolerance',1.e-14,'Algorithm','sqp','StepTolerance',1.000000e-12,'MaxIterations',1000,'PlotFcn',{'optimplotfvalconstr', 'optimplotconstrviolation', 'optimplotfirstorderopt' },'SubproblemAlgorithm','cg')
driverHPMIopt(3,3, 2,myoptions)
driverHPMIopt(3,3, 5,myoptions)
driverHPMIopt(3,3, 8,myoptions)
driverHPMIopt(3,3,10,myoptions)
driverHPMIopt(3,3,12,myoptions)
driverHPMIopt(3,3,15,myoptions)
driverHPMIopt(3,3,20,myoptions)
driverHPMIopt(3,3,22,myoptions)
driverHPMIopt(3,3,25,myoptions)
driverHPMIopt(4,3, 2,myoptions)
driverHPMIopt(4,3, 5,myoptions)
driverHPMIopt(4,3, 8,myoptions)
driverHPMIopt(4,3,10,myoptions)
driverHPMIopt(4,3,12,myoptions)
driverHPMIopt(4,3,15,myoptions)
driverHPMIopt(4,3,20,myoptions)
driverHPMIopt(4,3,22,myoptions)
driverHPMIopt(4,3,25,myoptions)
%driverHPMIopt(5,3,10,myoptions)
% monitor memory: while [ -e /proc/3291925 ] ; do  top -b -n 1 -p 3291925 >>process.txt ;sleep 60; done  

function driverHPMIopt(NGauss,NumberUncertain,modelSNR,myoptions)

  NGauss,NumberUncertain,modelSNR,myoptions.Algorithm
  close all
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
  FAType = {'WarmStart'};
  FAType = {'Const'};
  %% HACK- @cmwalker code for initial conditions - https://github.com/fuentesdt/TumorHPMRI/blob/master/models/gPC/walker/ShowMxyPub.m
  for i = 1:numel(FAType)
      switch (FAType{i})
          case('Const') % Nagashima for lactate const 10 pyruvate
              tic
              E1(1) = exp(-TR*(1/T1pmean+kplmean));
              E1(2) = exp(-TR/T1lmean);
              for n = 1:Ntime
                  %flips(2,n) = acos(sqrt((E1(2)^2-E1(2)^(2*(N-n+1)))/(1-E1(2)^(2*(N-n+1)))));
                  flips(2,n) = 15*pi/180;
                  flips(1,n) = 20*pi/180;
              end
              dbgoptflips  =[ 0.3873    0.3876    0.4299    0.5295    0.5512    0.5510    0.5618    0.5783    0.5772    0.5890    0.5853 0.5847    0.5657    0.5724    0.6042    0.5874    0.5597    0.5299    0.5020    0.4777    0.4578    0.4420    0.4297; 0.2236    0.2236    0.2297    0.2691    0.3527    0.4070    0.4150    0.4132    0.4151    0.4201    0.4267 0.4346    0.4433    0.4534    0.4657    0.4825    0.5071    0.5418    0.5779    0.5999    0.6040    0.6058    0.6068 ];
  
              params.FaList = dbgoptflips ;
              params.FaList = flips ;
          case('WarmStart') % Nagashima for lactate const 10 pyruvate
              warmstart = load('poptNG4Nu3sqpSNR10.mat');
              params.FaList = warmstart.popt.FaList;
      end
  
      
      
  
      tic
      %% Fitting
      [t_axis,Mxy,Mz] = model.compile(M0.',params);
      toc
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
      plot(TR_list,params.FaList(1,:)*180/pi,'b',TR_list,params.FaList(2,:)*180/pi,'k')
      ylabel('Const FA (deg) ')
      xlabel('sec')
      figure(3)
      plot(TR_list,Mz(1,:),'b--',TR_list,Mz(2,:),'k--')
      hold
      plot(TR_list,Mz(1,:)./cos(params.FaList(1,:)),'b',TR_list,Mz(2,:)./cos(params.FaList(2,:)),'k')
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
      % setup optimization variables
      Nspecies = 2
      FaList = optimvar('FaList',Nspecies,Ntime,'LowerBound',0, 'UpperBound',35*pi/180);
      TRList = TR_list;
      diffTR = diff(TRList);
  
      % noise calc for signal sum - assume same noise in pyruvate and lactate image
      %    modelSNR = (maxsignallac + maxsignalpyr)/(stdsignalpyr +stdsignallac  )
      %             = (maxsignallac + maxsignalpyr)/(2 * signuImage )
      %
      % image level ==> signuImage = (maxsignallac + maxsignalpyr)/2/modelSNR 
      % signal sum  ==> signu = Ntime * signuImage 
      % noise calc max signal assuming total signal is sum of gaussian RV
      signuImage = (max(Mxy(1,:))+max(Mxy(2,:)))/2/modelSNR;
      signu = Ntime * signuImage;
      [x2,xn2,xm2,w2,wn2]=GaussHermiteNDGauss(NGauss,0,signu);
      lqp2=length(xn2{1}(:));
  
      switch (NumberUncertain)
         case(3)
           [x,xn,xm,w,wn]=GaussHermiteNDGauss(NGauss,[tisinput(5:2:9)],[tisinput(6:2:10)]);
           
           T1Lqp   = T1lmean;
           kplqp   = xn{1}(:);
           klpqp   =    0 ;     % @cmwalker where do I get this from ? 
           kveqp   = xn{2}(:);
           t0qp    = xn{3}(:); 
         case(4)
           [x,xn,xm,w,wn]=GaussHermiteNDGauss(NGauss,[tisinput(1:2:7)],[tisinput(2:2:8)]);
           T1Pqp   = xn{1}(:);
           T1Lqp   = xn{2}(:);
           kplqp   = xn{3}(:);
           klpqp   =    0 ;     % @cmwalker where do I get this from ? 
           kveqp   = xn{4}(:);
           t0qp    = t0mean(1); 
         case(5)
           [x,xn,xm,w,wn]=GaussHermiteNDGauss(NGauss,[tisinput(1:2:9)],[tisinput(2:2:10)]);
           T1Pqp   = xn{1}(:);
           T1Lqp   = xn{2}(:);
           kplqp   = xn{3}(:);
           klpqp   =    0 ;     % @cmwalker where do I get this from ? 
           kveqp   = xn{4}(:);
           t0qp    = xn{5}(:); 
      end
      %alphaqp = xn{6}(:); 
      %betaqp  = xn{7}(:); 
  
      lqp=length(xn{1}(:));
      %statevariable    = optimvar('state',Ntime,Nspecies,lqp,'LowerBound',0,'UpperBound',.1);
      statevariableraw  = optimvar('state',Ntime,Nspecies,lqp,'LowerBound',0);
      stateconstraint  = optimconstr(    [Ntime,Nspecies,lqp]);
  
      % scaling important for the optimizaiton step length update
      scalestate = 1.e-2;
      statevariable =scalestate * statevariableraw;
  
      disp('build state variable')
      
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
      expATRoneone = exp(-currentTR*(kplqp + kveqp + T1Pqp.^(-1)));
      expATRtwoone = (kplqp.*exp(-currentTR*T1Lqp.^(-1)) - kplqp.*exp(-currentTR*(kplqp + kveqp + T1Pqp.^(-1)))).* (kplqp + kveqp - T1Lqp.^(-1) + T1Pqp.^(-1)).^(-1);
      expATRtwotwo = exp(-currentTR * T1Lqp.^(-1));
       
      % IC
      stateconstraint(1,:,:)  = statevariable(1,:,:) ==0;
      for iii = 1:Ntime-1
          currentTR = diffTR(iii);
          nsubstep = 5;
          deltat = currentTR /nsubstep ;
          % setup AIF
          integratedt =TRList(iii)+ [1:2:2*nsubstep+1]*deltat/2;
          %integrand = jmA0 * my_gampdf(integratedt(1:nsubstep )'-t0qp,jmalpha,jmbeta) ;
          integrand = jmA0 * gampdf(repmat(integratedt(1:nsubstep )',1,lqp)'- repmat(t0qp,1,nsubstep),jmalpha,jmbeta) ;
          aiftermpyr = deltat * kveqp.*  [ exp(- T1Pqp.^(-1) - kplqp - kveqp)*deltat*[.5:1:nsubstep]  ].* integrand ; 
          aiftermlac = deltat * kveqp.*  ([ (-kplqp.*exp((-T1Pqp.^(-1) - kplqp - kveqp) ) + kplqp.*exp(-T1Lqp.^(-1) )).* ((T1Pqp.^(-1) + kplqp + kveqp) - T1Lqp.^(-1) ).^(-1)] *deltat*[.5:1:nsubstep]  ).* integrand ; 
  
          % setup state as linear constraint
          stateconstraint(iii+1,1,:)  = statevariable(iii+1,1,:) ==  reshape(cos(FaList(1,iii))*expATRoneone.* squeeze( statevariable(iii,1,: ) ),1,1,lqp ) +  reshape( sum(aiftermpyr,2 ),1,1,lqp) ;
          stateconstraint(iii+1,2,:)  = statevariable(iii+1,2,:) ==  reshape(cos(FaList(2,iii))*expATRtwotwo.* squeeze( statevariable(iii,2,: ) ),1,1,lqp ) + reshape( sum(aiftermlac,2 ),1,1,lqp) +reshape( cos(FaList(1,iii))*expATRtwoone.* squeeze( statevariable(iii,1,: )  ),1,1,lqp) ; 
      end
  
      disp('build objective function')
      % TODO - repmat does not work well with AD
      % TODO - replace repmat with matrix
      % sumstatevariable = squeeze(sum(repmat(sin(FaList)',1,1,lqp).*statevariable,1));
      sumstatevariable = optimexpr([Nspecies,lqp]);
      for jjj = 1:lqp
         sumstatevariable(:,jjj) =  sum(sin(FaList)'.*statevariable(:,:,jjj),1)';
      end 
      %statematrix = optimexpr([lqp,lqp]);
      %lqpchoosetwo = nchoosek(1:lqp,2);
      %arraypermutationsjjj = repmat([1:lqp]',1,lqp) ;
      %arraypermutationsiii = repmat([1:lqp] ,lqp,1) ;
      %lqpchoosetwo = [arraypermutationsiii(:), arraypermutationsjjj(:)];
      %diffsummone = sumstatevariable(1,lqpchoosetwo(:,1)) - sumstatevariable(1,lqpchoosetwo(:,2));
      %diffsummtwo = sumstatevariable(2,lqpchoosetwo(:,1)) - sumstatevariable(2,lqpchoosetwo(:,2));
      %diffsummone = repmat(sumstatevariable(1,:)',1,lqp) - repmat(sumstatevariable(1,:) ,lqp,1);
      %diffsummtwo = repmat(sumstatevariable(2,:)',1,lqp) - repmat(sumstatevariable(2,:) ,lqp,1);
      expandvar  = ones(1,lqp);
      diffsummone = sumstatevariable(1,:)' * expandvar   - expandvar' * sumstatevariable(1,:);
      diffsummtwo = sumstatevariable(2,:)' * expandvar   - expandvar' * sumstatevariable(2,:);
  
      Hz = 0;
      for jjj=1:lqp2
        znu=xn2{1}(jjj) ;
        %Hz = Hz + wn2(jjj) * (wn(lqpchoosetwo(:,1))' * log(exp(-(znu + diffsummone').^2/sqrt(2)/signu   - (znu + diffsummtwo').^2/sqrt(2)/signu  ).* wn(lqpchoosetwo(:,2))));
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
      problem = prob2struct(convprob,'ObjectiveFunctionName','generatedObjective');
      %% 
      % Solve the new problem. The solution is essentially the same as before.
      
      x0.FaList = params.FaList;
      x0.state  = repmat(1/scalestate * ( Mz./cos(params.FaList))',1,1,lqp);
      %'HessianApproximation', 'lbfgs'
  
      % truthconstraint = infeasibility(stateconstraint,x0);
      [popt,fval,exitflag,output] = solve(convprob,x0,'Options',myoptions, 'ConstraintDerivative', 'auto-reverse', 'ObjectiveDerivative', 'auto-reverse' )
  
      toc;
      % save convergence history
      handle = figure(5)
      saveas(handle,sprintf('historyNG%dNu%d%sSNR%02d',NGauss,NumberUncertain,myoptions.Algorithm,modelSNR ),'png')
      % save solution
      optparams = params;
      optparams.FaList = popt.FaList;
      [t_axisopt,Mxyopt,Mzopt] = model.compile(M0.',optparams);
      save(sprintf('poptNG%dNu%d%sSNR%02d.mat',NGauss,NumberUncertain,myoptions.Algorithm,modelSNR) ,'popt','params','Mxy','Mz','Mxyopt','Mzopt','signu','signuImage')
      handle = figure(10)
      plot(params.TRList,Mxyopt(1,:),'b',params.TRList,Mxyopt(2,:),'k')
      ylabel('MI Mxy')
      xlabel('sec'); legend('Pyr','Lac')
      saveas(handle,sprintf('OptMxyNG%dNu%d%sSNR%02d',NGauss,NumberUncertain,myoptions.Algorithm,modelSNR),'png')
      handle = figure(11)
      plot(params.TRList,popt.FaList(1,:)*180/pi,'b',params.TRList,popt.FaList(2,:)*180/pi,'k')
      ylabel('MI FA (deg)')
      xlabel('sec'); legend('Pyr','Lac')
      saveas(handle,sprintf('OptFANG%dNu%d%sSNR%02d',NGauss,NumberUncertain,myoptions.Algorithm,modelSNR),'png')
      handle = figure(12)
      plot(params.TRList,Mzopt(1,:),'b--',params.TRList,Mzopt(2,:),'k--')
      hold
      plot(params.TRList,popt.state(:,1, 1),'b',params.TRList,popt.state(:,2, 1),'k')
      if(lqp > 1)
        plot(params.TRList,popt.state(:,1, 5),'b',params.TRList,popt.state(:,2, 5),'k')
        plot(params.TRList,popt.state(:,1,10),'b',params.TRList,popt.state(:,2,10),'k')
        plot(params.TRList,popt.state(:,1,15),'b',params.TRList,popt.state(:,2,15),'k')
      end
      ylabel('MI Mz ')
      xlabel('sec'); legend('Pyr','Lac')
      saveas(handle,sprintf('OptMzNG%dNu%d%sSNR%02d',NGauss,NumberUncertain,myoptions.Algorithm,modelSNR),'png')
  end 


end
