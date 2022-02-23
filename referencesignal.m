  clear all
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
  Ntime = 30;
  currentTR = 3;
  TR_list = (0:(Ntime-1))*currentTR ;
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
  paramstwo = struct('t0',[t0mean(1);0],'gammaPdfA',[alphamean(1)  ;1],'gammaPdfB',[betamean(1);1],...
      'scaleFactor',VIF_scale_fact,'T1s',[43,33],'ExchangeTerms',[0,.02 ;0,0],...
      'TRList',TR_list ,'PerfusionTerms',[.1,0],'volumeFractions',.90,...
      'fitOptions', opts);
  model = HPKinetics.NewMultiPoolTofftsGammaVIF();
  
  
  %% Choose Excitation Angle
  FAType = {'Const'};
  %% HACK- @cmwalker code for initial conditions - https://github.com/fuentesdt/TumorHPMRI/blob/master/models/gPC/walker/ShowMxyPub.m
  for i = 1:numel(FAType)
      switch (FAType{i})
          case('Const') % Nagashima for lactate const 10 pyruvate
              tic
              E1(1) = exp(-currentTR *(1/T1pmean+kplmean));
              E1(2) = exp(-currentTR /T1lmean);
              for n = 1:Ntime
                  % 20deg for pyruvate 30deg for lactate - currently used in brain
                  flips(2,n) = 30*pi/180;
                  flips(1,n) = 20*pi/180;
              end
              params.FaList = flips ;
              paramstwo.FaList = flips ;
      end
  
      
      
  
      tic
      %% Fitting
      [t_axistwo,Mxytwo,Mztwo] = model.compile(M0.',paramstwo);
      [t_axis,Mxy,Mz] = model.compile(M0.',params);
      toc
  end

     modelSNR=20;
     signuImage = max(Mxy(1,:))/modelSNR
     signuImagetwo = max(Mxytwo(1,:))/modelSNR
     signuImagetwo/signuImage 

