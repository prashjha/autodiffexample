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

% quad points
NGauss = 5;

% noise for data distribution
SignalNoiseMI = 10;

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
    'fitOptions', opts, ...
    'NumQP', NGauss, 'SignalNoiseMI', SignalNoiseMI);
model = HPKinetics.NewMultiPoolTofftsGammaVIF();


%% Get true Mz
%% Choose Excitation Angle
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
            params.FaList = flips;
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
    plot(TR_list,flips(1,:)*180/pi,'b',TR_list,flips(2,:)*180/pi,'k')
    ylabel('Const FA (deg) ')
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
    % setup optimization variables
    Nspecies = 2;
    FaList = optimvar('FaList',Nspecies,Ntime,'LowerBound',0, 'UpperBound',35*pi/180);
    TRList = TR_list;
    diffTR = diff(TRList);

    signu = SignalNoiseMI ; % TODO - FIXME
    [x2,xn2,xm2,w2,wn2]=GaussHermiteNDGauss(NGauss,0,signu);
    lqp2=length(xn2{1}(:));

    NumberUncertain = 3;
    switch (NumberUncertain)
       case(3)
         [x,xn,xm,w,wn]=GaussHermiteNDGauss(NGauss,[tisinput(5:2:9)],[tisinput(6:2:10)]);
         T1Pqp   = T1pmean + zeros(size(xn{1}(:)));
         T1Lqp   = T1lmean + zeros(size(xn{1}(:)));
         kplqp   = xn{1}(:);
         klpqp   =    0 + zeros(size(xn{1}(:)));
         kveqp   = xn{2}(:);
         t0qp    = xn{3}(:); 
       case(4)
         [x,xn,xm,w,wn]=GaussHermiteNDGauss(NGauss,[tisinput(1:2:7)],[tisinput(2:2:8)]);
         T1Pqp   = xn{1}(:);
         T1Lqp   = xn{2}(:);
         kplqp   = xn{3}(:);
         klpqp   =    0 + zeros(size(xn{1}(:)));
         kveqp   = xn{4}(:);
         t0qp    = t0mean(1) + zeros(size(xn{1}(:)));
       case(5)
         [x,xn,xm,w,wn]=GaussHermiteNDGauss(NGauss,[tisinput(1:2:9)],[tisinput(2:2:10)]);
         T1Pqp   = xn{1}(:);
         T1Lqp   = xn{2}(:);
         kplqp   = xn{3}(:);
         klpqp   =    0 + zeros(size(xn{1}(:)));
         kveqp   = xn{4}(:);
         t0qp    = xn{5}(:); 
    end

    lqp=length(xn{1}(:));

    %% optimizing
    % compute TR's from TR_list
    TRi = getTR(TR_list); % vector of 'size(TR_list) - 1'

    InitialGuessType = FAType{1};
    InitialGuess =  [flips(:)];   
    % Pulse Sequence Bounds
    pmin =  [flips(:)*0];     
    pmax =  [flips(:)*0+35*pi/180];
    findiffrelstep=1.e-6;
    tolx=1.e-9;%1.e-5;
    tolfun=1.e-9;%1.e-5;QALAS_synphan_MIcalc.m
    maxiter=400;

    % fn
    Fx = @(x) MI_GHQuad_HPTofts_With_Der_Parallel_old_method(x, M0, params, model, NumberUncertain, xn, wn, xn2, wn2, T1Pqp, T1Lqp, kplqp, klpqp, kveqp, t0qp);

    tic;
    [popt,fval,exitflag,output,lambda,grad,hessian] ...
     =fmincon(Fx, InitialGuess ,[],[],[],[],pmin,pmax,[],...
        optimset('TolX',tolx,'TolFun',tolfun,'MaxIter', ...
        maxiter,'Display','iter-detailed',... 
        'GradObj','on','PlotFcn',{'optimplotfvalconstr', 'optimplotconstrviolation', 'optimplotfirstorderopt' }));
    toc;

    % save convergence history
    handle = figure(5)
    saveas(handle,sprintf('historyNG%dNu%dadjSNR%02d',NGauss,NumberUncertain,SignalNoiseMI ),'png')
    params.FaList = reshape(popt(1:end),size(params.FaList ));
    save(sprintf('poptNG%dNu%dadjSNR%02d.mat',NGauss,NumberUncertain,SignalNoiseMI) ,'params')
    [t_axisopt,Mxyopt,Mzopt] = model.compile(M0.',params);
    handle = figure(10)
    plot(params.TRList,Mxyopt(1,:),'b',params.TRList,Mxyopt(2,:),'k')
    ylabel('adj MI Mxy')
    xlabel('sec'); legend('Pyr','Lac')
    saveas(handle,sprintf('OptMxyNG%dNu%dadjSNR%02d',NGauss,NumberUncertain,SignalNoiseMI),'png')
    handle = figure(11)
    plot(params.TRList,params.FaList(1,:)*180/pi,'b',params.TRList,params.FaList(2,:)*180/pi,'k')
    ylabel('adj MI FA (deg)')
    xlabel('sec'); legend('Pyr','Lac')
    saveas(handle,sprintf('OptFANG%dNu%dadjSNR%02d',NGauss,NumberUncertain,SignalNoiseMI),'png')
    handle = figure(12)
    plot(params.TRList,Mzopt(1,:),'b',params.TRList,Mzopt(2,:),'k')
    ylabel('adj MI Mz ')
    xlabel('sec'); legend('Pyr','Lac')
    saveas(handle,sprintf('OptMzNG%dNu%dadjSNR%02d',NGauss,NumberUncertain,SignalNoiseMI),'png')
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

function t = getTime(TR)
% compute time sequence from TR
N = size(TR,2) + 1;
t = zeros(1,N);
for i=2:N
    t(i) = 0;
    for j=1:(i-1)
        t(i) = t(i) + TR(j);
    end
end
end

%% to save the function values during optimization iteration
function stop = outfun(x,optimValues,state)
stop = false;
phistory.fval=[];
phistory.x=[];
phistory.iteration=[];
grad=[];
figmarkstr={'bo','gx','r+','c*','ms','yd','kv','w^'};
switch state
    case 'init'
        disp('init')
    case 'iter'
        % Concatenate current point and objective function
        % value with history. x must be a row vector.
        phistory.fval = [phistory.fval; optimValues.fval];
        phistory.x = [phistory.x; x];
        phistory.iteration = [phistory.iteration; optimValues.iteration];
        % Concatenate current search direction with
        % searchdir.
        grad = [grad;...
            optimValues.gradient'];
    case 'done'
        disp('done')
        save('opt_hisotry.mat', 'phistory')
    otherwise
end
end

function stop = outfun2(x,optimValues,state,fileID)
    stop = false;
    fprintf(fileID,'%u,%u,%16.14e \r\n',optimValues.iteration,optimValues.funccount, optimValues.fval);
end



