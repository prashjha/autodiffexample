clear all
close all
clc


%% Choose Excitation Angle
FAType = {'OED','Const'};
idoed = 1;
idcst = 2;

% load synthetic data
Nspecies = 2
Ntime = 23

% extract timehistory info
num_trials = 25;

TR = 2;
TR_list = (0:(Ntime-1))*TR;

% plot gamma
jmA0    = 1.
jmalpha = 2.5
jmbeta  = 4.5
jmt0    = 4
%jmaif   = jmA0 * (TR_list- jmt0  ).^jmalpha .* exp(-(TR_list- jmt0 )/jmbeta);
jmaif   = jmA0  * gampdf(TR_list - jmt0  , jmalpha , jmbeta);
figure(5)
plot(TR_list,jmaif ,'b')
ylabel('aif')
xlabel('sec')

% init params
initT1a = 30;
initT1b = 25;
initKpl = 0.15;
kve = 0.05;
ve = 1.00;
VIF_scale_fact = [jmA0;0];
opts = optimset('lsqcurvefit');
opts.TolFun = 1e-09;
opts.TolX = 1e-09;
opts.Display = 'off';

%for idesign = 1:numel(FAType)
for idesign = 2:2
    switch (FAType{idesign})
    case('Const') 
            params = struct('t0',[jmt0;0],'gammaPdfA',[jmalpha;1],'gammaPdfB',[jmbeta;1],...
                'scaleFactor',VIF_scale_fact,'T1s',[initT1a,initT1b],'ExchangeTerms',[0,initKpl;0,0],...
                'TRList',TR_list,'PerfusionTerms',[kve,0],'volumeFractions',ve,...
                'fitOptions', opts);
            params.FaList = 20*pi/180*ones(2,Ntime);
            params.FaList(2,:) = 15*pi/180*ones(1,Ntime);
    case('OED') 
            % load oed results
            params = struct('t0',[jmt0;0],'gammaPdfA',[jmalpha;1],'gammaPdfB',[jmbeta;1],...
                'scaleFactor',VIF_scale_fact,'T1s',[initT1a,initT1b],'ExchangeTerms',[0,initKpl;0,0],...
                'TRList',prashantoedresults.opt_res.TRList_opt,'PerfusionTerms',[kve,0],'volumeFractions',ve,...
                'fitOptions', opts);
            params.FaList = prashantoedresults.opt_res.FA_opt;
    end
    %% evaluate walker model
    model = HPKinetics.NewMultiPoolTofftsGammaVIF();
    M0 = [0,0];
    [t_axis,walkerMxy,walkerMz] = model.compile(M0.',params);
    
    % setup optimization variables
    kpl = initKpl;
    kveqp =   kve/ ve ;
    T1P = initT1a;
    T1L = initT1b;
    A0  = jmA0;
    
    statevariable  = zeros([Nspecies,Ntime]);


    disp('build constraints')
    disp('expm not available for AD.... ')
    disp('TODO: verify matrix exponent impelmented with sylvester formula')
    statevariable(:,1 )=0;
    for iii = 1:Ntime-1
      disp([Nspecies*(iii-1)+1 , Nspecies*(iii-1)+2 ,Ntime* Nspecies])
      klpqp =    0 ;     % @cmwalker where do I get this from ? 
      %A = [-1/T1P - kpl - kveqp,  klpqp; kpl, -1/T1L - klpqp];
      %[V,D] = eig(A);
      %traceA = A(1,1) + A(2,2);
      %detA = A(1,1)*A(2,2) - A(1,2)*A(1,2)  ;
      %lambda1 =traceA  + sqrt(traceA*traceA - 4 * detA) ;
      %lambda2 =traceA  - sqrt(traceA*traceA - 4 * detA) ;
      %V = [ A(1,2), A(2,2) - lambda2 ; A(1,1) - lambda1 ,-A(2,1)]; 
      %Vinv = 1/(V(2,2)* V(1,1)  -V(2,1)* V(1,2))*[V(2,2), -V(2,1); -V(1,2), V(1,1)];
      %expATR = V*[exp(lambda1 * TR) 0 ; 0 exp(lambda2 * TR)]*Vinv;
      % 
      currentTR = params.TRList(iii+1) - params.TRList(iii);
      nsubstep = 5;
      deltat = currentTR /nsubstep ;
      integratedt = [params.TRList(iii):deltat:params.TRList(iii+1)] +deltat/2  ;
      integrand = A0 *gampdf(integratedt(1:nsubstep )'-jmt0,jmalpha,jmbeta) ;
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
      expATR = [ exp(-currentTR*(kpl + kveqp + 1/T1P)),                   0; (kpl*exp(-currentTR/T1L) - kpl*exp(-currentTR*(kpl + kveqp + 1/T1P)))/(kpl + kveqp - 1/T1L + 1/T1P), exp(-currentTR/T1L)];
      % mid-point rule integration
      aifterm = kveqp * deltat * [ exp((-1/T1P - kpl - kveqp)*deltat*[.5:1:nsubstep] );
    (kpl*exp((-1/T1P - kpl - kveqp)*deltat*[.5:1:nsubstep] ) - kpl*exp(-1/T1L *deltat*[.5:1:nsubstep] ))/((-1/T1P - kpl - kveqp) + 1/T1L )] * integrand ;
      statevariable(:,iii+1) =  expATR *( statevariable(:,iii ))   + aifterm     ;
      statevariable(:,iii+1) =  cos(params.FaList(:,iii+1)).* statevariable(:,iii+1);
      %statevariable(:,iii+1) =  expATR *( cos(params.FaList(:,iii)).* statevariable(:,iii ))   + aifterm     ;
    end
end 
    figure(3)
    plot(params.TRList,walkerMz(1,:),'b--',params.TRList,walkerMz(2,:),'k--')
    hold
    plot(params.TRList,statevariable(1,:),'b',params.TRList,statevariable(2,:),'k')
