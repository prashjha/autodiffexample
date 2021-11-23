clear all
close all
clc


solverType = {'adj'};
solverType = {'adj','sqp','interior-point'};
gpList = [3, 4,5]
gpList = [4]
uncertainList = [3]
snrList = [2,10]
snrList = [2,10,25]
numsolves = numel(solverType) * length(gpList) * length(uncertainList) * length(snrList) +1
solnList(numsolves) = struct('gp',[],'snr',[],'numberuncertain',[],'FaList',[],'solver',[], 'params', [], 'Mxy', [], 'Mz', [],'signuImage',[],'signu',[]);
icount  = 0;
for isolver = 1:numel(solverType)
 for igp = 1:length(gpList)
  for inu = 1:length(uncertainList)
   for isnr = 1:length(snrList)
      worktmp = load(sprintf('poptNG%dNu%d%sSNR%02d.mat',gpList(igp),uncertainList(inu),solverType{isolver},snrList(isnr)));
      icount= icount+1;
      solnList (icount) = struct('gp',gpList(igp),'snr',snrList(isnr),'numberuncertain',uncertainList(inu),'FaList',worktmp.popt.FaList,'solver',solverType{isolver},'params',worktmp.params, 'Mxy',worktmp.Mxyopt, 'Mz',worktmp.Mzopt,'signuImage',worktmp.signuImage,'signu',worktmp.signu);
   end
  end
 end
end

solnList (numsolves) = struct('gp',-1,'snr',-1,'numberuncertain',-1,'FaList',worktmp.params.FaList,'solver','const','params',worktmp.params, 'Mxy',worktmp.Mxy, 'Mz',worktmp.Mz,'signuImage',worktmp.signuImage,'signu',worktmp.signu);

% extract timehistory info
num_trials = 25;
Ntime = 23;
Nspecies = 2;
timehistory  = zeros(Ntime,Nspecies,num_trials+1,length(solnList) );

for jjj =1:length(solnList)
  % NOTE - image noise is at the single image for single species - signu is for the sum over time for both species ==> divide by Ntime and Nspecies
  imagenoise = solnList(jjj).signuImage;
  %disp([xroi(jjj),yroi(jjj),zroi(jjj)]);
  timehistory(:,1,num_trials+1,jjj) = solnList(jjj).Mz(1,:);
  timehistory(:,2,num_trials+1,jjj) = solnList(jjj).Mz(2,:);
  % add noise for num_trials
  for kkk = 1:num_trials
      realchannel = imagenoise *randn(Ntime,1);
      imagchannel = imagenoise *randn(Ntime,1);
      timehistory(:,1,kkk,jjj)= timehistory(:,1,num_trials+1,jjj) + sqrt( realchannel.^2 + imagchannel.^2); 
      realchannel = imagenoise *randn(Ntime,1);
      imagchannel = imagenoise *randn(Ntime,1);
      timehistory(:,2,kkk,jjj)= timehistory(:,2,num_trials+1,jjj) + sqrt( realchannel.^2 + imagchannel.^2);  
  end
end

%% TR_list = solnList(1).params.TRList 
%% figure(1)
%% plot( TR_list, timehistory(:,1,1,1), 'b', ...
%%       TR_list, timehistory(:,2,1,1), 'b-.', ...
%%       TR_list, timehistory(:,1,1,2), 'g', ...
%%       TR_list, timehistory(:,2,1,2), 'g-.', ...
%%       TR_list, timehistory(:,1,1,3), 'r', ...
%%       TR_list, timehistory(:,2,1,3), 'r-.', ...
%%       TR_list, timehistory(:,1,1,4), 'c', ...
%%       TR_list, timehistory(:,2,1,4), 'c-.')
%% title('noise added')
%% figure(2)
%% plot( TR_list, timehistory(:,1,num_trials+1,1), 'b', ...
%%       TR_list, timehistory(:,2,num_trials+1,1), 'b-.', ...
%%       TR_list, timehistory(:,1,num_trials+1,2), 'g', ...
%%       TR_list, timehistory(:,2,num_trials+1,2), 'g-.', ...
%%       TR_list, timehistory(:,1,num_trials+1,3), 'r', ...
%%       TR_list, timehistory(:,2,num_trials+1,3), 'r-.', ...
%%       TR_list, timehistory(:,1,num_trials+1,4), 'c', ...
%%       TR_list, timehistory(:,2,num_trials+1,4), 'c-.')
%% title('original data')

%% store soln
storekplopt   = zeros(num_trials+1,length(solnList));
storekveqpopt = zeros(num_trials+1,length(solnList));
storeT1Popt   = zeros(num_trials+1,length(solnList));
storeT1Lopt   = zeros(num_trials+1,length(solnList));
storeA0opt    = zeros(num_trials+1,length(solnList));
for idesign = 1:length(solnList)
   % setup optimization variables
   numberParameters = 1
   switch (numberParameters)
       case(1) 
         kpl   = optimvar('kpl','LowerBound',0);
         kveqp = solnList(idesign ).params.PerfusionTerms(1) / solnList(idesign ).params.volumeFractions;
         T1P   = solnList(idesign ).params.T1s(1);
         T1L   = solnList(idesign ).params.T1s(2);
         A0    = solnList(idesign ).params.scaleFactor(1);
         % ground truth 
         xstar.kpl        = solnList(idesign).params.ExchangeTerms(1,2);
       case(2) 
         kpl   = optimvar('kpl','LowerBound',0);
         kveqp = optimvar('kveqp','LowerBound',0);
         T1P   = solnList(idesign ).params.T1s(1);        
         T1L   = solnList(idesign ).params.T1s(2);        
         A0    = solnList(idesign ).params.scaleFactor(1);
         % ground truth 
         xstar.kpl        = solnList(idesign).params.ExchangeTerms(1,2); 
         xstar.kveqp      = solnList(idesign ).params.PerfusionTerms(1) / solnList(idesign ).params.volumeFractions; 
       case(4) 
         kpl   = optimvar('kpl','LowerBound',0);
         kveqp = optimvar('kveqp','LowerBound',0);
         T1P   = optimvar('T1P','LowerBound',0);
         T1L   = optimvar('T1L','LowerBound',0);
         A0    = solnList(idesign ).params.scaleFactor(1);
         % ground truth 
         xstar.kpl        = solnList(idesign).params.ExchangeTerms(1,2);
         xstar.kveqp      = solnList(idesign ).params.PerfusionTerms(1) / solnList(idesign ).params.volumeFractions;
         xstar.T1P        = solnList(idesign ).params.T1s(1);
         xstar.T1L        = solnList(idesign ).params.T1s(2);
       case(5) 
         kpl = optimvar('kpl','LowerBound',0);
         kveqp = optimvar('kveqp','LowerBound',0);
         T1P = optimvar('T1P','LowerBound',0);
         T1L = optimvar('T1L','LowerBound',0);
         A0  = optimvar('A0' ,'LowerBound',0);
         % ground truth 
         xstar.kpl        = solnList(idesign).params.ExchangeTerms(1,2); 
         xstar.kveqp      = solnList(idesign ).params.PerfusionTerms(1) / solnList(idesign ).params.volumeFractions;
         xstar.T1P        = solnList(idesign ).params.T1s(1);        
         xstar.T1L        = solnList(idesign ).params.T1s(2);        
         xstar.A0         = solnList(idesign ).params.scaleFactor(1);
   end
   statevariable  = optimexpr([Nspecies,Ntime]);


   disp('build constraints')
   statevariable(:,1 )=0;
   for iii = 1:Ntime-1
     disp([Nspecies*(iii-1)+1 , Nspecies*(iii-1)+2 ,Ntime* Nspecies])
     klpqp =    0 ;     % @cmwalker where do I get this from ? 
     % 
     currentTR = solnList(idesign ).params.TRList(iii+1) - solnList(idesign ).params.TRList(iii);
     nsubstep = 5;
     deltat = currentTR /nsubstep ;
     integratedt = [solnList(idesign ).params.TRList(iii):deltat:solnList(idesign ).params.TRList(iii+1)] +deltat/2  ;
     integrand = A0 *gampdf(integratedt(1:nsubstep )'-solnList(idesign ).params.t0(1),solnList(idesign ).params.gammaPdfA(1),solnList(idesign ).params.gammaPdfB(1)) ;
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
     statevariable(:,iii+1) =  cos(solnList(idesign ).FaList(:,iii+1)).* statevariable(:,iii+1);
   end
   truthstate = evaluate(statevariable ,xstar);

   %% 
   % Create an optimization problem using these converted optimization expressions.
   disp('build objective function')
   %  for idtrial = num_trials+1:num_trials+1
   for idtrial = 1:num_trials+1
       tic;
       disp(sprintf('design = %d, trial = %d',idesign, idtrial )); solnList(idesign)
       mycostfcn = sum( (statevariable(1,:)'- timehistory(:,1,idtrial,idesign ) ).^2) + sum( (statevariable(2,:)'- timehistory(:,2,idtrial,idesign ) ).^2);

       disp('create optim prob')
       convprob = optimproblem('Objective',mycostfcn );
       %% 
       % View the new problem.
       
       %show(convprob)
       % problem = prob2struct(convprob,'ObjectiveFunctionName','generatedObjectiveRecover');
       %% 
       % Solve the new problem. The solution is essentially the same as before.
       

       %myoptions = optimoptions(@fmincon,'Display','iter-detailed','SpecifyObjectiveGradient',true, 'SpecifyConstraintGradient', true)
       %myoptions = optimoptions(@fmincon,'Display','iter-detailed','MaxFunctionEvaluations' , 3.000000e+04)
       myoptions = optimoptions(@lsqnonlin,'Display','iter-detailed');
           
       %[popt,fval,exitflag,output] = solve(convprob,x0,'Options',myoptions, 'ObjectiveDerivative', 'auto-reverse' , 'ConstraintDerivative', 'auto-reverse', 'solver', 'fmincon' )
       %[popt,fval,exitflag,output] = solve(convprob,x0,'Options',myoptions, 'ObjectiveDerivative', 'finite-differences','solver', 'fmincon' )
       % random initial guess
       switch (numberParameters)
          case(1) 
            x0.kpl        = unifrnd(.01,10);
          case(2) 
            x0.kpl        = unifrnd(.01,10);
            x0.kveqp      = unifrnd(.02,10);
          case(4) 
            x0.kpl        = unifrnd(.01,10);
            x0.kveqp      = unifrnd(.02,10);
            x0.T1P        = unifrnd(20 ,40);
            x0.T1L        = unifrnd(15 ,35);
          case(5) 
            x0.kpl        = unifrnd(.01,10);
            x0.kveqp      = unifrnd(.02,10);
            x0.T1P        = unifrnd(20 ,40);
            x0.T1L        = unifrnd(15 ,35);
            x0.A0         = unifrnd(1  ,100);
       end
       initialstate = evaluate(statevariable ,x0);
       [popt,fval,exitflag,output] = solve(convprob,x0,'Options',myoptions, 'solver','lsqnonlin' , 'ObjectiveDerivative', 'finite-differences')
       switch (numberParameters)
          case(1) 
            storekplopt(  idtrial,idesign)  = popt.kpl;
            storekveqpopt(idtrial,idesign)  = solnList(idesign ).params.PerfusionTerms(1) / solnList(idesign ).params.volumeFractions;
            storeT1Popt(  idtrial,idesign)  = solnList(idesign ).params.T1s(1);
            storeT1Lopt(  idtrial,idesign)  = solnList(idesign ).params.T1s(2);
            storeA0opt (  idtrial,idesign)  = solnList(idesign ).params.scaleFactor(1);                                               
          case(2) 
            storekplopt(  idtrial,idesign)  = popt.kpl;
            storekveqpopt(idtrial,idesign)  = popt.kveqp;
            storeT1Popt(  idtrial,idesign)  = solnList(idesign ).params.T1s(1); 
            storeT1Lopt(  idtrial,idesign)  = solnList(idesign ).params.T1s(2); 
            storeA0opt (  idtrial,idesign)  = solnList(idesign ).params.scaleFactor(1);
          case(4) 
            storekplopt(  idtrial,idesign)  = popt.kpl;
            storekveqpopt(idtrial,idesign)  = popt.kveqp;
            storeT1Popt(  idtrial,idesign)  = popt.T1P;
            storeT1Lopt(  idtrial,idesign)  = popt.T1L;
            storeA0opt (  idtrial,idesign)  = solnList(idesign ).params.scaleFactor(1);
          case(5) 
            storekplopt(  idtrial,idesign)  = popt.kpl;
            storekveqpopt(idtrial,idesign)  = popt.kveqp;
            storeT1Popt(  idtrial,idesign)  = popt.T1P;
            storeT1Lopt(  idtrial,idesign)  = popt.T1L;
            storeA0opt (  idtrial,idesign)  = popt.A0 ; 
       end
       slnstate = evaluate(statevariable ,popt);

       if idtrial >= num_trials 
       % plot
       idplot = (num_trials+1)*(idesign-1) + idtrial ;
       figure(idplot )
       plot(solnList(idesign ).params.TRList , timehistory(:,1,idtrial,idesign), 'b', ...
            solnList(idesign ).params.TRList , timehistory(:,2,idtrial,idesign), 'b-.', ...
            solnList(idesign ).params.TRList , timehistory(:,1,num_trials+1,idesign), 'r', ...
            solnList(idesign ).params.TRList , timehistory(:,2,num_trials+1,idesign), 'r-.', ...
            solnList(idesign ).params.TRList , initialstate(1,:), 'g', ...
            solnList(idesign ).params.TRList , initialstate(2,:), 'g-.', ...
            solnList(idesign ).params.TRList , truthstate(1,:), 'm', ...
            solnList(idesign ).params.TRList , truthstate(2,:), 'm-.', ...
            solnList(idesign ).params.TRList , slnstate(1,:), 'k', ...
            solnList(idesign ).params.TRList , slnstate(2,:), 'k-.')
       ylabel('Mz')
       xlabel('sec')
       title('curvefit ')
       legend('walker+rice','','walker','','ic','','truth','','df','')
       pause(.1)
       end 

       toc;
   end 

end 
%%
%%
idplot = (num_trials+1)* length(solnList)+1;
figure(idplot )
labellist = sprintfc('snr%02d',snrList); labellist{end+1} = 'const'
boxplot( [ storekplopt(1:num_trials,1:length(snrList)), storekplopt(1:num_trials,end)], labellist  )

idplot = idplot+1
figure(idplot )
inversevar = var(storekplopt(1:num_trials,:),0,1)
plot( snrList , inversevar(0*length(snrList)+1:1*length(snrList)), 'b',...
      snrList , inversevar(1*length(snrList)+1:2*length(snrList)), 'r',... 
      snrList , inversevar(2*length(snrList)+1:3*length(snrList)), 'k') 
hold
yline(inversevar (end))
legend('walker+rice','','walker','','ic','','truth','','df','')

%%figure(6)
%%disp([ mean(storekplopt(:,:,1),2) var(storekplopt(:,:,1),0,2) mean(storekplopt(:,:,2),2) var(storekplopt(:,:,2),0,2) ])
%%disp([ var(storekplopt(nroipixel:nroipixel+1,1:num_trials,1),0,2) var(storekplopt(nroipixel:nroipixel+1,1:num_trials,2),0,2) ])
    

