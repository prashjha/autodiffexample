clear all
close all
clc


solverType = {'adj'};
solverType = {'constDirect','sqp','interior-point'};
ObjectiveType = {'TotalSignal'}
gpList = [3, 4,5]
gpList = [3]
uncertainList = [3]
snrList = [2,10]
snrList = [2,10,25]
snrList = [2,10,20,25]
numsolves = numel(solverType) * length(gpList) * length(uncertainList) * length(snrList) + length(snrList)
solnList(numsolves) = struct('gp',[],'snr',[],'numberuncertain',[],'FaList',[],'solver',[], 'params', [], 'Mxy', [], 'Mz', [],'signuImage',[],'signu',[]);
icount  = 0;
for isolver = 1:numel(solverType)
 for igp = 1:length(gpList)
  for inu = 1:length(uncertainList)
   for isnr = 1:length(snrList)
      worktmp = load(sprintf('poptNG%dNu%d%s%sSNR%02d.mat',gpList(igp),uncertainList(inu),solverType{isolver},ObjectiveType{1},snrList(isnr)));
      icount= icount+1;
      solnList (icount) = struct('gp',gpList(igp),'snr',snrList(isnr),'numberuncertain',uncertainList(inu),'FaList',worktmp.popt.FaList,'solver',solverType{isolver},'params',worktmp.params, 'Mxy',worktmp.Mxyref, 'Mz',worktmp.Mzref,'signuImage',worktmp.signuImage,'signu',worktmp.signu);
   end
  end
 end
end

% compute variance for each SNR
for isnr = 1:length(snrList)
   icount= icount+1;
   solnList (icount) = struct('gp',-1,'snr',snrList(isnr),'numberuncertain',-1,'FaList',worktmp.params.FaList,'solver','const','params',worktmp.params, 'Mxy',worktmp.Mxy, 'Mz',worktmp.Mz,'signuImage',solnList(isnr).signuImage,'signu',solnList(isnr).signu);
end


% extract timehistory info
num_trials = 25;
Ntime    = size(solnList(1).Mz,2);
Nspecies = size(solnList(1).Mz,1);
timehistory  = zeros(Ntime,Nspecies,num_trials+1,length(solnList) );

for jjj =1:length(solnList)
  % NOTE - image noise is at the single image for single species - signu is for the sum over time for both species ==> divide by Ntime and Nspecies
  imagenoise = solnList(jjj).signuImage;
  %disp([xroi(jjj),yroi(jjj),zroi(jjj)]);
  timehistory(:,1,num_trials+1,jjj) = solnList(jjj).Mz(1,:);
  timehistory(:,2,num_trials+1,jjj) = solnList(jjj).Mz(2,:);
  % add noise for num_trials
  for kkk = 1:num_trials
      pyrnoise = imagenoise *randn(Ntime,1);
      timehistory(:,1,kkk,jjj)= timehistory(:,1,num_trials+1,jjj) + pyrnoise;
      lacnoise = imagenoise *randn(Ntime,1);
      timehistory(:,2,kkk,jjj)= timehistory(:,2,num_trials+1,jjj) + lacnoise;
  end
end

TR_list = solnList(1).params.TRList 
handle = figure(1)
plot( TR_list, timehistory(:,1,num_trials+1,2), 'b', ...
      TR_list, timehistory(:,2,num_trials+1,2), 'k', ...
      TR_list, timehistory(:,1,num_trials,2), 'b--', ...
      TR_list, timehistory(:,2,num_trials,2), 'k--')
ylabel('MI Mz ')
xlabel('sec'); legend('Pyr','Lac')
set(gca,'FontSize',16)
saveas(handle,'examplenoise','png')

