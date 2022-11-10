clear all
close all


hackuncertainList=3
hackgpList=5
hacksolvertype='interior-point'
ObjectiveType = {'TotalSignal','SumQuad'}
snrList = [2,5,10,15,20]
idplot = 0
NGauss = 5,NumberUncertain=3, QuadratureRule = 'Hermite';
myoptions.Algorithm = 'interior-point';
for isnr = 1:length(snrList)
   worktmp(isnr) = load(sprintf('optim_variable_TR_FA_results/poptNG%dNu%d%s%sSNR%02dHermite-OptFAandTR.mat', hackgpList,hackuncertainList,hacksolvertype,ObjectiveType{1},snrList( isnr )) ) 
   idplot = idplot+1
   handle = figure(idplot )
   %plot(params.TRList,Mxyopt(1,:),'b',params.TRList,Mxyopt(2,:),'k')
   plot( worktmp(isnr).popt.TRList, worktmp(isnr).popt.FaList(1,:)*180/pi,'b', worktmp(isnr).popt.TRList, worktmp(isnr).popt.FaList(2,:)*180/pi,'k')
   ylabel('MI FA (deg)')
   xlabel('sec'); legend('Pyr','Lac')
   ylim([0 40])
   set(gca,'FontSize',16)
   saveas(handle,sprintf('OptFANG%dNu%d%sTR%sSNR%02d%s' ,NGauss,NumberUncertain,myoptions.Algorithm,ObjectiveType{1},snrList( isnr ),QuadratureRule),'png')
   
   idplot = idplot+1
   handle = figure(idplot )
   M0 = [0,0];
   worktmp(isnr).params.TRList = worktmp(isnr).popt.TRList;
   worktmp(isnr).params.FaList = worktmp(isnr).popt.FaList;

   model = HPKinetics.NewMultiPoolTofftsGammaVIF();
   [t_axisopt,Mxyopt,Mzopt] = model.compile(M0.',worktmp(isnr).params);
   plot( worktmp(isnr).popt.TRList, Mxyopt(1,:),'b', worktmp(isnr).popt.TRList,Mxyopt(2,:),'k' )
   ylabel('MI Mxy')
   xlabel('sec'); legend('Pyr','Lac')
   set(gca,'FontSize',16)
   saveas(handle,sprintf('OptMxyNG%dNu%d%sTR%sSNR%02d%s',NGauss,NumberUncertain,myoptions.Algorithm,ObjectiveType{1},snrList( isnr ),QuadratureRule),'png')
end


 
