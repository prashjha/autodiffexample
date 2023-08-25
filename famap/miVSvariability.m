
pre1img = load('20200482-006-Pre-1.mat')
pre2img = load('20200482-006-Pre-2.mat')

niftiwrite(double(pre1img.masks),'pre1mask.nii');
niftiwrite(double(pre2img.masks),'pre2mask.nii');
niftiwrite(double(pre1img.kpla),'pre1kpla.nii');
niftiwrite(double(pre2img.kpla),'pre2kpla.nii');
plot( pre1img.kpla(:), pre2img.kpla(:),'x' )


SNRList = [ 2 5 10 15 20];
getMIdata = {};
for kkk = 1:length(SNRList)
midata = load(sprintf('../brutesearchNG5Nu3constDirectTotalSignalSNR%02dHermite.mat',SNRList(kkk)));
tmpMIdata = nan(size(pyruvateFA));
for iii =1:length(pyruvateFA(:)  ) 
    faindex = find( midata.pyrgrid(:)== round(pyruvateFA(iii))  & midata.lacgrid(:)== round(lactateFA(iii)));
    if(~isempty(faindex) )
        tmpMIdata(iii) = midata.brutesearch(faindex);
    end
end
niftiwrite(tmpMIdata,sprintf('miimage%02d.nii',SNRList(kkk)));
getMIdata{kkk} = tmpMIdata;
end

lacFA = [min(round(lactateFA(:))):max(round(lactateFA(:)))];
pyrFA = [min(round(pyruvateFA(:))):max(round(pyruvateFA(:)))];
mivalues = {nan(length( lacFA),length( pyrFA)), nan(length( lacFA),length( pyrFA)), nan(length( lacFA),length( pyrFA)), nan(length( lacFA),length( pyrFA)), nan(length( lacFA),length( pyrFA))};
varvalues= {nan(length( lacFA),length( pyrFA)), nan(length( lacFA),length( pyrFA)), nan(length( lacFA),length( pyrFA)), nan(length( lacFA),length( pyrFA)), nan(length( lacFA),length( pyrFA))};
    for kkk = 1:length(SNRList)
for iii = 1:length(lacFA)
  for jjj = 1:length(pyrFA)
    myindexone =  find(round(lactateFA(:))==lacFA(iii) & round(pyruvateFA(:)) == pyrFA(jjj) & pre1img.masks(:)==1 );
    myindextwo =  find(round(lactateFA(:))==lacFA(iii) & round(pyruvateFA(:)) == pyrFA(jjj) & pre2img.masks(:)==1 );
    splitdatakpl= [pre1img.kpla(myindexone );pre2img.kpla(myindextwo )];
      splitdataMI = [getMIdata{kkk}(myindexone );getMIdata{kkk}(myindextwo )];
      %if(~isempty(splitdataMI ))
      if( length(splitdataMI ) > 4)
        disp(sprintf('%d %d %d',lacFA(iii),pyrFA(jjj),SNRList(kkk)))
        myindexone
        myindextwo
        splitdataMI
        splitdatakpl
        var(splitdatakpl)
        mivalues{kkk}(iii,jjj) = mean(splitdataMI);
        varvalues{kkk}(iii,jjj)= var(splitdatakpl);
      end
    end
  end
end
%coercoef(mivalues{1}(:),varvalues{1}(:))
handle = figure(1);plot(mivalues{1}(:),varvalues{1}(:),'rx');xlabel('MI');ylabel('var'); title(sprintf('SNR%02d',SNRList(1)));set(gca,'FontSize',16);lsline;saveas(handle,sprintf('kplvarSNR%02d',SNRList(1)),'png')
handle = figure(2);plot(mivalues{2}(:),varvalues{2}(:),'gx');xlabel('MI');ylabel('var'); title(sprintf('SNR%02d',SNRList(2)));set(gca,'FontSize',16);lsline;saveas(handle,sprintf('kplvarSNR%02d',SNRList(2)),'png')
handle = figure(3);plot(mivalues{3}(:),varvalues{3}(:),'bx');xlabel('MI');ylabel('var'); title(sprintf('SNR%02d',SNRList(3)));set(gca,'FontSize',16);lsline;saveas(handle,sprintf('kplvarSNR%02d',SNRList(3)),'png')
handle = figure(4);plot(mivalues{4}(:),varvalues{4}(:),'kx');xlabel('MI');ylabel('var'); title(sprintf('SNR%02d',SNRList(4)));set(gca,'FontSize',16);lsline;saveas(handle,sprintf('kplvarSNR%02d',SNRList(4)),'png')
handle = figure(5);plot(mivalues{5}(:),varvalues{5}(:),'mx');xlabel('MI');ylabel('var'); title(sprintf('SNR%02d',SNRList(5)));set(gca,'FontSize',16);lsline;saveas(handle,sprintf('kplvarSNR%02d',SNRList(5)),'png')






