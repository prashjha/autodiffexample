
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
tmpMIdata = zeros(size(pyruvateFA));
for iii =1:length(pyruvateFA(:)  ) 
    faindex = find( midata.pyrgrid(:)== round(pyruvateFA(iii))  & midata.lacgrid(:)== round(lactateFA(iii)));
    if(~isempty(faindex) )
        tmpMIdata(iii) = midata.brutesearch(faindex);
    end
end
niftiwrite(tmpMIdata,sprintf('miimage%02d.nii',SNRList(kkk)));
getMIdata{kkk} = tmpMIdata;
end

minlac = min(round(lactateFA(:)))
maxlac = max(round(lactateFA(:)))
minpyr = min(round(pyruvateFA(:)))
maxpyr = max(round(pyruvateFA(:)))
mivalues = {nan(length( minlac:maxlac),length( minlac:maxlac)),nan(length( minlac:maxlac),length( minlac:maxlac)),nan(length( minlac:maxlac),length( minlac:maxlac)),nan(length( minlac:maxlac),length( minlac:maxlac)),nan(length( minlac:maxlac),length( minlac:maxlac)),};
varvalues= {nan(length( minlac:maxlac),length( minlac:maxlac)),nan(length( minlac:maxlac),length( minlac:maxlac)),nan(length( minlac:maxlac),length( minlac:maxlac)),nan(length( minlac:maxlac),length( minlac:maxlac)),nan(length( minlac:maxlac),length( minlac:maxlac)),};
for iii = minlac:maxlac
  for jjj = minpyr:maxpyr
    myindexone =  find(round(lactateFA(:))==iii & round(pyruvateFA(:)) == jjj & pre1img.masks(:)==1 );
    myindextwo =  find(round(lactateFA(:))==iii & round(pyruvateFA(:)) == jjj & pre2img.masks(:)==1 );
    splitdatakpl= [pre1img.kpla(myindexone );pre2img.kpla(myindextwo )];
    for kkk = 1:length(SNRList)
      splitdataMI = [getMIdata{kkk}(myindexone );getMIdata{kkk}(myindextwo )];
      if(~isempty(splitdataMI ))
        disp(sprintf('%d %d %d',iii,jjj,kkk))
        mivalues{kkk}(iii,jjj) = mean(splitdataMI);
        varvalues{kkk}(iii,jjj)= var(splitdatakpl);
      end
    end
  end
end
plot(mivalues{1}(:),varvalues{1}(:),'rx',...
     mivalues{2}(:),varvalues{2}(:),'gx',...
     mivalues{3}(:),varvalues{3}(:),'bx',...
     mivalues{4}(:),varvalues{4}(:),'kx',...
     mivalues{5}(:),varvalues{5}(:),'cx')
xlabel('MI')
ylabel('var')

