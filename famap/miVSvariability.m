

pre1img = load('20200482-006-Pre-1.mat')
pre2img = load('20200482-006-Pre-2.mat')

maskdata = niftiread('isofa.nii');
plot( pre1img.kpla(:), pre2img.kpla(:),'x' )

midata02 = load('brutesearchNG5Nu3constDirectTotalSignalSNR20Hermite.mat')

getMIdata = zeros(size(pyruvateFA));
for iii =1:length(pyruvateFA(:)  ) 
    faindex = find( midata02.pyrgrid(:)== round(pyruvateFA(iii))  & midata02.lacgrid(:)== round(lactateFA(iii)));
    if(~isempty(faindex) )
        getMIdata(iii) = midata02.brutesearch(faindex);
    end
end

splitdataMI = {}
splitdatakpl = {}
mivalues= nan(12,1)
varvalues= nan(12,1)
for iii = 1:12
    myindexone =  find(maskdata(:) ==iii & pre1img.masks(:)==1 );
    myindextwo =  find(maskdata(:) ==iii & pre2img.masks(:)==1 );
    splitdataMI{iii} = [getMIdata(myindexone );getMIdata(myindextwo )];
    splitdatakpl{iii}= [pre1img.kpla(myindexone );pre2img.kpla(myindextwo )];
    if(~isempty(splitdataMI{iii} ))
      mivalues(iii)=  mean(splitdataMI{iii});
      varvalues(iii)=  mean(splitdatakpl{iii});
    end
end
plot(mivalues,varvalues,'x')

