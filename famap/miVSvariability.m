

pre1img = load('20200482-006-Pre-1.mat')
pre2img = load('20200482-006-Pre-2.mat')

midata02 = load('brutesearchNG5Nu3constDirectTotalSignalSNR02Hermite.mat')

maskdata = niftiread('isofa.nii');

plot( pre1img.kpla(:), pre2img.kpla(:),'x' )

getMIdata = zeros(size(pyruvateFA));
for iii =1:length(pyruvateFA(:)  ) 
    faindex = find( midata02.pyrgrid(:)== round(pyruvateFA(iii))  & midata02.lacgrid(:)== round(lactateFA(iii)));
    if(~isempty(faindex) )
        getMIdata(iii) = midata02.brutesearch(faindex);
    end
end

splitdataMI = {}
splitdatakpl = {}
for iii = 1:12
    myindex =  find(maskdata(:) ==iii);
    splitdataMI{iii} = [getMIdata(myindex );getMIdata(myindex )];
    splitdatakpl{iii}= [pre1img.kpla(myindex );pre2img.kpla(myindex )];
end
