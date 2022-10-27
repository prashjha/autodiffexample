clear all; clc; clf;

snr_list = [2, 5, 10, 15, 20];
NGauss = 5;
NumberUncertain = 3;
Ntime = 30; 

% David base data file
f1 = 'results/poptNG5Nu3interior-pointTotalSignalSNR02Hermite/poptNG5Nu3interior-pointTotalSignalSNR02Hermite.mat';

addpath ../.                       % to read new model
addpath ../../autodiffexample_df/. % to read old model

% load model
model = HPModel();
model_old = HPKinetics.NewMultiPoolTofftsGammaVIF();
M0 = [0;0];

% load base data
d1 = load(f1); 
p1 = d1.params;
snr_oed_d1 = 2;
FaList1 = d1.popt.FaList; 
TRList1 = p1.TRList;
TR1 = model.getTR(TRList1);
Mxy1 = d1.Mxyref;
Mz1 = d1.Mzref;
p_new_model1 = p1;
p_new_model1.TRList = TRList1;
p_new_model1.FaList = FaList1;
[t1, Mxy_new_model1, Mz_new_model1] = model.compile(M0, p_new_model1);

% loop over snr and create a copy
for isnr = 1:length(snr_list)
    
    f2_path = sprintf('results/SNR_%d_NGauss_5_NumberUncertain_3_Nscans_30_optim_TA_FA', snr_list(isnr));
    f2 = sprintf('%s/%s.mat', f2_path, 'opt_res');

    % load second data
    d2 = load(f2);
    p2 = d2.params;
    snr_oed_d2 = snr_list(isnr); 
    FaList2 = p2.FaList; 
    TRList2 = p2.TRList;
    TR2 = model.getTR(TRList2);
    Mxy2 = d2.Mxy_opt;
    Mz2 = d2.Mz_opt;
    p_new_model2 = p2;
    p_new_model2.TRList = TRList2;
    p_new_model2.FaList = FaList2;
    [t2, Mxy_new_model2, Mz_new_model2] = model.compile(M0, p_new_model2);

    %% copy
    % create a copy of data 'p1' and update 'TRList' in that copy from data 'p2'
    % and verify that this new copy produces same output as the data 'd2' and 'p2'
    copy_model_par2 = p1;
    copy_model_par2.TRList = TRList2;
    copy_model_par2.FaList = FaList2;
    [copy_model2_t, copy_model2_Mxy, copy_model2_Mz] ...
        = model_old.compile(M0, copy_model_par2);

    copy_model_data2 = d1;
    % update Fa and reference Mxyref and Mzref, keep Mxy and Mz unchanged!
    copy_model_data2.popt.FaList = FaList2;
    copy_model_data2.popt.TRList = TRList2;
    copy_model_data2.Mxyref = copy_model2_Mxy;
    copy_model_data2.Mzref = copy_model2_Mz;
    % update signu and signuImage (code below is how these computed in d1 data)
    % > signuImage = max(Mxy(1,:))/modelSNR;
    % > % variance for Gauss RV is sum. sqrt for std
    % > signu = sqrt(2* Ntime) * signuImage;
    signuImage_d1 = d1.signuImage;
    signu_d1 = d1.signu;
    copy_model_data2.signuImage = (signuImage_d1 * snr_oed_d1) / snr_oed_d2;
    copy_model_data2.signu = sqrt(2 * Ntime) * copy_model_data2.signuImage;

    fval = copy_model_data2.fval;
    popt = copy_model_data2.popt;
    params = copy_model_data2.params;
    Mxy = copy_model_data2.Mxy;
    Mz = copy_model_data2.Mz;
    Mxyref = copy_model_data2.Mxyref;
    Mzref = copy_model_data2.Mzref;
    signu = copy_model_data2.signu;
    signuImage = copy_model_data2.signuImage;
    save(sprintf('%s/poptNG%dNu%d%s%sSNR%02d%s-OptFAandTR.mat', ...
        f2_path, NGauss,NumberUncertain,'interior-point','TotalSignal',...
        snr_oed_d2,'Hermite'),...
            'fval','popt','params','Mxy',...
            'Mz','Mxyref','Mzref','signu','signuImage')


    %% plot
    % verify that p_new_model2 and copy_model2 produce same results
    lw = 1;

    rad_deg_fac = 180./pi;

    figure(isnr)
    sgtitle(sprintf('Compare two different results'))

    subplot(1,3,1);
    plot(copy_model_par2.FaList(1,:)*rad_deg_fac, ...
        'DisplayName', sprintf('%s %s', 'Copy', 'pyr fa'), 'LineWidth', lw)
    hold on
    plot(copy_model_par2.FaList(2,:)*rad_deg_fac, ...
        'DisplayName', sprintf('%s %s', 'Copy', 'lac fa'), 'LineWidth', lw)
    hold on
    plot(p_new_model2.FaList(1,:)*rad_deg_fac, ...
        'DisplayName', sprintf('%s %s', 'Orig', 'pyr fa'), 'LineWidth', lw)
    hold on
    plot(p_new_model2.FaList(2,:)*rad_deg_fac, ...
        'DisplayName', sprintf('%s %s', 'Orig', 'lac fa'), 'LineWidth', lw)
    hold on
    xlabel('Acquisition steps')
    ylabel('Values')
    title('Flip angle values')
    legend
    drawnow

    subplot(1,3,2);
    plot(copy_model_par2.TRList, copy_model2_Mxy(1,:), ...
        'DisplayName', sprintf('%s %s', 'Copy', 'pyr Mxy'), 'LineWidth', lw)
    hold on
    plot(copy_model_par2.TRList, copy_model2_Mxy(2,:), ...
        'DisplayName', sprintf('%s %s', 'Copy', 'lac Mxy'), 'LineWidth', lw)
    hold on
    plot(p_new_model2.TRList, Mxy_new_model2(1,:), ...
        'DisplayName', sprintf('%s %s', 'Orig', 'pyr Mxy'), 'LineWidth', lw)
    hold on
    plot(p_new_model2.TRList, Mxy_new_model2(2,:), ...
        'DisplayName', sprintf('%s %s', 'Orig', 'lac Mxy'), 'LineWidth', lw)
    hold on
    xlabel('Time')
    ylabel('Values')
    title('Mxy')
    legend
    drawnow

    subplot(1,3,3);
    plot(copy_model_par2.TRList, copy_model2_Mz(1,:), ...
        'DisplayName', sprintf('%s %s', 'Copy', 'pyr Mz'), 'LineWidth', lw)
    hold on
    plot(copy_model_par2.TRList, copy_model2_Mz(2,:), ...
        'DisplayName', sprintf('%s %s', 'Copy', 'lac Mz'), 'LineWidth', lw)
    hold on
    plot(p_new_model2.TRList, Mz_new_model2(1,:), ...
        'DisplayName', sprintf('%s %s', 'Orig', 'pyr Mz'), 'LineWidth', lw)
    hold on
    plot(p_new_model2.TRList, Mz_new_model2(2,:), ...
        'DisplayName', sprintf('%s %s', 'Orig', 'lac Mz'), 'LineWidth', lw)
    hold on
    xlabel('Time')
    ylabel('Values')
    title('Mz')
    legend
    drawnow

end




