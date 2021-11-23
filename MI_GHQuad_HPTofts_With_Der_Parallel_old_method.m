function [MIobjfun, MIobjfun_Der]=MI_GHQuad_HPTofts_With_Der_Parallel(xopt,M0,params, model, NumberUncertain, xn, wn, xn2, wn2, T1Pqp, T1Lqp, kplqp, klpqp, kveqp, t0qp)

% Update optimization variables
optim_TR = false;
if optim_TR
    TRs = xopt(1:size(params.TRList,2)-1)';
    % compute time sequence from TRs
    N = size(params.TRList,2);
    params.TRList = zeros(1,N);
    for i=2:N
        t = 0;
        for j=1:(i-1)
            t = t + TRs(j);
        end
        params.TRList(i) = t;
    end
    params.FaList = reshape(xopt(size(params.TRList,2):end),size(params.FaList ));
else
    params.FaList = reshape(xopt(1:end),size(params.FaList ));
end

% precompute model solutions at all quad points
lqp=length(xn{1}(:));
disp('computing signal and derivatives')
N = size(params.TRList,2);
Mmodel = zeros(lqp,1);
MmodeldGdTR = zeros(size(params.TRList,1),size(params.TRList,2),lqp);
MmodeldGdFA = zeros(size(params.FaList,1),size(params.FaList,2),lqp);
parfor qp=1:lqp
    params_qp = params;
    params_qp.T1s(1) = T1Pqp(qp);
    params_qp.T1s(2) = T1Lqp(qp);
    params_qp.ExchangeTerms(1,2) = kplqp(qp);
    params_qp.ExchangeTerms(2,1) = klpqp(qp);
    params_qp.PerfusionTerms(1) = kveqp(qp);
    params_qp.t0 = t0qp(qp);
    [~,~,~,G,dGdTR,dGdFA] = model.compile_der(M0.',params_qp);
    Mmodel(qp,1) = G;
    MmodeldGdTR(:,:,qp) = dGdTR(:,:);
    MmodeldGdFA(:,:,qp) = dGdFA(:,:);
end

% Compute quadrature points for nu (omega_k1, ..., omega_kP)
lqp2=length(xn2{1}(:));
signu = params.SignalNoiseMI; 
covnu=signu*signu;

% Compute the ln term as a function of q1,...,qN QPs and k1,...,kP QPs by
% summing over s1,...,sN QPs
Hz = zeros(lqp,1);
dHzdFA = zeros(size(params.FaList,1),size(params.FaList,2), lqp);
dHzdTR = zeros(size(params.TRList,1),size(params.TRList,2), lqp);
parfor iii=1:lqp
    
    Hziii = 0.;
    dHzdFAiii = zeros(size(params.FaList,1),size(params.FaList,2));
    dHzdTRiii = zeros(size(params.TRList,1),size(params.TRList,2));
    
    Giii = Mmodel(iii)';
    
    for jjj=1:lqp2
        znu=xn2{1}(jjj) ;
        lntermtmp=0;
        dFA = zeros(size(params.FaList,1),size(params.FaList,2));
        dTR = zeros(size(params.TRList,1),size(params.TRList,2));
        for kkk=1:lqp
            
            Gkkk = Mmodel(kkk)';
            
            f = wn(kkk) * mvnpdf(znu*sqrt(2)*signu+Giii, Gkkk,covnu);
            lntermtmp=lntermtmp + f;
            % for derivative
            fact_der = (Gkkk - Giii - znu*sqrt(2)*signu)/covnu;
            dTR(:,:) = dTR(:,:) + f * fact_der * (MmodeldGdTR(:,:,iii) - MmodeldGdTR(:,:,kkk));
            dFA(:,:) = dFA(:,:) + f * fact_der * (MmodeldGdFA(:,:,iii) - MmodeldGdFA(:,:,kkk));
        end
        % for function
        lnterm = log(lntermtmp);
        Hziii = Hziii + wn(iii) * wn2(jjj) * lnterm;
        % for derivative
        inv_fact = 0.;
        if abs(lntermtmp) > 0.
            inv_fact = 1/lntermtmp;
        end
        dHzdTRiii = dHzdTRiii + wn(iii) * wn2(jjj) * inv_fact * dTR(:,:);
        dHzdFAiii = dHzdFAiii + wn(iii) * wn2(jjj) * inv_fact * dFA(:,:);
    end
    Hz(iii,1) = -Hziii; 
    dHzdTR(:,:,iii) = -dHzdTRiii;
    dHzdFA(:,:,iii) = -dHzdFAiii;
end

% sum
Hz_red = 0.;
dHzdFA_red = zeros(size(params.FaList,1),size(params.FaList,2));
dHzdTR_red = zeros(size(params.TRList,1),size(params.TRList,2));
for iii=1:lqp
    Hz_red = Hz_red + Hz(iii,1);
    %if iii == 1
    %    size(dHzdFA_red), size(dHzdFA_red(:,:)), size(dHzdFA), size(dHzdFA(:,:,iii))
    %end
    dHzdFA_red(:,:) = dHzdFA_red(:,:) + dHzdFA(:,:,iii);
    dHzdTR_red(:,:) = dHzdTR_red(:,:) + dHzdTR(:,:,iii);
end

%Hz_red, dHzdTR_red, dHzdFA_red

%%
% Compute H(z|eta), which is a function of measurement noise
Hzeta=.5*log((2*pi*2.7183)^5.*signu.^5);

% Compute mutual information as the difference in entropy of the
% evidence, H(z), and entropy of the likelihood, H(z|eta).
MI=Hz_red-Hzeta;

% Make objective function negative, because optimization is minimizing
% objective function and maximizing mutual information.
MIobjfun=-MI;

% return derivative (stack them to get single vector)
% note that TRs params are 1 less as we set first one to be fixed at 0
if optim_TR
    MIobjfun_Der = [-dHzdTR_red(2:end)'; -dHzdFA_red(:)];
else
    MIobjfun_Der = [-dHzdFA_red(:)];
end

end
