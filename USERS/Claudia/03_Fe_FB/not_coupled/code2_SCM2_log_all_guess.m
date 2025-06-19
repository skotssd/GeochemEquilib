function varargout = code2_SCM2_log_all_guess(logH, Caq, xHFO, HFOsi, HFOwi, ASFs, ASFw, varargin)
% code2_SCM2_log_all_guess:
% Solve the adsorption model at two types of sites (Hfos, Hfow)
% using a mininewton in log(Hfos_free), log(Hfow_free).
%
%   Returns 9 outputs if requested:
%      [Caq_out, xSurf, spcNames, Xfinal, err_s, err_w, ...
%         Site_s, Site_w, P_ads]
%
%   Here, the subfunction checkAdsorbedPO4 is used to calculate P_ads.

    iH = 1;  iP = 3;
    H_free = max(Caq(iH),1e-30);
    P_free = max(Caq(iP),1e-30);
    logH_in = log10(H_free);
    logP_in = log10(P_free);

    Site_s_total = max(HFOsi * ASFs * xHFO, 1e-30);
    Site_w_total = max(HFOwi * ASFw * xHFO, 1e-30);

    reacTable = {
       [1 0 1 0],  8.93, 'HfosH'
       [2 0 1 0], 16.22, 'HfosH2+'
       [1 0 0 1],  8.93, 'HfowH'
       [2 0 0 1], 16.22, 'HfowH2+'
       [0 1 0 1], 26.65, 'HfowPO4-4'
       [1 1 0 1], 34.32, 'HfowHPO4-3'
       [2 1 0 1], 40.32, 'HfowH2PO4-2'
       [0 1 1 0], 26.65, 'HfosPO4-4'
       [1 1 1 0], 34.32, 'HfosHPO4-3'
       [2 1 1 0], 40.22, 'HfosH2PO4-2'
    };
    nSp = size(reacTable,1);

    A_scm= zeros(nSp,4);
    K_scm= zeros(nSp,1);
    spcNames= cell(nSp,1);
    for i=1:nSp
        A_scm(i,:)  = reacTable{i,1};
        K_scm(i)    = reacTable{i,2};
        spcNames{i} = reacTable{i,3};
    end

    % indexes "Hfos" y "Hfow"
    idxHfos = false(nSp,1);
    idxHfow = false(nSp,1);
    for i=1:nSp
        if ~isempty(strfind(spcNames{i},'Hfos'))
            idxHfos(i)=true;
        end
        if ~isempty(strfind(spcNames{i},'Hfow'))
            idxHfow(i)=true;
        end
    end

    % guess
    if (nargin>7) && ~isempty(varargin{1})
        X = varargin{1};  % [logHfos, logHfow]
    else
        X = [log10(0.5*Site_s_total); log10(0.5*Site_w_total)];
    end

    tol=1e-12;  
    maxIt=200;  % un poco más alto
    ln10=log(10);

    for it=1:maxIt
        Hfos_free = 10^(X(1));
        Hfow_free = 10^(X(2));

        comp= [logH_in; logP_in; X(1); X(2)];
        logC= A_scm*comp + K_scm;
        C_s= 10.^logC;

        sum_s = sum(C_s(idxHfos));
        sum_w = sum(C_s(idxHfow));

        val_s= max(Hfos_free+sum_s,1e-30);
        val_w= max(Hfow_free+sum_w,1e-30);
        R = [ log(val_s) - log(Site_s_total);
              log(val_w) - log(Site_w_total) ];

        if norm(R,inf)<tol
            break
        end

        dval_s_dHfos = ln10*Hfos_free; 
        dval_s_dHfow = 0;
        dval_w_dHfos = 0;
        dval_w_dHfow = ln10*Hfow_free;

        for k=1:nSp
            cVal= C_s(k);
            a3= A_scm(k,3);  a4= A_scm(k,4);
            if idxHfos(k)
                dval_s_dHfos = dval_s_dHfos + ln10*cVal*a3;
                dval_s_dHfow = dval_s_dHfow + ln10*cVal*a4;
            end
            if idxHfow(k)
                dval_w_dHfos = dval_w_dHfos + ln10*cVal*a3;
                dval_w_dHfow = dval_w_dHfow + ln10*cVal*a4;
            end
        end

        J= [ dval_s_dHfos/val_s, dval_s_dHfow/val_s
             dval_w_dHfos/val_w, dval_w_dHfow/val_w ];

        dX= -J\R;
        alpha=0.1;  % damping
        X= X + alpha*dX;
        X= max(X, -30);
    end

    % salidas
    comp= [logH_in; logP_in; X(1); X(2)];
    logC= A_scm*comp + K_scm;
    C_s= 10.^logC;

    err_s = (10^X(1) + sum(C_s(idxHfos))) - Site_s_total;
    err_w = (10^X(2) + sum(C_s(idxHfow))) - Site_w_total;

    P_ads = checkAdsorbedPO4(C_s, spcNames);

    varargout= {Caq, C_s, spcNames, X, err_s, err_w, Site_s_total, Site_w_total, P_ads};
end

function P_ads = checkAdsorbedPO4(xSurf, spcNames)
    mask = false(size(xSurf));
    for i=1:numel(spcNames)
        if ~isempty(strfind(spcNames{i},'PO4'))
            mask(i)=true;
        end
    end
    P_ads= sum(xSurf(mask));
end