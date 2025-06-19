function R = code1_FB_main_debugged2(pHvals, FeT, PT, NaT, ClT, ...
                                  HFOsi, HFOwi, ASFs, ASFw,peIN)
% code1_FB_main_debug2:
%   Version 
%   without damping and with maxIterSpec=100 
%
%   R = code1_FB_main_debug2(pHvals, FeT, PT, NaT, ClT, HFOsi, HFOwi, ASFs, ASFw)
%
% Inputs:
%   pHvals, FeT, PT, NaT, ClT, HFOsi, HFOwi, ASFs, ASFw
%
% Outputs:
%   R: struct array con un elemento por cada pH, conteniendo:
%       .pH
%       .Caq
%       .xHFO, .xFePO4, .xFe3PO4
%       .P_ads
%       .errFe, .errP
%
% Use in jupyter:
%   results = code1_FB_main_debug2(pHvals, FeT, PT, NaT, ClT, HFOsi, HFOwi, ASFs, ASFw);

    clc;

    % == (1) The 'tableau' (aqueous phase + precip) ==
    Tableau = {
      'H+'   'e-'   'Fe+3'  'PO4-3'  'logK'  'phase' 'species'
       1      0      0       0        0       0      'H+'
       0      1      0       0        0       0      'e-'
       0      0      1       0        0       0      'Fe+3'
       0      0      0       1        0       0      'PO4-3'
       0      0      0       0        1       0      'Cl-'
       0      0      0       0        0       0      'Na+'
      -1      0      0       0      -14       0      'OH-'
       0      1      1       0      13.92     0      'Fe+2'
      -1      0      1       0      -2.19     0      'Fe(OH)+2'
      -2      0      1       0      -5.67     0      'Fe(OH)2+'
      -3      0      1       0     -12.56     0      'Fe(OH)3'
      -4      0      1       0     -21.6      0      'Fe(OH)4-'
      -3      0      1       0      -4.891    1      'Fe(OH)3s'
       0      0      1       1      23.2      1      'Fe(PO4)s'
       0      3      3       2      75.06     1      'Fe3(PO4)2s'
       1      0      0       1     12.346     0      'HPO4-2'
       2      0      0       1     19.553     0      'H2PO4-'
       3      0      0       1     21.721     0      'H3PO4'
       1      0      1       1     17.776     0      'FeHPO4+'
       2      0      1       1     24.983     0      'FeH2PO4+2'
    };

    % == (2) Parse tableau ==
    [tSto, kCol, phCol, spCol] = identifyColumns(Tableau(1,:));
    [Aaq, Kaq, spAq, As, Ks, spS] = parseTableauAuto(Tableau, tSto, kCol, phCol, spCol);

    % == (3) Global parameters ==
    %  Check if NaT, ClT, etc., were passed, they are not used
    tol_global    = 1e-20;
    maxIterGlobal = 100;
    maxIterSpec   = 100;  % sin damping

    nPts = numel(pHvals);

    % Prealloc out R
    R(nPts) = struct();

    % Spaces to store errors, etc.
    massErrFe_store = zeros(1,nPts);
    massErrP_store  = zeros(1,nPts);

    % Warm-start
    Xguess_global = [-3; -5; 0; 0; 0];
    Xsol_prev     = Xguess_global;

    tic();
    for ip = 1:nPts
        pH = pHvals(ip);
        pe = peIN - pH;
        logH = -pH;
        logE = -pe;

        if ip==1
            Xguess = Xguess_global;
        else
            Xguess = Xsol_prev;
        end

        PT_free = PT;

        for itOut=1:maxIterGlobal
            % Aqueous phase + precip: 
            [Xsol, Caq_sol, ok] = subsolve_aqueous_noDamp(Aaq, Kaq, As, Ks, ...
                FeT, PT_free, logH, logE, Xguess, maxIterSpec);

            if ~ok
                % If it fails, try again with a "softer" guess
                Xguess = [-3; -5; 0; 0; 0];
                [Xsol, Caq_sol, ok] = subsolve_aqueous_noDamp(Aaq, Kaq, As, Ks, ...
                    FeT, PT_free, logH, logE, Xguess, maxIterSpec);
            end

            xHFO    = Xsol(3);
            xFePO4  = Xsol(4);
            xFe3PO4 = Xsol(5);

            % Adsorption
            if xHFO < 1e-12
                P_ads = 0;
            else
                [Caq_ads, xSurf, spcNamesSurf, Xfinal, err_s, err_w, site_s, site_w] = ...
                   code2_SCM2_log_all_guess(logH, Caq_sol, xHFO, HFOsi, HFOwi, ASFs, ASFw);

                P_ads = checkAdsorbedPO4(xSurf, spcNamesSurf);

                % Actualizamos PTfree
                new_PT_free = PT - P_ads;
                if new_PT_free < 1e-30
                    new_PT_free = 1e-30;
                end

                % guess: to adjust logPO4 according to Caq_ads
                iPO4 = find(strcmp(spAq,'PO4-3'));
                if isempty(iPO4), iPO4=4; end
                guess_logPO4 = log10( max(Caq_ads(iPO4),1e-30) );
                Xguess = [ Xsol(1); guess_logPO4; xHFO; xFePO4; xFe3PO4 ];
            end

            if xHFO<1e-12
                new_PT_free = PT;  % If there is no HFO it does not adsorb
            end

            if abs(new_PT_free - PT_free) < tol_global
                break;
            end
            PT_free = new_PT_free;
        end

        % At the end of the outer convergence, we save
        [sumFe, sumP] = sumFeP_Aq(Caq_sol, Aaq);
        Fe_solid = xHFO + xFePO4 + 3*xFe3PO4;
        P_solid  =          xFePO4 + 2*xFe3PO4;

        massErrFe_store(ip)= abs(sumFe + Fe_solid - FeT);
        massErrP_store(ip) = abs(sumP + P_solid - PT);

        R(ip).pH     = pH;
        R(ip).Caq    = Caq_sol;
        R(ip).xHFO   = xHFO;
        R(ip).xFePO4 = xFePO4;
        R(ip).xFe3PO4= xFe3PO4;
        R(ip).P_ads  = (xHFO>1e-12)*P_ads;  % 0 si HFO<1e-12
        R(ip).errFe  = massErrFe_store(ip);
        R(ip).errP   = massErrP_store(ip);

        Xsol_prev = Xsol;  % warm-start
    end
    elapsedT = toc();
    fprintf('End code1_FB_main_debug2. t=%.3f s\n', elapsedT);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunc: subsolve_aqueous_noDamp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xsol, Caq_sol, successFlag] = subsolve_aqueous_noDamp( ...
    A_sol, K_sol, A_solid, K_solid, ...
    FeT, PTfree, logH, logE, Xguess, maxIterSpec)

    % without damping. 
    X = Xguess;
    tolNR=1e-14;
    successFlag=false;

    for iter=1:maxIterSpec
        [R, logC] = residualFunFeP_3solids(X, A_sol, K_sol, ...
                                           logH, logE, FeT, PTfree, ...
                                           A_solid, K_solid);
        J = jacobianFunFeP_3solids(X, A_sol, logC, logH, logE, ...
                                   A_solid, K_solid);

        dX  = -J \ R;   % no damping
        Xnew= X + dX;

        % force xSolid>=0
        for k=3:5
            if Xnew(k)<0
                Xnew(k)=0;
            end
        end

        if norm(Xnew - X,inf)<tolNR
            X= Xnew;
            successFlag=true;
            break
        end
        X= Xnew;
    end

    Xsol= X;
    % Solve Caq
    logCfinal= K_sol + A_sol*[logH; logE; Xsol(1); Xsol(2)];
    Caq_sol  = 10.^logCfinal;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunc: residualFunFeP_3solids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R, logC] = residualFunFeP_3solids(X, A_sol, K_sol, ...
                                            logH, logE, FeT, PTfree, ...
                                            A_solid, K_solid)
    % X= [logFe3, logPO4, xHFO, xFePO4, xFe3PO4]
    logFe3= X(1);
    logPO4= X(2);
    xHFO   = X(3);
    xFePO4 = X(4);
    xFe3PO4= X(5);

    compVec= [logH; logE; logFe3; logPO4];
    logC   = K_sol + A_sol*compVec;
    Caq    = 10.^logC;

    [sumFe, sumP] = sumFeP_Aq(Caq, A_sol);

    Fe_solid= xHFO + xFePO4 + 3*xFe3PO4;
    P_solid=          xFePO4 + 2*xFe3PO4;

    rFe= sumFe + Fe_solid - FeT;
    rP = sumP + P_solid - PTfree;

    % Fischer-Burmeister (3 solids)
    Rfb= zeros(3,1);
    for iSol=1:3
       Q_i= A_solid(iSol,:)*compVec + K_solid(iSol);
       SI_i= 10^Q_i;
       x_i= X(2 + iSol);

       [valFB,~,~] = fischerBurmeister(x_i, 1 - SI_i);
       Rfb(iSol)= valFB;
    end

    R= [rFe; rP; Rfb];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunc: jacobianFunFeP_3solids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = jacobianFunFeP_3solids(X, A_sol, logC, logH, logE, ...
                                    A_solid, K_solid)
    ln10 = log(10);

    logFe3= X(1);
    logPO4= X(2);
    xHFO  = X(3);
    xFePO4= X(4);
    xFe3PO4= X(5);

    Caq= 10.^logC;
    stoFe= A_sol(:,3);
    stoP = A_sol(:,4);

    dsumFe_dFe3= 0; dsumFe_dPO4= 0;
    dsumP_dFe3=  0; dsumP_dPO4=  0;

    for i=1:length(Caq)
        cVal= Caq(i);
        dC_dFe3= cVal*ln10*A_sol(i,3);
        dC_dPO4= cVal*ln10*A_sol(i,4);

        dsumFe_dFe3= dsumFe_dFe3 + stoFe(i)*dC_dFe3;
        dsumFe_dPO4= dsumFe_dPO4 + stoFe(i)*dC_dPO4;

        dsumP_dFe3=  dsumP_dFe3 + stoP(i)*dC_dFe3;
        dsumP_dPO4=  dsumP_dPO4 + stoP(i)*dC_dPO4;
    end

    drFe= [dsumFe_dFe3, dsumFe_dPO4, 1, 1, 3];
    drP = [dsumP_dFe3,  dsumP_dPO4, 0, 1, 2];

    compVec= [logH; logE; logFe3; logPO4];
    Jfb= zeros(3,5);

    for iSol=1:3
       Q_i= A_solid(iSol,:)*compVec + K_solid(iSol);
       SI_i= 10^Q_i;
       x_i= X(2 + iSol);

       [~, dFda, dFdb] = fischerBurmeister(x_i, 1 - SI_i);

       dSI_dFe3= SI_i*ln10*A_solid(iSol,3);
       dSI_dPO4= SI_i*ln10*A_solid(iSol,4);

       dFB= zeros(1,5);
       dFB(1)= dFdb * (- dSI_dFe3);
       dFB(2)= dFdb * (- dSI_dPO4);
       dFB(2 + iSol)= dFda;  

       Jfb(iSol,:)= dFB;
    end

    J= [drFe;
        drP;
        Jfb(1,:);
        Jfb(2,:);
        Jfb(3,:)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunc: sumFeP_Aq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sumFe, sumP] = sumFeP_Aq(Caq, A_sol)
    stoFe= A_sol(:,3);
    stoP = A_sol(:,4);
    sumFe= sum(stoFe.*Caq);
    sumP = sum(stoP .*Caq);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunc: fischerBurmeister
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [val, dval_da, dval_db] = fischerBurmeister(a,b)
    r = sqrt(a^2 + b^2);
    val= r - a - b;
    epsVal=1e-30;
    if r<epsVal
        val=0;
        dval_da= -1;
        dval_db= -1;
    else
        dval_da= (a/r)-1;
        dval_db= (b/r)-1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunc: checkAdsorbedPO4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pAds = checkAdsorbedPO4(xSurf, spcNames)
    pAds=0;
    for i=1:numel(xSurf)
        if ~isempty(strfind(spcNames{i}, 'PO4'))
            pAds = pAds + xSurf(i);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunc: parseTableauAuto
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A_sol, K_sol, spc_sol, A_solid, K_solid, spc_solid] = ...
    parseTableauAuto(Tableau, titleStoichCols, colLogK, colPhase, colSpec)
    data = Tableau(2:end,:);
    nStoich= numel(titleStoichCols);

    A_sol=[];  K_sol=[];  spc_sol={};
    A_solid=[];K_solid=[];spc_solid={};

    for i=1:size(data,1)
        rowStoich= zeros(1,nStoich);
        for s=1:nStoich
            rowStoich(s)= data{i,s};
        end
        valK= data{i,colLogK};
        ph = data{i,colPhase};
        sn = data{i,colSpec};

        if ph==0
            A_sol= [A_sol; rowStoich];
            K_sol= [K_sol; valK];
            spc_sol= [spc_sol; sn];
        else
            A_solid= [A_solid; rowStoich];
            K_solid= [K_solid; valK];
            spc_solid= [spc_solid; sn];
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunc: identifyColumns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [titleStoichCols, colLogK, colPhase, colSpecies] = identifyColumns(headerRow)
    colLogK   = find(strcmp(headerRow,'logK'));
    colPhase  = find(strcmp(headerRow,'phase'));
    colSpecies= find(strcmp(headerRow,'species'));
    colStoichEnd= colLogK - 1;
    titleStoichCols= headerRow(1:colStoichEnd);
end
