function adaptive_Fe_redox_adsorption_FB_solver9(pHvals, FeT, PT, peINT, HFOsi, HFOwi)
% -------------------------------------------------------------------------
% The script contains TWO versions of the Newton solver.
% A condition on peINT is used to decide which one to use:
% 1. If peINT <= 10 (reducing): It uses the "cautious" solver with a
%    step limiter to handle difficult transitions.
% 2. If peINT > 10 (oxic): It uses the "bold" solver without a step
%    limiter to allow for unrestricted initial precipitation of solids.
% -------------------------------------------------------------------------

%% == (1) Define the 'Tableau' ==
Tableau = {
    % Header (column names):
    'H+'  'e-'  'Fe+3' 'PO4-3' 'SITE_S' 'SITE_W' 'logK'     'phase' 'species'
    % Basis components
     1     0     0      0       0        0          0          0       'H+'
     0     1     0      0       0        0          0          0       'e-'
     0     0     1      0       0        0          0          0       'Fe+3'
     0     0     0      1       0        0          0          0       'PO4-3'
     0     0     0      0       1        0          0          0       'SITE_S'
     0     0     0      0       0        1          0          0       'SITE_W'
    % Kw
    -1     0     0      0       0        0        -14          0       'OH-'
    % Fe(II) / Fe(III) hydrolyzed
     0     1     1      0       0        0         13.92       0       'Fe+2'
    -1     0     1      0       0        0         -2.19       0       'Fe(OH)+2'
    -2     0     1      0       0        0         -5.67       0       'Fe(OH)2+'
    -3     0     1      0       0        0        -12.56       0       'Fe(OH)3'
    -4     0     1      0       0        0        -21.6        0       'Fe(OH)4-'
    -2     0     2      0       0        0         -2.95       0       'Fe2(OH)2+4'
    -4     0     3      0       0        0         -6.3        0       'Fe3(OH)4+5'
    -1     1     1      0       0        0         -3.52       0       'Fe(OH)+'
    -2     1     1      0       0        0         -7.55       0       'Fe(OH)2'
    -3     1     1      0       0        0        -17.98       0       'Fe(OH)3-'
    % Phosphoric acid
     1     0     0      1       0        0         12.346      0       'HPO4-2'
     2     0     0      1       0        0         19.553      0       'H2PO4-'
     3     0     0      1       0        0         21.721      0       'H3PO4'
    % Soluble Fe-P
     1     0     1      1       0        0         17.776      0       'FeHPO4+'
     2     0     1      1       0        0         24.983      0       'FeH2PO4+2'
     1     1     1      1       0        0         28.966      0       'FeHPO4'
     2     1     1      1       0        0         35.273      0       'FeH2PO4+'
    % Surface species
     1     0     0      0       1        0          8.93       0       'HfosH'
     2     0     0      0       1        0         16.22       0       'HfosH2+'
     1     0     0      0       0        1          8.93       0       'HfowH'
     2     0     0      0       0        1         16.22       0       'HfowH2+'
     0     0     0      1       0        1         26.65       0       'HfowPO4-4'
     1     0     0      1       0        1         34.32       0       'HfowHPO4-3'
     2     0     0      1       0        1         40.32       0       'HfowH2PO4-2'
     0     0     0      1       1        0         26.65       0       'HfosPO4-4'
     1     0     0      1       1        0         34.32       0       'HfosHPO4-3'
     2     0     0      1       1        0         40.22       0       'HfosH2PO4-2'
    % Solids
    -3     0     1      0       0        0         -4.891      1       'Fe(OH)3s'
     0     0     1      1       0        0         23.2        1       'FePO4s'
     0     3     3      2       0        0         75.06       1       'Fe3(PO4)2s'
};

%% == (2) Automatically parse the table ==
headerRow = Tableau(1,:);
colLogK_idx    = find(strcmp(headerRow,'logK'));
colPhase_idx   = find(strcmp(headerRow,'phase'));
colSpecies_idx = find(strcmp(headerRow,'species'));
stoichCols_indices = 1:(colLogK_idx-1);
componentNames = headerRow(stoichCols_indices);
dataRows = Tableau(2:end,:);
is_solid_phase = cell2mat(dataRows(:,colPhase_idx)) == 1;
A_aq_surf    = cell2mat(dataRows(~is_solid_phase, stoichCols_indices));
K_aq_surf    = cell2mat(dataRows(~is_solid_phase, colLogK_idx));
spc_aq_surf = dataRows(~is_solid_phase, colSpecies_idx);
idx_H_comp = find(strcmp(componentNames, 'H+'));
idx_e_comp = find(strcmp(componentNames, 'e-'));
idx_Fe3_comp = find(strcmp(componentNames, 'Fe+3'));
idx_PO4_comp = find(strcmp(componentNames, 'PO4-3'));
solid_comp_indices_in_A = [idx_H_comp, idx_e_comp, idx_Fe3_comp, idx_PO4_comp];
A_solid_full   = cell2mat(dataRows(is_solid_phase, stoichCols_indices));
A_solid        = A_solid_full(:, solid_comp_indices_in_A);
K_solid        = cell2mat(dataRows(is_solid_phase, colLogK_idx));
spc_solid      = dataRows(is_solid_phase, colSpecies_idx);

%% == (3) Global parameters and results storage ==
ASFs      = 1;
ASFw      = 1;
tolNR     = 1e-12;
maxIterNR = 500;
nPts      = length(pHvals);

num_aq_surf_species = size(A_aq_surf,1);
num_solid_species = size(A_solid,1);
num_unknowns = 4 + num_solid_species;
Caq_surf_final_store = zeros(num_aq_surf_species, nPts);
xsolid_final_store   = zeros(num_solid_species, nPts);
massErrFe_store      = zeros(1, nPts);
massErrP_store       = zeros(1, nPts);
massErrSiteS_store   = zeros(1, nPts);
massErrSiteW_store   = zeros(1, nPts);
X_solution_store     = zeros(num_unknowns, nPts);
SiteT_s_store        = zeros(1, nPts);
SiteT_w_store        = zeros(1, nPts);

if peINT <= 10
    disp('Reducing conditions detected. Using guess for vivianite.');
    Xguess_global = [-17; -7; -40; -40; 1e-20; 1e-20; FeT/3];
else
    disp('Oxic conditions detected. Using guess for HFO.');
    Xguess_global = [-5; -8; -20; -20; 1e-20; 1e-20; 1e-20];
end
X_prev = Xguess_global;

disp('Starting simulation...');
timeStart = tic;

%% == (4) Main loop over pH ==
for ip = 1:nPts
    pH = pHvals(ip);
    pe = peINT - pH;
    logH = -pH;
    logE = -pe;

    if ip == 1
        X = Xguess_global;
    else
        X = X_prev;
    end
    
    successFlag = false;
    
    if peINT <= 10
        for iter_nr = 1:maxIterNR
            [R, ~, ~, ~, ~] = residualFun_fully_coupled(X, A_aq_surf, K_aq_surf, logH, logE, FeT, PT, HFOsi, HFOwi, ASFs, ASFw, A_solid, K_solid, componentNames);
            if norm(R, inf) < tolNR, successFlag = true; break; end
            J = jacobianFun_fully_coupled(X, A_aq_surf, K_aq_surf, logH, logE, HFOsi, HFOwi, ASFs, ASFw, A_solid, K_solid, componentNames, X);
            dX = -J \ R;

            max_log_step = 2.0;
            max_solid_step_frac = 0.5;
            dX(1:4) = max(min(dX(1:4), max_log_step), -max_log_step);
            max_solid_change = max(FeT, PT) * max_solid_step_frac;
            dX(5:7) = max(min(dX(5:7), max_solid_change), -max_solid_change);

            alpha = 1.0;
            norm_R = norm(R,inf);
            while alpha > 1e-8
                Xnew = X + alpha * dX;
                Xnew(5:end) = max(Xnew(5:end), 1e-30);
                Xnew(1:4) = min(max(Xnew(1:4), -40), 10);
                R_new = residualFun_fully_coupled(Xnew, A_aq_surf, K_aq_surf, logH, logE, FeT, PT, HFOsi, HFOwi, ASFs, ASFw, A_solid, K_solid, componentNames);
                if any(~isfinite(R_new)), alpha = alpha / 2; continue; end
                if norm(R_new, inf) < norm_R, X = Xnew; break; else, alpha = alpha / 2; end
            end
            if alpha <= 1e-8, X = Xnew; end
        end
    else
        for iter_nr = 1:maxIterNR
            [R, ~, ~, ~, ~] = residualFun_fully_coupled(X, A_aq_surf, K_aq_surf, logH, logE, FeT, PT, HFOsi, HFOwi, ASFs, ASFw, A_solid, K_solid, componentNames);
            if norm(R, inf) < tolNR, successFlag = true; break; end
            J = jacobianFun_fully_coupled(X, A_aq_surf, K_aq_surf, logH, logE, HFOsi, HFOwi, ASFs, ASFw, A_solid, K_solid, componentNames, X);
            dX = -J \ R;
            
            alpha = 1.0;
            norm_R = norm(R,inf);
            while alpha > 1e-8
                Xnew = X + alpha * dX;
                Xnew(5:end) = max(Xnew(5:end), 1e-30);
                Xnew(1:4) = min(max(Xnew(1:4), -40), 10);
                R_new = residualFun_fully_coupled(Xnew, A_aq_surf, K_aq_surf, logH, logE, FeT, PT, HFOsi, HFOwi, ASFs, ASFw, A_solid, K_solid, componentNames);
                if any(~isfinite(R_new)), alpha = alpha / 2; continue; end
                if norm(R_new, inf) < (1 - alpha*0.0001) * norm_R, X = Xnew; break; else, alpha = alpha / 2; end
            end
            if alpha <= 1e-8, X = Xnew; end
        end
    end
    
    if ~successFlag
        fprintf('Newton-Raphson failed to converge at pH = %.2f after %d iterations.\n', pH, maxIterNR);
        break; 
    else
        [~, ~, C_final, SiteT_s_final, SiteT_w_final] = residualFun_fully_coupled(X, A_aq_surf, K_aq_surf, logH, logE, FeT, PT, HFOsi, HFOwi, ASFs, ASFw, A_solid, K_solid, componentNames);
        Caq_surf_final_store(:,ip) = C_final;
        xsolid_final_store(:,ip)   = X(5:end);
        X_solution_store(:,ip)     = X;
        SiteT_s_store(ip)          = SiteT_s_final;
        SiteT_w_store(ip)          = SiteT_w_final;
        
        idx_Fe_comp_in_A = find(strcmp(componentNames, 'Fe+3'));
        idx_PO4_comp_in_A= find(strcmp(componentNames, 'PO4-3'));
        idx_SIT_S_in_A = find(strcmp(componentNames, 'SITE_S'));
        idx_SIT_W_in_A = find(strcmp(componentNames, 'SITE_W'));
        sumFe_calc    = sum(A_aq_surf(:,idx_Fe_comp_in_A) .* C_final);
        sumP_calc     = sum(A_aq_surf(:,idx_PO4_comp_in_A) .* C_final);
        sumSiteS_calc = sum(A_aq_surf(:,idx_SIT_S_in_A) .* C_final);
        sumSiteW_calc = sum(A_aq_surf(:,idx_SIT_W_in_A) .* C_final);
        Fe_solid_calc = X(5) * A_solid_full(1,idx_Fe_comp_in_A) + X(6) * A_solid_full(2,idx_Fe_comp_in_A) + X(7) * A_solid_full(3,idx_Fe_comp_in_A);
        P_solid_calc  = X(5) * A_solid_full(1,idx_PO4_comp_in_A) + X(6) * A_solid_full(2,idx_PO4_comp_in_A) + X(7) * A_solid_full(3,idx_PO4_comp_in_A);
        massErrFe_store(ip)    = sumFe_calc + Fe_solid_calc - FeT;
        massErrP_store(ip)     = sumP_calc  + P_solid_calc  - PT;
        massErrSiteS_store(ip) = sumSiteS_calc - SiteT_s_final;
        massErrSiteW_store(ip) = sumSiteW_calc - SiteT_w_final;
    end
    X_prev = X;
end

elapsedTime = toc(timeStart);
fprintf('Simulation completed.\n');
fprintf('Total execution time: %.2f seconds\n', elapsedTime);

%% == (5) Plot and Save (Safe Version) ==
fprintf('Generating figures...\n');

% --- Preparing data for figures ---
idx_Fe3 = find(strcmp(spc_aq_surf,'Fe+3'));
Fe3_aq_conc = Caq_surf_final_store(idx_Fe3,:);
HFO_P_ads = zeros(1,nPts);
idx_P = find(strcmp(componentNames, 'PO4-3'));
idx_Ss = find(strcmp(componentNames, 'SITE_S')); 
idx_Sw = find(strcmp(componentNames, 'SITE_W'));

for ip=1:nPts
    if ~isnan(Caq_surf_final_store(1,ip))
        for k=1:num_aq_surf_species
            % Sum if it is a surface species containing P
            if (A_aq_surf(k,idx_Ss)>0 || A_aq_surf(k,idx_Sw)>0) && A_aq_surf(k,idx_P)>0
                HFO_P_ads(ip) = HFO_P_ads(ip) + A_aq_surf(k,idx_P)*Caq_surf_final_store(k,ip);
            end
        end
    end
end

% --- Internal Function to Save (Standard PNG Driver) ---
function save_safe(hFig, name)
    set(hFig, 'PaperUnits', 'inches');
    set(hFig, 'PaperPosition', [0 0 10 6]); % Size 10x6 inches
    set(hFig, 'PaperSize', [10 6]);
    
    % KEY: Use "-dpng" (without cairo) to ensure compatibility
    print(hFig, name, "-dpng", "-r120");
    close(hFig); % Close to prevent memory saturation
end

% --- Generate FIGURE 1 ---
h1 = figure(1, 'visible', 'off');
plot(pHvals, xsolid_final_store(1,:), '-or', 'DisplayName', 'Fe(OH)3s'); hold on;
plot(pHvals, xsolid_final_store(2,:), '-sb', 'DisplayName', 'FePO4(s)');
plot(pHvals, xsolid_final_store(3,:), '-^g', 'DisplayName', 'Fe3(PO4)2(s)');
plot(pHvals, Fe3_aq_conc, '--m', 'DisplayName', 'Fe^{3+} (aq)');
plot(pHvals, HFO_P_ads, '-xk', 'DisplayName', 'HFO-P (ads)');
title('Species Distribution'); xlabel('pH'); ylabel('Concentration (M)');
legend("Location", "best"); grid on;
save_safe(h1, "Fig1.png");

% --- Generate FIGURE 2 ---
h2 = figure(2, 'visible', 'off');
semilogy(pHvals, abs(massErrFe_store), '-r', 'DisplayName', 'Err Fe'); hold on;
semilogy(pHvals, abs(massErrP_store), '-b', 'DisplayName', 'Err P');
title('Mass Balance Error'); xlabel('pH'); ylabel('Error (M)');
legend("Location", "best"); grid on;
save_safe(h2, "Fig2.png");

% --- Generate FIGURE 3 ---
h3 = figure(3, 'visible', 'off');
subplot(2,1,1); plot(pHvals, massErrSiteS_store); hold on; plot(pHvals, massErrSiteW_store);
title('Surface Site Balance Error'); grid on;
subplot(2,1,2); plot(pHvals, SiteT_s_store); hold on; plot(pHvals, SiteT_w_store);
title('Total Sites'); grid on;
save_safe(h3, "Fig3.png");

% --- Generate FIGURE 4 ---
all_err = [abs(massErrFe_store); abs(massErrP_store); abs(massErrSiteS_store); abs(massErrSiteW_store)];
h4 = figure(4, 'visible', 'off');
semilogy(pHvals, max(all_err,[],1), '-k', 'LineWidth', 2);
title('Maximum System Error'); xlabel('pH'); grid on;
save_safe(h4, "Fig4.png");

fprintf('Done! Images saved: Fig1.png, Fig2.png, Fig3.png, Fig4.png\n');
fprintf('Use the Markdown cell below to view them.\n');

end

%% == (6) Helper Functions ==
function [R, logC_aq_surf, C_aq_surf, SiteT_s, SiteT_w] = residualFun_fully_coupled(X, A_aq_surf, K_aq_surf, logH, logE, FeT, PT, HFOsi, HFOwi, ASFs, ASFw, A_solid, K_solid, componentNames)
    logFe3 = X(1); logPO4 = X(2); logSITE_S_free = X(3); logSITE_W_free = X(4);
    xFeOH3 = X(5); xFePO4 = X(6); xFe3PO42 = X(7);
    compVec_aq_surf = [logH; logE; logFe3; logPO4; logSITE_S_free; logSITE_W_free];
    logC_aq_surf = K_aq_surf + A_aq_surf * compVec_aq_surf;
    C_aq_surf = 10.^logC_aq_surf;
    C_aq_surf(~isfinite(C_aq_surf)) = 0;
    idx_Fe_comp_in_A = find(strcmp(componentNames, 'Fe+3'));
    idx_PO4_comp_in_A = find(strcmp(componentNames, 'PO4-3'));
    idx_SIT_S_comp_in_A = find(strcmp(componentNames, 'SITE_S'));
    idx_SIT_W_comp_in_A = find(strcmp(componentNames, 'SITE_W'));
    sumFe_aq_surf = sum(A_aq_surf(:,idx_Fe_comp_in_A) .* C_aq_surf);
    sumP_aq_surf = sum(A_aq_surf(:,idx_PO4_comp_in_A) .* C_aq_surf);
    sumSITE_S_aq_surf = sum(A_aq_surf(:,idx_SIT_S_comp_in_A).* C_aq_surf);
    sumSITE_W_aq_surf = sum(A_aq_surf(:,idx_SIT_W_comp_in_A).* C_aq_surf);
    Fe_solid = xFeOH3 * 1 + xFePO4 * 1 + xFe3PO42 * 3;
    P_solid = xFeOH3 * 0 + xFePO4 * 1 + xFe3PO42 * 2;
    rFe = sumFe_aq_surf + Fe_solid - FeT;
    rP = sumP_aq_surf + P_solid - PT;
    SiteT_s = HFOsi * ASFs * xFeOH3;
    SiteT_w = HFOwi * ASFw * xFeOH3;
    rSiteS = sumSITE_S_aq_surf - SiteT_s;
    rSiteW = sumSITE_W_aq_surf - SiteT_w;
    Rfb = zeros(3,1);
    compVec_solids = [logH; logE; logFe3; logPO4];
    for iSol = 1:size(A_solid,1)
        log_SI_val = A_solid(iSol,:) * compVec_solids + K_solid(iSol);
        SI_i = 10^log_SI_val;
        x_i = X(4 + iSol);
        valFB = fischerBurmeister(x_i, 1 - SI_i);
        Rfb(iSol) = valFB;
    end
    R = [rFe; rP; rSiteS; rSiteW; Rfb];
end

function J = jacobianFun_fully_coupled(X, A_aq_surf, K_aq_surf, logH, logE, HFOsi, HFOwi, ASFs, ASFw, A_solid, K_solid, componentNames, X_full)
    logFe3 = X_full(1); logPO4 = X_full(2); logSITE_S_free = X_full(3); logSITE_W_free = X_full(4);
    compVec_aq_surf = [logH; logE; logFe3; logPO4; logSITE_S_free; logSITE_W_free];
    logC_aq_surf = A_aq_surf * compVec_aq_surf + K_aq_surf;
    num_unknowns = length(X_full);
    J = zeros(num_unknowns, num_unknowns);
    ln10 = log(10);
    C_aq_surf = 10.^logC_aq_surf;
    C_aq_surf(~isfinite(C_aq_surf)) = 0;
    idx_Fe_comp_in_A = find(strcmp(componentNames, 'Fe+3'));
    idx_PO4_comp_in_A = find(strcmp(componentNames, 'PO4-3'));
    idx_SIT_S_comp_in_A = find(strcmp(componentNames, 'SITE_S'));
    idx_SIT_W_comp_in_A = find(strcmp(componentNames, 'SITE_W'));
    map_X_to_A_col = [idx_Fe_comp_in_A, idx_PO4_comp_in_A, idx_SIT_S_comp_in_A, idx_SIT_W_comp_in_A];
    for k_X = 1:4
        dC_dlogXk = C_aq_surf .* ln10 .* A_aq_surf(:, map_X_to_A_col(k_X));
        J(1, k_X) = sum(A_aq_surf(:,idx_Fe_comp_in_A) .* dC_dlogXk);
        J(2, k_X) = sum(A_aq_surf(:,idx_PO4_comp_in_A) .* dC_dlogXk);
        J(3, k_X) = sum(A_aq_surf(:,idx_SIT_S_comp_in_A) .* dC_dlogXk);
        J(4, k_X) = sum(A_aq_surf(:,idx_SIT_W_comp_in_A) .* dC_dlogXk);
    end
    J(1, 5) = 1; J(1, 6) = 1; J(1, 7) = 3;
    J(2, 5) = 0; J(2, 6) = 1; J(2, 7) = 2;
    J(3, 5) = -(HFOsi * ASFs); J(3, 6) = 0; J(3, 7) = 0;
    J(4, 5) = -(HFOwi * ASFw); J(4, 6) = 0; J(4, 7) = 0;
    for iSol = 1:size(A_solid,1)
        current_row_in_J = 4 + iSol;
        x_i_solid = X_full(4 + iSol);
        log_SI_val = A_solid(iSol,:) * [logH; logE; logFe3; logPO4] + K_solid(iSol);
        SI_i = 10^log_SI_val;
        [~, dFB_dxSolid, dFB_dOneMinusSI] = fischerBurmeister(x_i_solid, 1 - SI_i);
        dSI_dlogFe3 = SI_i * ln10 * A_solid(iSol, 3);
        J(current_row_in_J, 1) = dFB_dOneMinusSI * (-dSI_dlogFe3);
        dSI_dlogPO4 = SI_i * ln10 * A_solid(iSol, 4);
        J(current_row_in_J, 2) = dFB_dOneMinusSI * (-dSI_dlogPO4);
        J(current_row_in_J, 3) = 0;
        J(current_row_in_J, 4) = 0;
        J(current_row_in_J, 4 + iSol) = dFB_dxSolid;
    end
end

function [val, dval_da, dval_db] = fischerBurmeister(a,b)
    r = sqrt(a^2 + b^2);
    val = r - a - b;
    epsVal = 1e-20;
    if nargout > 1
        if r < epsVal
            dval_da = a/(r+epsVal) - 1;
            dval_db = b/(r+epsVal) - 1;
        else
            dval_da = (a/r) - 1;
            dval_db = (b/r) - 1;
        end
    end
end