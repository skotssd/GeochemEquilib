function [Ag_diss, Cl_diss, x_solid, massErr] = complex_AgClTableau_FB(pH, pe, TOTALS)

% =========================================================================
%  MAIN FUNCTION: AgCltableau (simplest)
%    Returns [Ag, Cl, AgCls, MASSERR]
% =========================================================================


    % (A) Define the Ag–Cl “Tableau” with H⁺, e⁻, OH⁻ and the solid phase
    %     Columns: [H, e-, Ag, Cl, logK, phase, species]
Tableau = { ...
%  H   e-   Ag   Cl    logK      phase   species
   1   0    0    0     0.0000     0    'H+'        % 1
   0   1    0    0     0.0000     0    'e-'        % 2
   0   0    1    0     0.0000     0    'Ag+'       % 3
   0   0    0    1     0.0000     0    'Cl-'       % 4
  -1   0    0    0   -14.0000     0    'OH-'       % 5
   0   0    1    1     3.2971     0    'AgCl'      % 6
   0   0    1    2     5.2989     0    'AgCl2-'    % 7
   0   0    1    3     5.1310     0    'AgCl3-2'   % 8
   0   0    1    4     3.8050     0    'AgCl4-3'   % 9
   0   0    1    1     9.7453     1    'AgCl(s)'   %10
};

    % (B) Parse the tableau to separate aqueous and solid species
    [A_solution, K_solution, ~, ...
     A_solid,    K_solid,    ~] = parseAgClTableau(Tableau);

    % (C) TOT_ (vector containing AgT y ClT)
    %     (H and e⁻ are fixed by pH and pe; they do not appear in material balances)
    TOT_ = TOTALS;  % TOT_(1)=AgT, TOT_(2)=ClT

    % (D) Solve speciation using the Fischer–Burmeister formulation
    [Ag_diss, Cl_diss, x_solid, massErr] = solveAgCl_FBroutine( ...
        A_solution, K_solution, A_solid, K_solid, pH, pe, TOT_);
endfunction


% -------------------------------------------------------------------------
% parseAgClTableau
%   Extract aqueous rows (phase = 0) and solid rows (phase = 1).
% -------------------------------------------------------------------------
function [A_sol, K_sol, spc_sol, A_solid, K_solid, spc_solid] = parseAgClTableau(Tableau)
    colH     = 1;  % stoich H
    colE     = 2;  % stoich e-
    colAg    = 3;  % stoich Ag
    colCl    = 4;  % stoich Cl
    colLogK  = 5;  % logK
    colPhase = 6;  % phase
    colName  = 7;  % species name (cell array)

    nRows = size(Tableau,1);

    A_sol     = [];  K_sol     = [];  spc_sol   = {};
    A_solid   = [];  K_solid   = [];  spc_solid = {};

    for i = 1:nRows
        hVal     = Tableau(i, colH);
        eVal     = Tableau(i, colE);
        agVal    = Tableau(i, colAg);
        clVal    = Tableau(i, colCl);
        logKVal  = Tableau(i, colLogK);
        phaseVal = Tableau(i, colPhase);
        nameVal  = Tableau{i, colName};  % species string

        % Convert to double if stored in cell arrays
        if iscell(hVal),     hVal     = hVal{1};     end
        if iscell(eVal),     eVal     = eVal{1};     end
        if iscell(agVal),    agVal    = agVal{1};    end
        if iscell(clVal),    clVal    = clVal{1};    end
        if iscell(logKVal),  logKVal  = logKVal{1};  end
        if iscell(phaseVal), phaseVal = phaseVal{1}; end

        stoichVec = [hVal, eVal, agVal, clVal];

        if phaseVal == 0
            % Aqueous species
            A_sol   = [A_sol; stoichVec];
            K_sol   = [K_sol; logKVal];
            spc_sol = [spc_sol; {nameVal}];
        else
            % Solid phase
            A_solid   = [A_solid; stoichVec];
            K_solid   = [K_solid; logKVal];
            spc_solid = [spc_solid; {nameVal}];
        end
    end
endfunction


% -------------------------------------------------------------------------
% solveAgCl_FBroutine
%   Implements Newton–Raphson with the complementarity
%   of Fischer–Burmeister for a single solid.
% -------------------------------------------------------------------------
function [Ag_diss, Cl_diss, x_solid, massErr] = solveAgCl_FBroutine( ...
          A_solution, K_solution, A_solid, K_solid, pH, pe, TOT_)

    % TOT_ = [AgT, ClT].
    % A_solution Nx_sol x 4 => [H, e-, Ag, Cl]
    % A_solid   Nx_solid x 4
    % Adjust K according to pH and pe:
    stoichH_sol = A_solution(:,1);
    stoichE_sol = A_solution(:,2);
    stoichAg_sol= A_solution(:,3);
    stoichCl_sol= A_solution(:,4);

    Kcorr_sol   = K_solution + (stoichH_sol.*(-pH)) + (stoichE_sol.*(-pe));
    A_sol_use   = [stoichAg_sol, stoichCl_sol];

    stoichH_solid = A_solid(:,1);
    stoichE_solid = A_solid(:,2);
    stoichAg_solid= A_solid(:,3);
    stoichCl_solid= A_solid(:,4);

    Kcorr_solid  = K_solid + (stoichH_solid.*(-pH)) + (stoichE_solid.*(-pe));
    A_solid_use  = [stoichAg_solid, stoichCl_solid];

    % Initial guess
    logAg_guess = log10( TOT_(1)*0.8 );
    logCl_guess = log10( TOT_(2)*0.8 );
    xcp_guess   = 0;

    X = [logAg_guess; logCl_guess; xcp_guess];

    tol     = 1e-14;
    maxIter = 100;

    for iter = 1:maxIter
        [R, SI_val, logC] = residualFB_AgCl(X, A_sol_use, Kcorr_sol, TOT_, A_solid_use, Kcorr_solid);
        J = jacobianFB_AgCl(X, A_sol_use, Kcorr_sol, A_solid_use, Kcorr_solid, SI_val, logC);

        dX   = -J \ R;
        Xnew = X + dX;

        % Force xcp >= 0
        if Xnew(3) < 0
            Xnew(3) = 0;
        end

        if norm(Xnew - X) < tol
            X = Xnew;
            break;
        end
        X = Xnew;
    end

    % Final result
    logAg = X(1);
    logCl = X(2);
    xcp   = X(3);

    % Recompute aqueous species
    logC_sol = Kcorr_sol + A_sol_use*[logAg; logCl];
    C_sol    = 10.^logC_sol;

    sumAg = sum( A_sol_use(:,1).* C_sol );
    sumCl = sum( A_sol_use(:,2).* C_sol );

    Ag_diss  = sumAg;  % [Ag(aq)]
    Cl_diss  = sumCl;  % [Cl(aq)]
    x_solid  = xcp;    % [AgCl(s)]
    errAg    = sumAg + xcp - TOT_(1);
    errCl    = sumCl + xcp - TOT_(2);
    massErr  = [errAg, errCl];
endfunction


% -------------------------------------------------------------------------
% residualFB_AgCl
%   R = [Rmass(Ag); Rmass(Cl); RFB], where RFB = FB(xcp, 1 - SI)
% -------------------------------------------------------------------------
function [R, SI, logC] = residualFB_AgCl(X, A_sol, K_sol, TOT_, A_solid, K_solid)
    logAg = X(1);
    logCl = X(2);
    xcp   = X(3);

    logC = K_sol + A_sol*[logAg; logCl];
    C    = 10.^logC;

    sumAg = sum(A_sol(:,1) .* C);
    sumCl = sum(A_sol(:,2) .* C);

    RmassAg = sumAg + xcp - TOT_(1);
    RmassCl = sumCl + xcp - TOT_(2);

    % Saturation Index for solid
    Q  = A_solid*[logAg; logCl];
    SI = 10.^(Q + K_solid);

    % Complementarity (Fischer-Burmeister)
    [valFB, ~, ~] = fischerBurmeister(xcp, 1 - SI);
    RFB = valFB;

    R = [RmassAg; RmassCl; RFB];
endfunction


% -------------------------------------------------------------------------
% jacobianFB_AgCl
%   Jacobian of R using semismooth derivatives
% -------------------------------------------------------------------------
function J = jacobianFB_AgCl(X, A_sol, K_sol, A_solid, K_solid, SI, logC)
    ln10  = log(10);
    logAg = X(1);
    logCl = X(2);
    xcp   = X(3);

    C       = 10.^logC;
    dC_dAg  = ln10 * C .* A_sol(:,1);
    dC_dCl  = ln10 * C .* A_sol(:,2);

    dRmassAg_dAg  = sum( A_sol(:,1) .* dC_dAg );
    dRmassAg_dCl  = sum( A_sol(:,1) .* dC_dCl );
    dRmassAg_dxcp = 1.0;

    dRmassCl_dAg  = sum( A_sol(:,2) .* dC_dAg );
    dRmassCl_dCl  = sum( A_sol(:,2) .* dC_dCl );
    dRmassCl_dxcp = 1.0;

    % FB(xcp, 1 - SI)
    [~, dFB_da, dFB_db] = fischerBurmeister(xcp, 1 - SI);

    % dSI/d(logAg)
    dQ_dAg  = A_solid(1,1);
    dQ_dCl  = A_solid(1,2);
    dSI_dAg = SI * ln10 * dQ_dAg;   % chain rule
    dSI_dCl = SI * ln10 * dQ_dCl;

    % d(1-SI)/d(logAg) = -dSI_dAg
    d1mSI_dAg = - dSI_dAg;
    d1mSI_dCl = - dSI_dCl;

    dRFB_dAg  = dFB_db * d1mSI_dAg;
    dRFB_dCl  = dFB_db * d1mSI_dCl;
    dRFB_dxcp = dFB_da;  % wrt xcp

    J = [
        dRmassAg_dAg,   dRmassAg_dCl,   dRmassAg_dxcp
        dRmassCl_dAg,   dRmassCl_dCl,   dRmassCl_dxcp
        dRFB_dAg,       dRFB_dCl,       dRFB_dxcp
    ];
endfunction


% -------------------------------------------------------------------------
% fischerBurmeister
%   FB(a,b) = sqrt(a^2 + b^2) - a - b
% -------------------------------------------------------------------------
function [val, dFda, dFdb] = fischerBurmeister(a,b)
    r = sqrt(a^2 + b^2);
    val = r - a - b;

    epsVal = 1e-20;
    if r < epsVal
        % fallback -> avoids division by 0
        val  = 0;
        dFda = -1;
        dFdb = -1;
    else
        dFda = (a/r) - 1;
        dFdb = (b/r) - 1;
    end
endfunction