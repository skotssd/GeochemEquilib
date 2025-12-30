function [Ag_diss, Cl_diss, x_solid, massErrAg, massErrCl, iterUsed, timeMs, repMethod, repJac] = ...
         complex_AgCltableau_PW(USE_LOG, USE_ANALYTIC_JAC, pH, pe, TOT_)
% =====================================================================
%  Ag–Cl system with aqueous AgCl⁰, AgCl₂⁻, AgCl₃²⁻, AgCl₄³⁻ plus AgCl(s),
%  solved with the P-W complementarity formulation.
%
%  Inputs
%    USE_LOG          – logical, true → solve in log10 variables
%    USE_ANALYTIC_JAC – logical, true → use analytic Jacobian
%    pH, pe           – fixed proton and electron scales
%    TOT_             – [AgT; ClT]  (totals in mol L⁻¹)
%
%  Outputs
%    Ag_diss, Cl_diss – final dissolved Ag and Cl (mol L⁻¹)
%    x_solid          – mol L⁻¹ of AgCl(s) (assuming V = 1 L)
%    massErrAg/Cl     – mass-balance residuals
%    iterUsed         – Newton iterations
%    timeMs           – elapsed time (ms) for this run
%    repMethod        – 'log10'  |  'linear'
%    repJac           – 'analytic'  |  'numeric'
% =====================================================================

% ---------------------------------------------------------------------
% 1.  Reporting strings
% ---------------------------------------------------------------------
repMethod = ternary(USE_LOG,          'log10',   'linear');
repJac    = ternary(USE_ANALYTIC_JAC, 'analytic','numeric');

% ---------------------------------------------------------------------
% 2.  Totals
% ---------------------------------------------------------------------
AgT = TOT_(1);
ClT = TOT_(2);

% ---------------------------------------------------------------------
% 3.  Tableau 
% ---------------------------------------------------------------------
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

TableauNum  = cell2mat(Tableau(:,1:6));
phaseCol    = 6;
aqRows      = TableauNum(:,phaseCol) == 0;
solRows     = TableauNum(:,phaseCol) == 1;

A_sol       = TableauNum(aqRows,1:4);
K_sol       = TableauNum(aqRows,5);
A_solid     = TableauNum(solRows,1:4);
K_solid     = TableauNum(solRows,5);

% ---------------------------------------------------------------------
% 4.  pH/pe correction of log K
% ---------------------------------------------------------------------
corrHe   = @(A) A(:,1).*(-pH) + A(:,2).*(-pe);
K_sol_e  = K_sol   + corrHe(A_sol);
K_solid_e= K_solid + corrHe(A_solid);

A_sol_x   = A_sol(:,3:4);   % [Ag, Cl] in aqueous species
A_solid_x = A_solid(:,3:4); % [Ag, Cl] in solid

% ---------------------------------------------------------------------
% 5.  Initial guess
% ---------------------------------------------------------------------
if USE_LOG
    X = [log10(AgT*0.8);  log10(ClT*0.8);  0];  % [logAg, logCl, xcp]
else
    X = [AgT*0.8;         ClT*0.8;         0];  % [Ag, Cl, xcp]
end

% ---------------------------------------------------------------------
% 6.  Newton iteration
% ---------------------------------------------------------------------
tol = 1e-12;
maxIter = 100;
iterUsed = 0;

[R, SI_val, tmpC] = residualPW(X);
tStart = tic;

while norm(R) > tol && iterUsed < maxIter
    J   = jacobianPW(X, SI_val, tmpC);
    X   = X - J\R;
    [R, SI_val, tmpC] = residualPW(X);
    iterUsed = iterUsed + 1;
end

timeMs = toc(tStart)*1e3;

% ---------------------------------------------------------------------
% 7.  Results
% ---------------------------------------------------------------------
% --- dissolved-species vector
if USE_LOG
    Caq = 10.^tmpC;
else
    Caq = tmpC;
end
Ag_diss = sum(A_sol_x(:,1).*Caq);
Cl_diss = sum(A_sol_x(:,2).*Caq);
x_solid = X(3);

massErrAg = Ag_diss + x_solid - AgT;
massErrCl = Cl_diss + x_solid - ClT;

% ========================= nested functions ==========================

function [Rout, SI, Cvec] = residualPW(Xv)
    ln10 = log(10);
    if USE_LOG
        Cvec = K_sol_e + A_sol_x*Xv(1:2);    % log10 Ci
        C    = 10.^Cvec;
    else
        Ag = Xv(1); Cl = Xv(2);
        C  = 10.^K_sol_e .* (Ag.^A_sol_x(:,1)) .* (Cl.^A_sol_x(:,2));
        Cvec = C;
    end

    RmassAg = sum(A_sol_x(:,1).*C) + Xv(3) - AgT;
    RmassCl = sum(A_sol_x(:,2).*C) + Xv(3) - ClT;

    if USE_LOG
        Q  = A_solid_x*Xv(1:2);
        SI = 10.^(Q + K_solid_e);
    else
        SI = Xv(1)^A_solid_x(1) * Xv(2)^A_solid_x(2) * 10^K_solid_e;
    end

    u = 1 - SI;  v = Xv(3);
    phiPW = 0.5*(u + v - abs(u - v));

    Rout = [RmassAg; RmassCl; phiPW];
end

function J = jacobianPW(Xv, SI, Cvec)
    if ~USE_ANALYTIC_JAC
        J = numJac(@residualPW, Xv);
        return
    end
    ln10 = log(10);
    if USE_LOG
        C      = 10.^Cvec;
        dC_dAg = ln10*C .* A_sol_x(:,1);
        dC_dCl = ln10*C .* A_sol_x(:,2);
    else
        Ag = Xv(1); Cl = Xv(2); C = Cvec;
        dC_dAg = C .* A_sol_x(:,1) / Ag;
        dC_dCl = C .* A_sol_x(:,2) / Cl;
    end

    dR_Ag_dAg = sum(A_sol_x(:,1).*dC_dAg);
    dR_Ag_dCl = sum(A_sol_x(:,1).*dC_dCl);
    dR_Cl_dAg = sum(A_sol_x(:,2).*dC_dAg);
    dR_Cl_dCl = sum(A_sol_x(:,2).*dC_dCl);

    u = 1 - SI; v = Xv(3); s = sign(u - v);
    dphi_du = 0.5*(1 - s); dphi_dv = 0.5*(1 + s);

    if USE_LOG
        dSI_dAg = ln10*SI*A_solid_x(1);
        dSI_dCl = ln10*SI*A_solid_x(2);
    else
        dSI_dAg = SI*A_solid_x(1) / Xv(1);
        dSI_dCl = SI*A_solid_x(2) / Xv(2);
    end

    dphi_dAg  = -dphi_du*dSI_dAg;
    dphi_dCl  = -dphi_du*dSI_dCl;
    dphi_dxcp =  dphi_dv;

    J = [ dR_Ag_dAg  dR_Ag_dCl  1
          dR_Cl_dAg  dR_Cl_dCl  1
          dphi_dAg   dphi_dCl   dphi_dxcp ];
end

function Jn = numJac(fun,x)
    h = 1e-8; fx = fun(x);
    n = numel(x); m = numel(fx); Jn = zeros(m,n);
    for k = 1:n
        xh = x; xh(k) = xh(k)+h;
        Jn(:,k) = (fun(xh)-fx)/h;
    end
end

end  % ── end main function ─────────────────────────────────────────────

% Helper one-liner (ternary operator)
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
