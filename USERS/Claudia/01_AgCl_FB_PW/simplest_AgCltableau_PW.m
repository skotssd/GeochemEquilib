function [Ag_diss, Cl_diss, x_solid, massErrAg, massErrCl, iterUsed, timeMs, repMethod, repJac] = ...
          simplest_AgCltableau_PW(USE_LOG, USE_ANALYTIC_JAC, pH, pe, TOT_)
% =====================================================================
% Ag–Cl equilibrated with AgCl(s) using PW complementarity (no if).
% Fixed table → SINGLE Newton iteration until convergence.
%
% Inputs:
%    USE_LOG          : bool  (true => log10 variables, false => linear)
%    USE_ANALYTIC_JAC : bool  (true => jac analitic, false => numeric)
%    pH, pe           : scale log H+, e-
%    TOT_             : [AgT; ClT]
%
%  Ouputs:
%    Ag_diss, Cl_diss : final aqueous concentrations (M)
%    x_solid          : moles of solid, suppousing V=1 L (M)
%    massErrAg, massErrCl : mass errors in Ag, Cl
%    iterUsed         : number of Newton iterations
%    timeMs           : total time in ms
%    repMethod        : 'log10' or 'linear'
%    repJac           : 'analitic' or 'numeric'
% =====================================================================

% Determine strings for the final report (outside this function)
if USE_LOG
    repMethod = 'log10';
else
    repMethod = 'lineal';
end
if USE_ANALYTIC_JAC
    repJac = 'analitic';
else
    repJac = 'numeric';
end

% Extract totals
AgT = TOT_(1);
ClT = TOT_(2);

% Define Tableau
Tableau = [ ...
%  H   e-   Ag   Cl    logK      phase
   1   0    0    0      0         0    % H+
   0   1    0    0      0         0    % e-
   0   0    1    0      0         0    % Ag+
   0   0    0    1      0         0    % Cl-
  -1   0    0    0    -14         0    % OH-
   0   0    1    1   9.7453       1];  % AgCl(s)

phaseCol = 6;
aqRows   = (Tableau(:,phaseCol) == 0);
solRows  = (Tableau(:,phaseCol) == 1);

A_sol   = Tableau(aqRows,1:4);
K_sol   = Tableau(aqRows,5);
A_solid = Tableau(solRows,1:4);
K_solid = Tableau(solRows,5);

% Adjust logK according to pH and pe
corrHe    = @(A) A(:,1).*(-pH) + A(:,2).*(-pe);
K_sol_e   = K_sol   + corrHe(A_sol);
K_solid_e = K_solid + corrHe(A_solid);

A_sol_x   = A_sol(:,3:4);
A_solid_x = A_solid(:,3:4);

% Initial guess
if USE_LOG
    X = [log10(AgT*0.8); log10(ClT*0.8); 0];  % [logAg, logCl, xcp]
else
    X = [AgT*0.8; ClT*0.8; 0];               % [Ag, Cl, xcp]
end

tol     = 1e-12;
maxIter = 100;
iter    = 0;

% Evaluate initial residual
[R, SI_val, tmpC] = residualPW(X);

tStart = tic;
while (norm(R) > tol) && (iter < maxIter)
    J   = jacobianPW(X, SI_val, tmpC);
    dX  = -J \ R;
    X   = X + dX;
    [R, SI_val, tmpC] = residualPW(X);
    iter = iter + 1;
end
elapsed = toc(tStart);
timeMs  = elapsed * 1e3;  % miliseconds
iterUsed= iter;

% Calculate concentrations
if USE_LOG
    Caq = 10.^tmpC;
else
    Caq = tmpC;
end

Ag_diss = sum(A_sol_x(:,1).*Caq);
Cl_diss = sum(A_sol_x(:,2).*Caq);
x_solid = X(3);

massErrAg = (Ag_diss + x_solid - AgT);
massErrCl = (Cl_diss + x_solid - ClT);

% ---------------------------------------------------------------------
% Local functions
    function [Rout, SI, Cvec] = residualPW(Xvar)
        ln10 = log(10);

        if USE_LOG
            logAg = Xvar(1);
            logCl = Xvar(2);
            Cvec  = K_sol_e + A_sol_x*[logAg; logCl];  % log10(Ci)
            C     = 10.^Cvec;
        else
            Ag = Xvar(1);
            Cl = Xvar(2);
            C  = zeros(size(K_sol_e));
            for ii = 1:numel(C)
                C(ii) = 10^(K_sol_e(ii)) * Ag^A_sol_x(ii,1)* Cl^A_sol_x(ii,2);
            end
            Cvec = C;
        end

        % Mass Balances
        RmassAg = sum(A_sol_x(:,1).*C) + Xvar(3) - AgT;
        RmassCl = sum(A_sol_x(:,2).*C) + Xvar(3) - ClT;

        % Saturation
        if USE_LOG
            Q  = A_solid_x*[Xvar(1); Xvar(2)];
            SI = 10.^(Q + K_solid_e);
        else
            SI = Xvar(1)^A_solid_x(1,1)* Xvar(2)^A_solid_x(1,2)* 10^K_solid_e;
        end

        % PW complementarity
        u = 1 - SI;
        v = Xvar(3);
        phiPW = 0.5*(u + v - abs(u - v));

        Rout = [RmassAg; RmassCl; phiPW];
    end

    function Jout = jacobianPW(Xvar, SI_, Cvec_)
        if ~USE_ANALYTIC_JAC
            % Numneric Jacobian
            Jout = numJac(@residualPW, Xvar);
            return;
        end

        % Analytic Jacobian
        ln10 = log(10);

        if USE_LOG
            C       = 10.^Cvec_;
            dC_dAg  = ln10*C .* A_sol_x(:,1);
            dC_dCl  = ln10*C .* A_sol_x(:,2);
        else
            Ag      = Xvar(1);
            Cl      = Xvar(2);
            C       = Cvec_;
            dC_dAg  = C .* A_sol_x(:,1)./Ag;
            dC_dCl  = C .* A_sol_x(:,2)./Cl;
        end

        dR_Ag_dAg = sum(A_sol_x(:,1).* dC_dAg);
        dR_Ag_dCl = sum(A_sol_x(:,1).* dC_dCl);
        dR_Cl_dAg = sum(A_sol_x(:,2).* dC_dAg);
        dR_Cl_dCl = sum(A_sol_x(:,2).* dC_dCl);

        u = 1 - SI_;
        v = Xvar(3);
        s = sign(u - v);

        dphi_du = 0.5*(1 - s);
        dphi_dv = 0.5*(1 + s);

        if USE_LOG
            dSI_dAg = ln10 * SI_ * A_solid_x(1,1);
            dSI_dCl = ln10 * SI_ * A_solid_x(1,2);
        else
            Ag = Xvar(1);
            Cl = Xvar(2);
            dSI_dAg = SI_ * A_solid_x(1,1) / Ag;
            dSI_dCl = SI_ * A_solid_x(1,2) / Cl;
        end

        dphi_dAg  = dphi_du * (-dSI_dAg);
        dphi_dCl  = dphi_du * (-dSI_dCl);
        dphi_dxcp = dphi_dv;

        Jout = [ dR_Ag_dAg   dR_Ag_dCl   1
                 dR_Cl_dAg   dR_Cl_dCl   1
                 dphi_dAg    dphi_dCl    dphi_dxcp ];
    end

    function Jnum = numJac(funHandle, xvec)
        hh  = 1e-8;
        fx  = funHandle(xvec);
        nn  = numel(xvec);
        mm  = numel(fx);
        Jnum= zeros(mm, nn);
        for kk = 1:nn
            xh    = xvec;
            xh(kk)= xh(kk) + hh;
            Jnum(:,kk) = (funHandle(xh) - fx)/hh;
        end
    end

end
