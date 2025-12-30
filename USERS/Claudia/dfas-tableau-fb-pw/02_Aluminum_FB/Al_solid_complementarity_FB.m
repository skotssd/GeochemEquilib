function [results, results2D] = Al_solid_complementarity_FB(AlT, pHvals)
% al_solid_complementarity_FB
%
% Aluminum(III) speciation with a single Al(OH)3(s) solid phase,
% using the Fischer-Burmeister formulation (ONE solid phase).
% We ignore e- and Cl- columns, so A_solid and K_solid are scalars.
%
% Input:
%   - AlT    : total aluminum concentration (M)
%   - pHvals : vector of pH values for 1D sweep
%
% Output:
%   - results  : structure containing 1D sweep data (pH, species, SI, etc.)
%   - results2D: structure containing 2D “surface” data (pH vs. AlT)
%
% (Code by Claudia Santiviago Petzoldt)
% -------------------------------------------------------------------------

% (1) Define the Tableau with header row
Tableau = {...
  {'H+'},   {'e-'},   {'Al+3'},  {'Cl-'},   {'logK'},  {'phase'},   {'species'};
   1        0         0         0          0         0           {'H+'}
   0        1         0         0          0         0           {'e-'}
   0        0         1         0          0         0           {'Al^{3+}'}
   0        0         0         1          0         0           {'Cl-'}
  -1        0         0         0       -13.9951     0           {'OH-'}
  -2        0         1         0       -10.5945     0           {'Al(OH)2+'}
  -1        0         1         0       -4.9571      0           {'Al(OH)^{2+}'}
  -4        0         1         0       -22.7        0           {'Al(OH)4-'}
  -2        0         2         0       -7.7902      0           {'Al2(OH)2^{4+}'}
  -4        0         3         0       -13.8803     0           {'Al3(OH)4^{5+}'}
  -32       0         13        0       -98.73       0           {'Al13(OH)32^{7+}'}
  % SINGLE solid (phase=1)
  -3        0         1         0       -7.7560      1           {'Al(OH)3(s)'}
};

% (2) Parse the Tableau
[A_solution, K_solution, species_solution, A_solid, K_solid, species_solid] = ...
    parseAlTableau(Tableau);

% (3) Parameters for 1D pH sweep
tol     = 1e-12;   % Newton tolerance
maxIter = 100;     % max. Newton iterations

% (4) Arrays 1D
nPts    = length(pHvals);
nSol    = size(A_solution,1);
C_all   = zeros(nSol, nPts);
xcp_vec = zeros(1, nPts);
SI_vec  = zeros(1, nPts);
massErr = zeros(1, nPts);

% Initial guess
Xguess = [-3; 0];  % [logAl3 ; xcp]

% (5) 1D Loop over pH
for ip = 1:nPts
    pH   = pHvals(ip);
    logH = -pH;

    X = Xguess;
    for iter = 1:maxIter
        [R, SI_val, logC] = residualFB_Al(X, A_solution, K_solution, logH, AlT, A_solid, K_solid);
        J = jacobianFB_Al(X, A_solution, K_solution, logH, A_solid, K_solid, SI_val, logC);

        deltaX = - J \ R;
        Xnew   = X + deltaX;
        if Xnew(2) < 0
            Xnew(2) = 0; 
        end

        if norm(Xnew - X) < tol
            X = Xnew;
            break
        end
        X = Xnew;
    end

    logAl3  = X(1);
    xcpVal  = X(2);

    % Reconstruct aqueous concentrations
    logC_final  = K_solution + A_solution*[logH; logAl3];
    Caq         = 10.^logC_final;
    C_all(:,ip) = Caq;
    xcp_vec(ip) = xcpVal;
    SI_vec(ip)  = SI_val;

    % Mass balance error
    sumAl = sumAlAq(Caq, A_solution) + xcpVal;
    massErr(ip) = abs(sumAl - AlT);

    % Pass solution as guess for next pH
    Xguess = X;
end

% Create the results struct (1D)
results = struct();
results.pH         = pHvals;
results.species    = species_solution;
results.Caq        = C_all;    % (nSol x nPts)
results.xcp        = xcp_vec; 
results.SI         = SI_vec;
results.massErr    = massErr;

% (6) Build the 2D surface: pH vs. AlT (with "chaining")
pH_grid  = 2 : 0.4 : 12;
AlT_grid = logspace(-6, -4, 10);
[PHmesh, AlTmesh] = meshgrid(pH_grid, AlT_grid);

massErrSurf = zeros(size(PHmesh));

nSol2  = size(A_solution,1);
CaqSurf= zeros([size(PHmesh), nSol2]);
xcpSurf= zeros(size(PHmesh));
SISurf = zeros(size(PHmesh));

for irow = 1:size(PHmesh,1)
    Xloc = [-3; 0]; % initial guess for each row
    for icol = 1:size(PHmesh,2)
        pHtemp   = PHmesh(irow, icol);
        logHtemp = -pHtemp;
        AlTemp   = AlTmesh(irow, icol);

        for iter = 1:50
            [Rloc, SI_loc, logCloc] = residualFB_Al_local(Xloc, A_solution, K_solution, logHtemp, AlTemp, A_solid, K_solid);
            Jloc = jacobianFB_Al_local(Xloc, A_solution, K_solution, logHtemp, A_solid, K_solid, SI_loc, logCloc);

            dX   = - Jloc \ Rloc;
            Xnew = Xloc + dX;
            if Xnew(2) < 0
                Xnew(2) = 0; 
            end

            if norm(Xnew - Xloc) < 1e-10
                Xloc = Xnew;
                break
            end
            Xloc = Xnew;
        end

        % Mass balance error
        logCtmp = K_solution + A_solution*[logHtemp; Xloc(1)];
        Caqtmp  = 10.^logCtmp;
        sumAl   = sumAlAq(Caqtmp, A_solution) + Xloc(2);
        massErrSurf(irow, icol) = abs(sumAl - AlTemp);

        % Store for 2D
        CaqSurf(irow, icol, :) = Caqtmp;
        xcpSurf(irow, icol)    = Xloc(2);
        SISurf(irow, icol)     = SI_loc;
    end
end

threshold = 1e-20;
isTiny    = (massErrSurf < threshold);

results2D = struct();
results2D.PHmesh      = PHmesh;
results2D.AlTmesh     = AlTmesh;
results2D.massErrSurf = massErrSurf;
results2D.threshold   = threshold;
results2D.isTiny      = isTiny;
results2D.CaqSurf     = CaqSurf;
results2D.xcpSurf     = xcpSurf;
results2D.SISurf      = SISurf;

end % end of main function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Auxiliary subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A_sol, K_sol, spc_sol, A_solid, K_solid, spc_solid] = parseAlTableau(Tableau)
    data  = Tableau(2:end,:);
    nRows = size(data,1);

    A_sol    = [];
    K_sol    = [];
    spc_sol  = {};
    A_solid  = [];
    K_solid  = [];
    spc_solid= {};

    for i=1:nRows
        rowStoichH  = data{i,1};
        rowStoichAl = data{i,3};
        valLogK     = data{i,5};
        phase       = data{i,6};
        spcName     = data{i,7};

        % skip e- and Cl-
        if strcmp(spcName,'e-') || strcmp(spcName,'Cl-')
            continue;
        end

        stoich2 = [rowStoichH, rowStoichAl];

        if phase==0
            A_sol   = [A_sol; stoich2];
            K_sol   = [K_sol; valLogK];
            spc_sol = [spc_sol; spcName];
        else
            A_solid   = [A_solid; stoich2];
            K_solid   = [K_solid; valLogK];
            spc_solid = [spc_solid; spcName];
        end
    end
end

function [R, SI, logC] = residualFB_Al(X, A, K, logH, AlT, A_solid, K_solid)
    % X(1) = logAl3, X(2) = xcp
    logAl3 = X(1);
    xcp    = X(2);

    logC = K + A*[logH; logAl3];
    Caq  = 10.^logC;

    sumAlAcuoso = sumAlAq(Caq, A);
    Rmass = sumAlAcuoso + xcp - AlT;

    Q  = A_solid*[logH; logAl3];
    SI = 10^(Q + K_solid);

    % Fischer-Burmeister
    [valFB,~,~] = fischerBurmeister(xcp, 1 - SI);

    R = [Rmass; valFB];
end

function J = jacobianFB_Al(X, A, K, logH, A_solid, K_solid, SI, logC)
    ln10     = log(10);
    Caq      = 10.^logC;
    stoichAl = A(:,2);

    dC_dlf       = ln10.*Caq.*stoichAl; 
    dRmass_dlf   = sum(dC_dlf);
    dRmass_dxcp  = 1;

    xcp          = X(2);
    [~, dFda, dFdb] = fischerBurmeister(xcp, 1 - SI);

    a_solid2   = A_solid(1,2);
    dSI_dlf    = SI * ln10 * a_solid2;  
    d1mSI_dlf  = - dSI_dlf;

    dRFB_dlf   = dFdb * d1mSI_dlf;
    dRFB_dxcp  = dFda;

    J = [
      dRmass_dlf,   dRmass_dxcp
      dRFB_dlf,     dRFB_dxcp
    ];
end

function s = sumAlAq(Caq, A)
    stoichAl = A(:,2);
    s = sum(stoichAl .* Caq);
end

function [val, dFda, dFdb] = fischerBurmeister(a,b)
    % FB(a,b) = sqrt(a^2 + b^2) - a - b
    r = sqrt(a^2 + b^2);
    val = r - a - b;

    epsVal = 1e-20;
    if r < epsVal
        val  = 0;
        dFda = -1;
        dFdb = -1;
    else
        dFda = (a/r) - 1;
        dFdb = (b/r) - 1;
    end
end

function [R, SI, logC] = residualFB_Al_local(X, A_sol, K_sol, logH, AlT, A_solid, K_solid)
    logAl3 = X(1);
    xcp    = X(2);

    logC = K_sol + A_sol*[logH; logAl3];
    Caq  = 10.^logC;

    sumAlAcuoso = sumAlAq(Caq, A_sol);
    Rmass = sumAlAcuoso + xcp - AlT;

    Q  = A_solid*[logH; logAl3];
    SI = 10^(Q + K_solid);

    [valFB,~,~] = fischerBurmeister(xcp, 1 - SI);
    R = [Rmass; valFB];
end

function J = jacobianFB_Al_local(X, A_sol, K_sol, logH, A_solid, K_solid, SI, logC)
    ln10     = log(10);
    Caq      = 10.^logC;
    stoichAl = A_sol(:,2);

    dC_dlf       = ln10.*Caq.*stoichAl;
    dRmass_dlf   = sum(dC_dlf);
    dRmass_dxcp  = 1;

    xcp          = X(2);
    [~, dFda, dFdb] = fischerBurmeister(xcp, 1 - SI);

    a_solid2   = A_solid(1,2);
    dSI_dlf    = SI * ln10 * a_solid2;
    d1mSI_dlf  = - dSI_dlf;

    dRFB_dlf   = dFdb * d1mSI_dlf;
    dRFB_dxcp  = dFda;

    J = [
        dRmass_dlf,   dRmass_dxcp
        dRFB_dlf,     dRFB_dxcp
    ];
end
