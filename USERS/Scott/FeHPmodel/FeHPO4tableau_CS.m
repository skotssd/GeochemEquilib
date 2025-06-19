function [Fep3, FeLOHR3s, HFO_P_ads, FeLPO4Rs, Fe3LPO4R2s, MASSERR,ErrCompName] = ...
         FeHPO4tableau_CS(pH, pe, TotVec, ...
                          flag0, flag1, flag2, flag3, flag4, flag5)
% =========================================================================
% FeHPO4tableau_CS  –  Especiación Fe–P–HFO usando PHREEQC
% -------------------------------------------------------------------------
% Devuelve (en mol kg-w, salvo MASSERR):
%   1) Fep3         Fe3+ acuoso
%   2) FeLOHR3s     Fe(OH)3(s)  (HFO)
%   3) HFO_P_ads    PO4 adsorbido total (sitios fuertes + débiles)
%   4) FeLPO4Rs     Fe(PO4)(s)  (strengite)
%   5) Fe3LPO4R2s   Fe3(PO4)2(s) (vivianita)
%   6) MASSERR      Máx. error de balance (signado)  [mol kg-w]
%
% TotVec  = [FeT  PT  ClT  NaT  HFOs  HFOw]   (fila o columna)
% =========================================================================

%% 1. ---------- Tableau ---------------------------------------------------
Tableau = {
% H+ e- Fe3 PO4 Cl Na Hfos Hfow logK   phase  species
 'H+' 'e-' 'Fe+3' 'PO4-3' 'Cl-' 'Na+' 'Hfos-' 'Hfow-' 'logK' 'phase' 'species';
  1    0    0      0       0     0     0       0       0      0   'H+';
  0    1    0      0       0     0     0       0       0      0   'e-';
  0    0    1      0       0     0     0       0       0      0   'Fe+3';
  0    0    0      1       0     0     0       0       0      0   'PO4-3';
  0    0    0      0       1     0     0       0       0      0   'Cl-';
  0    0    0      0       0     1     0       0       0      0   'Na+';
  0    0    0      0       0     0     1       0       0      0   'Hfos-';
  0    0    0      0       0     0     0       1       0      0   'Hfow-';
 -1    0    0      0       0     0     0       0    -14      0    'OH-';
  0    1    1      0       0     0     0       0    13.92     0   'Fe+2';
 -1    0    1      0       0     0     0       0    -2.19     0   'Fe(OH)+2';
 -2    0    1      0       0     0     0       0    -5.67     0   'Fe(OH)2+';
 -3    0    1      0       0     0     0       0   -12.56     0   'Fe(OH)3';
 -4    0    1      0       0     0     0       0   -21.60     0   'Fe(OH)4-';
 -2    0    2      0       0     0     0       0    -2.95     0   'Fe2(OH)2+4';
 -4    0    3      0       0     0     0       0    -6.30     0   'Fe3(OH)4+5';
 -1    1    1      0       0     0     0       0    -3.52     0   'Fe(OH)+';
 -2    1    1      0       0     0     0       0    -7.55     0   'Fe(OH)2';
 -3    1    1      0       0     0     0       0   -17.98     0   'Fe(OH)3-';
  1    0    0      1       0     0     0       0    12.346    0   'HPO4-2';
  2    0    0      1       0     0     0       0    19.553    0   'H2PO4-';
  3    0    0      1       0     0     0       0    21.721    0   'H3PO4';
  1    0    1      1       0     0     0       0    17.776    0   'FeHPO4+';
  2    0    1      1       0     0     0       0    24.983    0   'FeH2PO4+2';
  1    1    1      1       0     0     0       0    28.966    0   'FeHPO4';
  2    1    1      1       0     0     0       0    35.273    0   'FeH2PO4+';
  1    0    0      0       0     0     1       0     8.93     0   'HfosH';
  2    0    0      0       0     0     1       0    16.22     0   'HfosH2+';
  1    0    0      0       0     0     0       1     8.93     0   'HfowH';
  2    0    0      0       0     0     0       1    16.22     0   'HfowH2+';
  0    0    0      1       0     0     0       1    26.65     0   'HfowPO4-4';
  1    0    0      1       0     0     0       1    34.32     0   'HfowHPO4-3';
  2    0    0      1       0     0     0       1    40.32     0   'HfowH2PO4-2';
  0    0    0      1       0     0     1       0    26.65     0   'HfosPO4-4';
  1    0    0      1       0     0     1       0    34.32     0   'HfosHPO4-3';
  2    0    0      1       0     0     1       0    40.22     0   'HfosH2PO4-2';
 -3    0    1      0       0     0     0       0    -4.891    1   'Fe(OH)3s';
  0    0    1      1       0     0     0       0     23.2     1   'Fe(PO4)s';
  0    3    3      2       0     0     0       0     75.06    1   'Fe3(PO4)2s'
};

%% 2. ---------- Base de datos PHREEQC ------------------------------------
[totalnames, speciesexport, minerals] = writePHREEQCdatabase(Tableau);
totalnames    = totalnames';                      % vector columna
speciesexport = speciesexport(:,3:end)';          % quitar columna logK y fase

%% 3. ---------- Ejecutar PHREEQC -----------------------------------------
Temp = 25;                           % °C
acid = 'HCl'; pHfix = 1; show = 0;
[SPEC, SPECNAMES, SOLID, SOLIDNAMES] = ...
    writerunPHREEQC(Temp, pH, pe, ...
                    totalnames, TotVec, minerals, ...
                    speciesexport, 'DATABASE.dat', ...
                    show, acid, pHfix);
SPEC  = SPEC(:);        %  ← fuerza vector-columna n×1
SOLID = SOLID(:);       %  ← lo mismo para los sólidos

%% 4. ---------- Concentraciones de interés -------------------------------
Fep3        = concOf('Fep3');
FeLOHR3s    = concSolid('FeLOHR3s');     % Fe(OH)3(s)
FeLPO4Rs    = concSolid('FeLPO4Rs');     % Fe(PO4)(s)
Fe3LPO4R2s  = concSolid('Fe3LPO4R2s');   % Fe3(PO4)2(s)

ads_names   = {'HfowPO4m4','HfowHPO4m3','HfowH2PO4m2', ...
               'HfosPO4m4','HfosHPO4m3','HfosH2PO4m2'};
HFO_P_ads   = 0;
for ii = 1:numel(ads_names)
    HFO_P_ads += concOf(ads_names{ii});
end

%% 5. ---------- Balance de masa (6 componentes, descarta Cl) -------------
%   Usar las mismas filas que usa writePHREEQCdatabase:
phases   = Tableau(4:end,end-1);                 % filas 4..end
Aall     = cell2mat(Tableau(4:end,3:end-3));     % este bloque es 6 comp.

sol_idx  = find(cell2mat(phases)==0);
soli_idx = find(cell2mat(phases)==1);
Asol     = Aall;  Asol(soli_idx,:) = [];
Asolid   = Aall;  Asolid(sol_idx,:) = [];

mass_soln  = Asol'   * SPEC;        % (6×34)*(34×1)
mass_solid = Asolid' * SOLID;       % (6×2 )*(2 ×1)
Tcalc      = mass_soln + mass_solid;

Totals         = TotVec(:);         % columna
Totals(3)      = [];                % quitar Cl
Tcalc(3)       = [];

% --- Identificar componente del error máximo ---
% Define el orden de los componentes en el vector 'delta'
compNames = {'Fe', 'P', 'Na', 'Hfos', 'Hfow'}; % Coincide con el orden después de quitar Cl

delta        = Totals - Tcalc;
[abs_err, k] = max(abs(delta)); % 'k' es el índice del error máximo

MASSERR      = abs_err * sign(delta(k)); % Mantiene el error con signo
ErrCompName  = compNames{k}; % Obtiene el nombre del componente usando el índice 'k'

%% ------------------ funciones internas ----------------------------------
function v = concOf(nm)
    idx = strcmp(SPECNAMES,nm);
    if any(idx), v = SPEC(idx); else, v = 0; end
end
function v = concSolid(nm)
    idx = strcmp(SOLIDNAMES,nm);
    if any(idx), v = SOLID(idx); else, v = 0; end
end

end  % ================================ FIN FUNCIÓN ========================
