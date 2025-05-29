clear; clc; close all;

%% 1) Inputs
AlT = 1.0e-5;
pHvals = 2 : 0.1 : 12;

%% 2) Llamar a la funcion que realiza los calculos
[results, results2D] = Al_solid_complementarity_FB(AlT, pHvals);

%% 3) Figuras (parte 1): Barrido 1D con subplots

figure('Name','1D: pH vs. speciation, 4 subplots');

% (a) Algunas especies acuosas vs pH
subplot(2,2,1)
iAl3   = find(strcmp(results.species,'Al^{3+}'));
iAlOH4= find(strcmp(results.species,'Al(OH)4-'));
iAl13 = find(strcmp(results.species,'Al13(OH)32^{7+}'));

if ~isempty(iAl3)
    plot(results.pH, results.Caq(iAl3,:)./AlT, '-b','LineWidth',2,'DisplayName','Al^{3+}'); 
    hold on
end
if ~isempty(iAlOH4)
    plot(results.pH, results.Caq(iAlOH4,:)./AlT, '-m','LineWidth',2,'DisplayName','Al(OH)4^-');
end
if ~isempty(iAl13)
    plot(results.pH, results.Caq(iAl13,:)./AlT, '-r','LineWidth',2,'DisplayName','Al_{13}');
end
xlabel('pH'); ylabel('fraction Al');
legend('Location','northeastoutside'); grid on
title('(a) Main aqueous species');

% (b) Solido
subplot(2,2,2)
plot(results.pH, results.xcp./AlT, '-o','LineWidth',1.5);
xlabel('pH'); ylabel('frac. Al in solid'); 
title('(b) Al(OH)_3(s)'); grid on

% (c) Error de masa
subplot(2,2,3)
semilogy(results.pH, results.massErr, '-o','LineWidth',1.5);
xlabel('pH'); ylabel('Mass balance error (M)');
title('(c) mass error'); grid on

% (d) Saturation Index
subplot(2,2,4)
plot(results.pH, results.SI, '-o','LineWidth',1.5);
xlabel('pH'); ylabel('SI(Al(OH)3)');
title('(d) saturation index'); grid on


%% 4) Figura semilogy de concentraciones acuosas y el sólido

figure('Name','Semilogy comparison of species');
iAl3plus    = find(strcmp(results.species,'Al^{3+}'));
iAlOH3plus  = find(strcmp(results.species,'Al(OH)^{2+}'));
iAlOH4minus = find(strcmp(results.species,'Al(OH)4-'));
iAl13       = find(strcmp(results.species,'Al13(OH)32^{7+}'));
iAl2        = find(strcmp(results.species,'Al2(OH)2^{4+}'));
iAl3        = find(strcmp(results.species,'Al3(OH)4^{5+}'));

semilogy(results.pH, results.Caq(iAl3plus,:), ...
         results.pH, results.Caq(iAlOH3plus,:), ...
         results.pH, results.Caq(iAlOH4minus,:), ...
         results.pH, results.Caq(iAl13,:), ...
         results.pH, results.Caq(iAl2,:), ...
         results.pH, results.Caq(iAl3,:), ...
         results.pH, results.xcp, 'LineWidth',1.2);
legend('Al^{3+}','Al(OH)^{2+}','Al(OH)4^-','Al_{13}','Al_{2}','Al_{3}','Al(OH)_3(s)', ...
       'Location','northeastoutside');
xlabel('pH'); ylabel('Concentration (M)');
title('Semilogy: Comparison of species concentrations');
ylim([1E-30 1e-4])
grid on;


%% 5) Superficie 3D con log10(Error)

figure('Name','Surface 3D: pH vs AlT vs log10(error)');
surf(results2D.PHmesh, log10(results2D.AlTmesh), log10(results2D.massErrSurf));
xlabel('pH');
ylabel('log_{10}(AlT)');
zlabel('log_{10}(massErr)');
title('Surface pH vs AlT (error)'); 
grid on

%% 6) Superficie con zona gris para errores < 1e-20
threshold = results2D.threshold;  % 1e-20
isTiny = results2D.isTiny;

Z_forLog = results2D.massErrSurf;
Z_forLog(isTiny) = threshold;
Zplot = log10(Z_forLog);

Xmesh = results2D.PHmesh;          % 2D
Ymesh = log10(results2D.AlTmesh);  % 2D, convert AlT to log10

figure('Name','Surface log10(error) w/ gray zone <1e-20');
surf(Xmesh, Ymesh, Zplot);
shading interp
colormap('jet')
colorbar
xlabel('pH');
ylabel('log_{10}(AlT)');
zlabel('log_{10}(error)');
title('Cells <1e-20 in gray');
hold on

% Celdas tiny en gris
Z2 = Zplot;
Z2(~isTiny) = NaN; 
surf(Xmesh, Ymesh, Z2, ...
     'FaceColor',[0.6 0.6 0.6], ...
     'EdgeColor','none');

grid on