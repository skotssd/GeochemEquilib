function [solutionspeciesconcs, speciesnames, SOLIDconcs, SOLIDnames]=...
    writerunPHREEQC(T,pH,pe,totalnames,totalvector,minerals,speciesexport,database,show,acid,pHfixedcheck);

NOOFSOLIDS=size(minerals,1); %also CO2(g) if present

% Construct the text file to run PHREEQC from MATLAB-------------------
fileID=fopen('runphreeqc.txt','w');
fclose(fileID);
fileID=fopen('runphreeqc.txt','a');

% add species that website says can help with convergence
fprintf(fileID,'SOLUTION_SPECIES\n');
fprintf(fileID,'H2O + 0.01e- = H2O-0.01\n');
fprintf(fileID,'log_k   -9.0\n');
fprintf(fileID,'\n');

% define solution  -------------------------------------------------------
fprintf(fileID,'SOLUTION 1\n');
fprintf(fileID,['       pe      ' num2str(pe), '\n']); % the redox value
fprintf(fileID,['       pH      ' num2str(pH), '\n']); % pH value
fprintf(fileID,['       temp      ' num2str(T), '\n']); % the temperature
fprintf(fileID,'-units mol/kgw\n'); % the unit of the input; usually mol/L is used
% put in the totals

for i=1:size(totalnames,1)
    % tst means no charge balance species, but could set Na or Cl as the charge balance spcies. no idea the advantage.
    tf=strcmp('tst',cell2mat(totalnames(i))); 
    if tf==1
        totaltxt=[cell2mat(totalnames(i)),' ', num2str(totalvector(i)), ' charge\n'];
        fprintf(fileID,totaltxt);
    end
    if tf==0
        element=cell2mat(totalnames(i)); element= regexprep(element,'[O]','');
    totaltxt=[element,' ', num2str(totalvector(i)), ' \n'];
    fprintf(fileID,totaltxt);
    end
end
fprintf(fileID,'END');
fprintf(fileID,'\n');
fprintf(fileID,'USE solution 1');
fprintf(fileID,'\n');

% define numerical solver options
fprintf(fileID,'KNOBS\n');
fprintf(fileID,'-iterations 200\n');
fprintf(fileID,'\n');

% define equilibrium phases (solids, fix pH, pe) -------------------------
fprintf(fileID,'EQUILIBRIUM_PHASES 1\n'); 
for i=1:NOOFSOLIDS
    tst=cell2mat(minerals(i)); T=startsWith(minerals(i),"CO2(g)"); 
    if T==1
    phasestxt=[cell2mat(minerals(i)), '   -3.5\n'];
    else
    phasestxt=[cell2mat(minerals(i)), '   0.0   0\n'];
    end
    fprintf(fileID,phasestxt);
 end


if pHfixedcheck==1; 
    pHfixline=['       Fix_H+ -',num2str(pH),'          ',acid,' 10.0\n']; 
    fprintf(fileID,pHfixline);
    fprintf(fileID,'-force_equality true\n');
    pefixline=['       Fix_pe ',num2str(-1*pe),'          O2\n']; 
    fprintf(fileID,pefixline);
    fprintf(fileID,'-force_equality true\n');
    fprintf(fileID,'\n');
end

% define outputs of model ------------------------------------------------
fprintf(fileID,'SELECTED_OUTPUT\n');
fprintf(fileID,'-file selected.out\n');
fprintf(fileID,'-selected_out true\n');
fprintf(fileID,'-user_punch true\n');
fprintf(fileID,'-high_precision true\n');
fprintf(fileID,'-reset false\n');
fprintf(fileID,'-simulation false\n');
fprintf(fileID,'-state false\n');
fprintf(fileID,'-distance false\n');
fprintf(fileID,'-time false\n');
fprintf(fileID,'-step false\n');
fprintf(fileID,'-ph false\n');
fprintf(fileID,'-pe false\n');
fprintf(fileID,'-reaction false\n');
fprintf(fileID,'-temperature false\n');
fprintf(fileID,'-alkalinity false\n');
fprintf(fileID,'-ionic_strength false\n');
fprintf(fileID,'-water false\n');
fprintf(fileID,'-charge_balance false\n');
fprintf(fileID,'-percent_error false\n');

% define USERPUNCH outputs of model ------------------------------------------------
fprintf(fileID,'USER_PUNCH\n'); writestr1=['-headings ']; writestr2=['10 PUNCH '];
for i=1:length(speciesexport)
    writestr1=[writestr1, cell2mat(speciesexport(i)),"\t"];
    writestr2=[writestr2, "MOL(""",cell2mat(speciesexport(i)),"""),"];
end
for i=1:length(minerals)
    writestr1=[writestr1, cell2mat(minerals(i)),"\t"];
    writestr2=[writestr2, "EQUI_DELTA(""",cell2mat(minerals(i)),"""),"];
end
writestr1=[writestr1,'\n']; writestr2=[writestr2,'\n'];
fprintf(fileID,writestr1);
fprintf(fileID,'-start\n');
fprintf(fileID,writestr2);
fprintf(fileID,'-end\n');
fclose(fileID);
 
% run the model -----------------------------------------------------------
str=['system("phreeqc runphreeqc.txt out.txt ', database,'");'];
if show==1;eval(str); end % output to the screen 
if show==0;evalc(str); end % so no screen output

%parse the outputs of the model
%evalc(str); % so no screen output
fid = fopen('selected.out','rt');
hdr = strtrim(regexp(fgetl(fid),'\t','split'));
hdr=hdr(1:end-1)';
mat = cell2mat(textscan(fid,repmat('%f',1,numel(hdr))));

speciesnames= regexprep(speciesexport,'[(]','L');
speciesnames = regexprep(speciesnames,'[)]','R');
speciesnames = regexprep(speciesnames,'[+]','p');
speciesnames = regexprep(speciesnames,'[-]','m');
[n,m]=size(speciesexport);
solutionspeciesconcs=mat(1:n); SOLIDconcs=mat(n+1:end);
SOLIDnames = regexprep(minerals,'[(]','L');
SOLIDnames = regexprep(SOLIDnames,'[)]','R');

end % end of function
