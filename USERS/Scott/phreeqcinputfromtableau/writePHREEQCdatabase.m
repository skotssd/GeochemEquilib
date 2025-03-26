function [totalnames,speciesnames,minerals]=writePHREEQCdatabase(Tableau)

%delete existing phreeqc files
system("rm DATABASE.dat");
%create empty phreeqc files
system("touch DATABASE.dat");

totalnames=[]; speciesnames=[]; % initialize outputs

% parse the tableau
components=Tableau(1,1:size(Tableau,2)-3); noofcomponents=length(components);
elements=Tableau(1,3:size(Tableau,2)-3);
elements = regexprep(elements,'[\d"]','');
elements = regexprep(elements,'[+]',''); totalnames=elements;
elements = regexprep(elements,'[-]',''); totalnames=elements;
species=Tableau(2:end,size(Tableau,2)); 
phases=Tableau(2:end,size(Tableau,2)-1);
logK=Tableau(2:end,size(Tableau,2)-2);
c1=0; c2=0; %two counters
for i=1:length(species)
    tst=cell2mat(phases(i));
    if tst==0; c1=c1+1; solutionspecies(c1)=species(i); solnindex(c1)=i; end
    if tst==1; c2=c2+1; solidspecies(c2)=species(i); solidindex(c2)=i; end
end
minerals=solidspecies;
speciesnames=solutionspecies;
Aall=cell2mat(Tableau(2:end,1:end-3));
Asoln=Aall; Asolid=Aall;
Asoln(solidindex,:)=[]; 
Asolid(solnindex,:)=[];
logKsoln=logK; logKsolid=logK;
logKsoln(solidindex,:)=[];
logKsolid(solnindex,:)=[];

% ---------- write masterspecies -----------------------------------------

MASTERSPECIES=[...
{'SOLUTION_MASTER_SPECIES\n'}
{'  \n'}
{'#element species        alk     gfw_formula     element_gfw\n'}
{'  \n'}
{'H        H+             -1.     H               1.008\n'}
{'H(0)     H2             0.0     H\n'}
{'H(1)     H+             -1.     0.0\n'}
{'E        e-             0.0     0.0             0.0\n'}
{'O        H2O            0.0     O               16.00\n'}
{'O(0)     O2             0.0     O\n'}
{'O(-2)    H2O            0.0     0.0\n'}
{'Cl       Cl-            0.0     Cl              35.4527\n'}
{'Cl(-1)   Cl-            0       Cl\n'}
{'  \n'}
];
fileID=fopen('DATABASE.dat','a');
for i=1:size(MASTERSPECIES,1) % H+ and e- already written
    line=cell2mat(MASTERSPECIES(i,:));
    fprintf(fileID,line)
end
% add components to the master species

for i=3:length(components) % H+ and e- already written
    component=cell2mat(components(:,i));
    %extractbefore.m would work great but not in Octave yet
    for j=1:length(component)
        tst=component(j);
        if tst=="+"; break; end
        if tst=="-"; break; end
        componentnocharge(:,j)=tst; 
    end
    line=[componentnocharge,'\t\t',component,'\t\t\t  0 \t 1 \t\t\t\t 1'];
    % keep track of the element names
    fprintf(fileID,line)
    fprintf(fileID,' \n')
end
fprintf(fileID,' \n');
fprintf(fileID,' \n');
fclose(fileID);

%------------- write solution species ------------------------------------------

SOLUTIONSPECIES=[...
{'SOLUTION_SPECIES\n'}
{'# use -gamma 1e10 0.0 to counteract automatic activity correction\n'}
{'  \n'}
{'H+ = H+\n'}
{'      log_k           0.000\n'}
{'      -gamma          1e10     0.0\n'}
{'  \n'}
{'e- = e-\n'}
{'        log_k           0.000\n'}
{'  \n'}
{'H2O = H2O\n'}
{'        log_k           0.000\n'}
{'  \n'}
{'Cl- =  Cl- \n'}
{'    log_k 0\n'}
{'  \n'}
{'H2O = OH- + H+\n'}
{'        log_k           -14.000\n'}
{'        -gamma          1e10     0.0\n'}
{'  \n'}
{'2 H2O = O2 + 4 H+ + 4 e-\n'}
{'        log_k           -86.08\n'}
{'        -gamma          1e10     0.0\n'}
{'  \n'}
{'2 H+ + 2 e- = H2\n'}
{'        log_k           -3.15\n'}
{'        -gamma          1e10     0.0\n'}
{'  \n'}
];
fileID=fopen('DATABASE.dat','a');
for i=1:size(SOLUTIONSPECIES,1)
    line=cell2mat(SOLUTIONSPECIES(i,:));
    fprintf(fileID,line)
end
% add components rxns the solutionspecies
for i=3:length(components)
    component=cell2mat(components(:,i));
    line=[component,'=',component,'\n'];
    fprintf(fileID,line)
    line=['log_k \t 0.000\n'];
    fprintf(fileID,line)
    line=['-gamma  1e10 0.0\n'];
    fprintf(fileID,line)
    fprintf(fileID,' \n');
end
% add species rxns
for i=noofcomponents+2:size(Asoln,1) % adjust so to not include identity at top, or OH-
    test=Asoln(i,:);
    writestr=[];
    for j=1:length(test)
        %only negative will be hydroxide in normal tableau writing because use formation constants
        if test(j)<0; Astr=num2str(abs(Asoln(i,j)));
        writestr=[Astr, " H2O "];
        end
        if test(j)>0; Astr=num2str(abs(Asoln(i,j))); writestr=[writestr, " + ", Astr, cell2mat(components(j))];
        end
    end
    %for j=1:length(test)
    writestr=[writestr, " = ", cell2mat(solutionspecies(i))];
    %only negative will be hydroxide in normal tableau writing because use formation constants
    if test(1)<=-1; Astr=num2str(abs(Asoln(i,1)));
        writestr=[writestr," + ",Astr, " H+ "];
    end
    fprintf(fileID,writestr);
    fprintf(fileID,' \n');
    Kstr=num2str(cell2mat(logKsoln(i)));
    line=['log_k\t',Kstr,'\n'];
    fprintf(fileID,line)
    line=['-gamma  1e10 0.0\n'];
    fprintf(fileID,line)
    fprintf(fileID,' \n')
end
fclose(fileID);

% ---- write solids --------------------------------------------------------------------------

% add solids
fileID=fopen('DATABASE.dat','a');
fprintf(fileID,'PHASES')
fprintf(fileID,' \n')
fprintf(fileID,'Fix_H+\n')
fprintf(fileID,'\t H+ = H+ ; log_k 0.0\n')
fprintf(fileID,'Fix_pe\n')
fprintf(fileID,'\t e- = e- ; log_k 0.0\n')

for i=1:size(Asolid,1)

    writestr=[cell2mat(solidspecies(i))];
    fprintf(fileID,writestr);
    fprintf(fileID,' \n');

    writestr=[];
    solidstr=cell2mat(solidspecies(i)); solidstr=solidstr(1:end-1);
    writestr=solidstr;
    %only negative will be hydroxide in normal tableau writing because use formation constants
    if test(1)<-1; Astr=num2str(abs(Asolid(i,1)));
        writestr=[writestr," + ",Astr, " H+ "];
    end
    writestr=[writestr,' = '];
     
    test=Asolid(i,:);
    for j=1:length(test)
        %only negative will be hydroxide in normal tableau writing because use formation constants
        if test(j)<0; Astr=num2str(abs(Asolid(i,j)));
        writestr=[writestr, Astr, "H2O "];
        end
        if test(j)>0; writestr=[writestr, " + ", cell2mat(components(j))];
        end
    end
    fprintf(fileID,writestr);    
    fprintf(fileID,' \n');
    % put in as dissosciation but tableau written as association
    logKval=-1*cell2mat(logKsolid(i)); Kstr=num2str((logKval)); 
    line=['log_k\t',Kstr,'\n'];
    fprintf(fileID,line)
    fprintf(fileID,' \n')
end

% write output file
%punchboxtext=[...
%{'USER_PUNCH\n'}
%{'        -headings seconds  Ca  pH  \n'}
%{'  10 PUNCH SIM_TIME, TOT("Al"), -LA("H+")\n'}
%];
%[nolinesPUNCH,length]=size(punchboxtext);

fclose(fileID);

end