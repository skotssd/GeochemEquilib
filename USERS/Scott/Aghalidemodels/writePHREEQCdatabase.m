function [totalnames,speciesnames,minerals]=writePHREEQCdatabase(Tableau)

%delete existing phreeqc files
system("rm DATABASE.dat");
%create start of DATABASE (the possible components part) from DEFAULT.dat
system("cp DEFAULT.dat DATABASE.dat");

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
% ----------- not necessary, already in DEFAULT.dat ----------------------
%------------- write solution species ------------------------------------------

fileID=fopen('DATABASE.dat','a');
% add species rxns
for i=noofcomponents+2:size(Asoln,1) % adjust so to not include identity at top, or OH-
    test=Asoln(i,:);
    LSstr=[]; RSstr=[]; % just write the reactants if positive, the products if negative
    % if first entry is negative treat as hydroxide
    if test(1)<0; Astr=num2str(abs(test(1))); LSstr=[Astr, "H2O "]; RSstr=[Astr," H+ "]; 
        for j=2:length(test)
            if test(j)>0; Astr=num2str(Asoln(i,j)); LSstr=[LSstr, " + ", Astr, cell2mat(components(j))]; end
        end
    end
    if test(1)>=0;
        for j=2:length(test)
            if test(j)>0; Astr=num2str(Asoln(i,j));
                LSstr=[LSstr, " + ", Astr, cell2mat(components(j))]; end
        end
    end
    RSstr=[" + ", cell2mat(solutionspecies(i)), " + ",RSstr];
    % remove leading + on LS and RS
    if LSstr(1:3)==" + "; LSstr=LSstr(4:end); end
    if RSstr(1:3)==" + "; RSstr=RSstr(4:end); end
    writestr=[LSstr, " = ", RSstr];
    if writestr(end-2:end)==" + "; writestr=writestr(1:end-3); end
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

    writestr=[]; LSstr=[]; RSstr=[];
    solidstr=cell2mat(solidspecies(i)); solidstr=solidstr(1:end-1);
    LSstr=solidstr;
    %only negative will be hydroxide in normal tableau writing because use formation constants
    if test(1)<-1; Astr=num2str(abs(Asolid(i,1)));
        LSstr=[LSstr," + ",Astr, " H+ "];
        RSstr=[RSstr, Astr, "H2O "];
    end
    test=Asolid(i,:);
    for j=2:length(test) % ignore H+ b/c taken care of already with test(1)<0 test
        %only negative will be hydroxide in normal tableau writing because use formation constants
        if test(j)>0; Astr=num2str(abs(Asolid(i,j))); RSstr=[RSstr, " + ", Astr,cell2mat(components(j))];
        end
    end
    if RSstr(1:3)==" + "; RSstr=RSstr(4:end); end
    writestr=[LSstr," = ", RSstr];
    fprintf(fileID,writestr);    
    fprintf(fileID,' \n');
    % put in as dissosciation but tableau written as association
    logKval=-1*cell2mat(logKsolid(i)); Kstr=num2str((logKval)); 
    line=['log_k\t',Kstr,'\n'];
    fprintf(fileID,line)
    fprintf(fileID,' \n')
end

fclose(fileID);

end % end of function