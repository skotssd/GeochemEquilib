function [Alp3,AlLOHR3s,Al13LOHR32p7,AlLOHR4m,MASSERR]=AltableauPICKER(pH,pe,T,flag0,flag1,flag2,flag3,flag4,flag5)

% input tableau.  change this part % ----------------------------------------------

Tableau=[...
{'H+'} {'e-'} {'Al+3'}     {'Cl-'}  {'logK'}                    {'phase'}    {'species'} 
1       0        0          0        0                           0          {'H+'}
0       1        0          0        0                           0          {'e-'}
0       0        1          0        0                           0          {'Al+3'}
0       0        0          1        0                           0          {'Cl-'}
-1      0        0          0        -13.9951                    0          {'OH-'}
-2      0        1          0        -10.5945                    0          {'Al(OH)2+'}
-1      0        1          0        -4.9571                     0          {'Al(OH)+2'}
-4      0        1          0        -22.7                       0          {'Al(OH)4-'}
-2      0        2          0        -7.7902                     0          {'Al2(OH)2+4'}
-4      0        3          0        -13.8803                    0          {'Al3(OH)4+5'}
-32     0        13         0        -98.73                      0          {'Al13(OH)32+7'}
%solids
-3       0       1          0        -7.7560                     1          {'Al(OH)3s'}
];


% end of tableau.  ------------------ % ----------------------------------------------

if flag0==0 % solve with tableau (using other flags for options)

    [KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES]=processtableau(Tableau(2:end,1:end),pH,pe);
    [SPECIESCONCS,SPECIATIONNAMES,MASSERR,X]=returnspeciationRE(KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES,T,flag1,flag2,flag3,flag4,flag5);
    for k=1:size(SPECIESCONCS,1)
          txt=[SPECIATIONNAMES(k,:),'=SPECIESCONCS(k);'];
          eval(txt)
    end

    MASSERR=max(MASSERR);
end

if flag0==1 %solve using PHREEQC

    %write databasefile

    [totalnames,speciesexport,minerals]=writePHREEQCdatabase(Tableau); %writes "DATABASE.dat"
    totalnames=totalnames'; speciesexport=speciesexport(:,3:end); speciesexport=speciesexport';
    totalvector=T; Temp=25; database=['DATABASE.dat']; show=0; acid=['HCl']; pHfixedcheck=1;
    
    [SPECIESCONCS, SPECIATIONNAMES, SOLIDconcs, SOLIDnames]=...
    runPHREEQCv3(Temp,pH,pe,totalnames,totalvector,minerals,speciesexport,database,show,acid,pHfixedcheck);

    SPECIATIONNAMES = regexprep(SPECIATIONNAMES,'[m_]','');
    SPECIATIONNAMES = regexprep(SPECIATIONNAMES,'[(]','L');
    SPECIATIONNAMES = regexprep(SPECIATIONNAMES,'[)]','R');
    SPECIATIONNAMES = regexprep(SPECIATIONNAMES,'[+]','p');
    SPECIATIONNAMES = regexprep(SPECIATIONNAMES,'[-]','m');

    SOLIDnames = regexprep(SOLIDnames,'[(]','L');
    SOLIDnames = regexprep(SOLIDnames,'[)]','R');
       
    for k=1:size(SPECIESCONCS,1)
          txt=[cell2mat(SPECIATIONNAMES(k,:)),'=SPECIESCONCS(k,:);'];
          eval(txt)
    end

    for k=1:size(SOLIDconcs,1)
          txt=[cell2mat(SOLIDnames(k,:)),'=SOLIDconcs(k,:);'];
          eval(txt)
    end

    %determine total masses
    %parse the tableau
    phases=Tableau(4:end,end-1); c1=0; c2=0;
    for i=1:length(phases)
        tst=cell2mat(phases(i));
        if tst==0; c1=c1+1; solnindex(c1)=i; end
        if tst==1; c2=c2+1; solidindex(c2)=i; end
    end
    Aall=cell2mat(Tableau(4:end,3:end-3));
    Asoln=Aall; Asolid=Aall;
    Asoln(solidindex,:)=[];
    Asolid(solnindex,:)=[];
    masssoln=Asoln'*SPECIESCONCS; 
    % ignore the chloride error b/c adjusted to fix pH
    MASSERR=T(1)-(masssoln(1)+AlLOHR3s);

    %Rmass=Asolution'*C+Asolid'*Xsolidzero-T;
    
end

end % end of function

% ---------------- SUBFUNCTIONS --------------------------------------------------------
