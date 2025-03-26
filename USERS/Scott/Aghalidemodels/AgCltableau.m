function [Agp,Clm,AgCls,MASSERR]=AgCltableau(pH,pe,T,flag0,flag1,flag2,flag3,flag4,flag5)

% input tableau.  change this part % ----------------------------------------------

Tableau=[...
{'H+'} {'e-'} {'Ag+'}  {'Na+'}   {'Cl-'}  {'logK'}                    {'phase'}    {'species'} 
1       0        0      0           0        0                           0          {'H+'}
0       1        0      0           0        0                           0          {'e-'}
0       0        1      0           0        0                           0          {'Ag+'}
0       0        0      1           0        0                           0          {'Na+'}
0       0        0      0           1        0                           0          {'Cl-'}
-1      0        0      0           0        -14                         0          {'OH-'}
0       0        1      0           1        -3                          0          {'AgCl'}
-1      0        1      0           0        -12                         0          {'Ag(OH)'}
%solids
0       0        1      0           1        9.25                        1          {'AgCls'}
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
    totalvector=T; Temp=25; database=['DATABASE.dat']; show=0; acid=['HX']; pHfixedcheck=1;
    
    [SPECIESCONCS, SPECIATIONNAMES, SOLIDconcs, SOLIDnames]=...
    writerunPHREEQC(Temp,pH,pe,totalnames,totalvector,minerals,speciesexport,database,show,acid,pHfixedcheck);
    SPECIESCONCS=SPECIESCONCS'; SOLIDconcs=SOLIDconcs';
       
    for k=1:size(SPECIESCONCS,1)
          txt=[cell2mat(SPECIATIONNAMES(k,:)),'=SPECIESCONCS(k,:);'];
          eval(txt)
    end

    for k=1:size(SOLIDconcs',1)
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
    masssolid=Asolid'*SOLIDconcs;
    Tcalc=masssoln+masssolid;
    MASSERR=max(T-Tcalc);
end

end % end of function

% ---------------- SUBFUNCTIONS --------------------------------------------------------
