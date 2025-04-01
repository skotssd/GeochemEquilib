function [Fep3,FeLOHR3s,HFOP,MASSERR]=FeHPO4tableau(pH,pe,T,flag0,flag1,flag2,flag3,flag4,flag5)

HFOsT=T(5); HFOwT=T(6);

% input tableau.  change this part % ----------------------------------------------

Tableau=[...
{'H+'} {'e-'} {'Fe+3'}  {'PO4-3'}   {'Cl-'}  {'Na+'}  {'Hfos-'}  {'Hfow-'}  {'logK'}                    {'phase'}    {'species'} 
1       0        0      0           0          0       0            0     0                           0          {'H+'}
0       1        0      0           0          0       0            0      0                           0          {'e-'}
0       0        1      0           0          0       0            0     0                           0          {'Fe+3'}
0       0        0      1           0          0       0            0     0                           0          {'PO4-3'}
0       0        0      0           1          0       0            0     0                           0          {'Cl-'}
0       0        0      0           0          1       0            0     0                           0          {'Na+'}
0       0        0      0           0          0       1            0     0                           0          {'Hfos-'}
0       0        0      0           0          0       0            1     0                           0          {'Hfos+'}
% Kw
-1      0        0      0           0          0       0            0   -14                           0          {'OH-'}
% hydroloysis/redox products
0       1        1      0           0          0       0            0     13.92                       0          {'Fe+2'}
-1      0        1      0           0          0       0            0     -2.19                       0          {'Fe(OH)+2'}
-2      0        1      0           0          0       0            0     -5.67                       0          {'Fe(OH)2+'}
-3      0        1      0           0          0       0            0     -12.56                      0          {'Fe(OH)3'}
-4      0        1      0           0          0       0            0     -21.6                       0          {'Fe(OH)4-'}
-2      0        2      0           0          0       0            0     -2.95                       0          {'Fe2(OH)2+4'}
-4      0        3      0           0          0       0            0     -6.3                        0          {'Fe3(OH)4+5'}
-1      1        1      0           0          0       0            0     -3.52                       0          {'Fe(OH)+'}
-2      1        1      0           0          0       0            0     -7.55                       0          {'Fe(OH)2'}
-3      1        1      0           0          0       0            0     -17.98                      0          {'Fe(OH)3-'}
%phosphoric acid 
1       0        0      1           0        0         0            0     12.346                      0          {'HPO4-2'}
2       0        0      1           0        0         0            0     19.553                      0          {'H2PO4-'}
3       0        0      1           0        0         0            0     21.721                      0          {'H3PO4'}
% Fe-P soluble complexes
1        0       1      1           0        0         0            0     17.776                      0          {'FeHPO4+'}
2        0       1      1           0        0         0            0     24.983                      0          {'FeH2PO4+2'}
1        1       1      1           0        0         0            0     28.966                      0          {'FeHPO4'}
2        1       1      1           0        0         0            0     35.273                      0          {'FeH2PO4+'}
% HFO  SCM reactions
1        0       0      0           0        0         1            0     8.93                        0          {'HfosH'}
2        0       0      0           0        0         1            0     16.22                       0          {'HfosH2+'}
1        0       0      0           0        0         0            1     8.93                        0          {'HfowH'}
2        0       0      0           0        0         0            1     16.22                       0          {'HfowH2+'}
0        0       0      1           0        0         0            1     26.65                       0          {'HfowPO4-4'}
1        0       0      1           0        0         0            1     34.32                       0          {'HfowHPO4-3'}
2        0       0      1           0        0         0            1     40.32                       0          {'HfowH2PO4-2'}
0        0       0      1           0        0         1            0     26.65                       0          {'HfosPO4-4'}
1        0       0      1           0        0         1            0     34.32                       0          {'HfosHPO4-3'}
2        0       0      1           0        0         1            0     40.22                       0          {'HfosH2PO4-2'}
%solids
-3      0        1      0           0        0         0            0     -4.891                      1          {'Fe(OH)3s'} %HFO
0       0        1      1           0        0         0            0     23.2                        1          {'Fe(PO4)s'} % strengite
0       3        3      2           0        0        0            0      75.06                       1          {'Fe3(PO4)2s'} %vivianite
];

Tableau=Tableau(2:end,1:end); INDEX=[];

if HFOsT==0 %7th entry
    Nentries=cell2mat(Tableau(:,7));
    indices = find(Nentries>=1);
    Tableau(indices,:) = [];
    INDEX=[INDEX 7];
    HFOP=0;
end

if HFOwT==0 %8th entry
    Nentries=cell2mat(Tableau(:,8));
    indices = find(Nentries>=1);
    Tableau(indices,:) = [];
    INDEX=[INDEX 8];
    HFOP=0;
end

Tableau(:,INDEX)=[];

T(T==0)=[];

% end of tableau.  ------------------ % ----------------------------------------------

if flag0==0 % solve with tableau (using other flags for options)

%[KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES]=processtableau(Tableau(2:end,1:end),pH,pe);
[KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES]=processtableau(Tableau,pH,pe);
    [SPECIESCONCS,SPECIATIONNAMES,masserr,X]=returnspeciationRE(KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES,...
    T,flag1,flag2,flag3,flag4,flag5);
    for k=1:size(SPECIESCONCS,1)
          txt=[SPECIATIONNAMES(k,:),'=SPECIESCONCS(k);'];
          eval(txt)
    end
    
    %[MASSERR,IDX]=max((MASSERR)); % oops. this would ignore negative errors
    %MASSERR=max(100*((MASSERR)./T)); % relative does not make sense when some components are so tiny (like HFOinitial conc)
    [MASSERR,INDX]=max(abs(masserr)); MASSERR=MASSERR*sign(masserr(INDX));
end

if flag0==1 %solve using PHREEQC

    %write databasefile

    [totalnames,speciesexport,minerals]=writePHREEQCdatabase(Tableau);
        totalnames=totalnames'; speciesexport=speciesexport(:,3:end); speciesexport=speciesexport';
    totalvector=T; Temp=25; database=['DATABASE.dat']; show=0; acid=['HCl']; pHfixedcheck=1;
    
    [SPECIESCONCS, SPECIATIONNAMES, SOLIDconcs, SOLIDnames]=...
    writerunPHREEQC(Temp,pH,pe,totalnames,totalvector,minerals,speciesexport,database,show,acid,pHfixedcheck);
    SPECIESCONCS=SPECIESCONCS'; SOLIDconcs=SOLIDconcs';
       
    for k=1:size(SPECIESCONCS,1)
          txt=[cell2mat(SPECIATIONNAMES(k,:)),'=SPECIESCONCS(k,:);'];
          eval(txt)
    end

    for k=1:size(SOLIDconcs',1)
          %txt=[cell2mat(SOLIDnames(k,:)),'=SOLIDconcs(k,:)']
          txt=[cell2mat(SOLIDnames(:,k)),'=SOLIDconcs(k,:);'];
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
    Tcalc=masssoln+masssolid; Tcalc(3)=[]; T(3)=[]; % chloride not meaningful because used for fixed pH
    MASSERR=max(100*((T-Tcalc)./T));% relative does not make sense when some components are so tiny (like HFOinitial conc)
    masserr=(T-Tcalc); [MASSERR,INDX]=max(abs(masserr)); MASSERR=MASSERR*sign(masserr(INDX));
end

if HFOsT>0; HFOP=HfowPO4m4+HfowHPO4m3+HfowH2PO4m2+HfosPO4m4+HfosHPO4m3+HfosH2PO4m2; end

end % end of function

% ---------------- SUBFUNCTIONS --------------------------------------------------------
