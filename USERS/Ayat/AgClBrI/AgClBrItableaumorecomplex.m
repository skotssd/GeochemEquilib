function [Ag,Br,Cl,I,AgBrs,AgCls,AgIs,MASSERR]=AgClBrItableaumorecomplex(pH,pe,TOTALS,flag1,flag2,flag3,flag4,flag5)
%function [Ag]=AgClBrtableaumorecomplex(pH,pe,T,flag1,flag2,flag3,flag4,flag5)
%function [Ag]=AgClBrtableaumorecomplex(pH,pe,TOTALS)

% input tableau.  change this part % ----------------------------------------------

Tableau=[...
%H      e        Ag        Br      Cl      I         logK                             phase    species1 
1       0        0         0       0       0         0                                0       {'H'}
0       1        0         0       0       0         0                                0       {'e'}
0       0        1         0       0       0         0                                0       {'Ag'}
0       0        0         1       0       0         0                                0       {'Br'}
0       0        0         0       1       0         0                                0       {'Cl'}
0       0        0         0       0       1         0                                0       {'I'}
-1      0        0         0       0       1         -14                              0       {'OH'}
0       0        1         1       0       0         4.6                              0       {'AgBr'}
0       0        1         2       0       0         7.5                              0       {'AgBr2'}
0       0        1         3       0       0         8.1                              0       {'AgBr3'}
0       0        1         4       0       0         8.7                              0       {'AgBr4'}
0       0        1         0       1       0         3.2971                           0       {'AgCl'}
0       0        1         0       2       0         5.2989                           0       {'AgCl2'}
0       0        1         0       3       0         5.1310                           0       {'AgCl3'}
0       0        1         0       4       0         3.8050                           0       {'AgCl4'}
0       0        1         0       0       1         6.6                              0       {'AgI'}
0       0        1         0       0       2         11.7                             0       {'AgI2'}
0       0        1         0       0       3         12.6                             0       {'AgI3'}
0       0        1         0       0       4         14.2                             0       {'AgI4'}
%solid
0       0        1         1       0       0         12.30                            1       {'AgBrs'}
0       0        1         0       1       0         9.7453                           1       {'AgCls'}
0       0        1         0       0       1         16.08                            1       {'AgIs'}                                     
];

% end of tableau.  ------------------ % ----------------------------------------------

[KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES]=processtableau(Tableau,pH,pe);

[SPECIESCONCS,SPECIATIONNAMES,MASSERR,X]=returnspeciationRE(KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES,TOTALS,flag1,flag2,flag3,flag4,flag5);

for k=1:size(SPECIESCONCS,1)
      txt=[SPECIATIONNAMES(k,:),'=SPECIESCONCS(k);'];
      eval(txt)
end

end

% ---------------- SUBFUNCTIONS --------------------------------------------------------

