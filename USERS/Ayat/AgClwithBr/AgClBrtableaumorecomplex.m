function [Ag,Br,Cl,AgBrs,AgCls,MASSERR]=AgClBrtableaumorecomplex(pH,pe,TOTALS,flag1,flag2,flag3,flag4,flag5)
%function [Ag]=AgClBrtableaumorecomplex(pH,pe,T,flag1,flag2,flag3,flag4,flag5)
%function [Ag]=AgClBrtableaumorecomplex(pH,pe,TOTALS)

% input tableau.  change this part % ----------------------------------------------

Tableau=[...
%H      e        Ag        Br      Cl       logK                             phase    species1 
1       0        0         0       0        0                                0       {'H'}
0       1        0         0       0        0                                0       {'e'}
0       0        1         0       0        0                                0       {'Ag'}
0       0        0         1       0        0                                0       {'Br'}
0       0        0         0       1        0                                0       {'Cl'}
-1      0        0         0       0        -14                              0       {'OH'}
0       0        1         1       0        4.6                              0       {'AgBr'}
0       0        1         2       0        7.5                              0       {'AgBr2'}
0       0        1         3       0        8.1                              0       {'AgBr3'}
0       0        1         4       0        8.7                              0       {'AgBr4'}
0       0        1         0       1        3.2971                           0       {'AgCl'}
0       0        1         0       2        5.2989                           0       {'AgCl2'}
0       0        1         0       3        5.1310                           0       {'AgCl3'}
0       0        1         0       4        3.8050                           0       {'AgCl4'}
%solid
0       0        1         1       0        12.30                            1       {'AgBrs'}
0       0        1         0       1        9.7453                           1       {'AgCls'}
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

