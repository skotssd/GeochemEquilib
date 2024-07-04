function [Fe,Fe(OH)3,MASSERR]=Fetableau(pH,pe,T',flag1,flag2,flag3,flag4,flag5,database)

% input tableau.  change this part % ----------------------------------------------

Tableau=[...
%H      e        Fe          logK                         phase    species1 
1       0        0           0                            0        {'H'}
0       1        0           0                            0        {'e'}
0       0        1           0                            0        {'Fe3'}
-1      0        0           -14                          0        {'OH'}
-2      0        1           -5.67                        0        {'Fe(OH)2+'}
-3      0        1           -12                          0        {'Fe(OH)3'}
-4      0        1           -21.6                        0        {'Fe(OH)4-'}
-4      0        1           -46                          0        {'Fe(OH)4-2'}
-2      0        2           -2.95                        0        {'Fe2(OH)2+4'}
-4      0        3           -6.3                         0        {'Fe3(OH)4+5'}
-1      0        1           -9.5                         0        {'FeOH+'}
-1      0        1           -2.19                        0        {'FeOH+2'}

%solid
-3      0        1           5.6556                       1        {'Fe(OH)3'}
];

% end of tableau.  ------------------ % ----------------------------------------------

[KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES]=processtableau(Tableau,pH,pe);

[SPECIESCONCS,SPECIATIONNAMES,MASSERR,X]=returnspeciationRE(KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES,T,flag1,flag2,flag3,flag4,flag5);

for k=1:size(SPECIESCONCS,1)
      txt=[SPECIATIONNAMES(k,:),'=SPECIESCONCS(k);'];
      eval(txt)
end

end

% ---------------- SUBFUNCTIONS --------------------------------------------------------

