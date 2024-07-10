function [Fe,FeOH3,MASSERR]=Fetableau(pH,pe,T,flag1,flag2,flag3,flag4,flag5)

% input tableau.  change this part % ----------------------------------------------

Tableau=[...
%H      e        Fe          logK                         phase    species1 
1       0        0           0                            0        {'H'}
0       1        0           0                            0        {'e'}
0       0        1           0                            0        {'Fe'}
-1      0        0           -14                          0        {'OH'}
-2      0        1           -5.67                        0        {'FeOH2'}
-3      0        1           -12                          0        {'FeOH3'}
-4      0        1           -21.6                        0        {'FeOH4'}
-2      0        2           -2.95                        0        {'Fe2OH2'}
-4      0        3           -6.3                         0        {'Fe3OH4'}
-1      0        1           -2.19                        0        {'FeOH'}

%solid
-3      0        1           5.6556                       1        {'FeOH3'}
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

