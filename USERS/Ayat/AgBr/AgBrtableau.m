function [Ag,Br,AgBrs,MASSERR]=AgBrtableau(pH,pe,T,flag1,flag2,flag3,flag4,flag5)

% input tableau.  change this part % ----------------------------------------------

Tableau=[...
%H      e        Ag        Br     logK                                phase    species1 
1       0        0         0      0                                   0    {'H'}
0       1        0         0      0                                   0    {'e'}
0       0        1         0      0                                   0    {'Ag'}
0       0        0         1      0                                   0    {'Br'}
-1      0        0         0      -14                                 0    {'OH'}
%solids
0       0        1         1      12.30                               1    {'AgBrs'}
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

