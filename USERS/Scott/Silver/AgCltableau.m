function [Ag,Cl,AgCls,MASSERR]=AgCltableau(pH,pe,T,flag1,flag2,flag3,flag4,flag5)

% input tableau.  change this part % ----------------------------------------------

Tableau=[...
%H      e        Ag        Cl     logK                                phase    species1 
1       0        0         0      0                                   0    {'H'}
0       1        0         0      0                                   0    {'e'}
0       0        1         0      0                                   0    {'Ag'}
0       0        0         1      0                                   0    {'Cl'}
-1      0        0         0      -14                                 0    {'OH'}
%solids
0       0        1         1      9.7453                              1    {'AgCls'} %9.74634 is used by PHREEQC
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

