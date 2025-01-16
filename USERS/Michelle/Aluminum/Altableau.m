function [Al,AlOH3s,Al13,AlOH4,MASSERR]=Altableau(pH,pe,T,flag1,flag2,flag3,flag4,flag5)

% input tableau.  change this part % ----------------------------------------------

Tableau=[...
%H      e        Al             logK                                phase    species1 
1       0        0               0                                   0    {'H'}
0       1        0               0                                   0    {'e'}
0       0        1               0                                   0    {'Al'}
-1      0        0               -13.9951                            0    {'OH'}
-2      0        1               -10.5945                            0    {'AlOH2'}
-1      0        1               -4.9571                             0    {'AlOH'}
-1      0        1               -4.9571                             0    {'AlOH'}
-4      0        1               -22.7                               0    {'AlOH4'}
-2      0        2               -7.7902                             0    {'Al2OH2'}
-4      0        3               -13.8803                            0    {'Al3OH4'}
-32     0        13              -98.73                              0    {'Al13'}
%solids
-3       0       1               -7.7560                             1    {'AlOH3s'}
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
