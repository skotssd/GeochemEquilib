% increase Ca versus time by dissolution of portlandite
% and model pH increase

function II=portlanditedissolution

Ksp=5.5e-6; k=4e-2; Cainf=Ksp^(1/3); p=[Cainf k];
CT=0.01; NaT=0.005;

time=0:1:100;
CaT=CaTvtime(p,time);
figure(1); plot(time,CaT,'ko')

% for each timestep solve for pH
for i=1:length(time)
    %pH(i)=pHfromCaT(CaT(i),NaT);
    %pHlowCO2(i)=pHfromCaTalmostclosed(CaT(i),NaT); % just OH from Ca(OH)2 dissolution
    pH(i)=pHfromCaT(CaT(i),0); % ignore initial alkalinity
    pHlowCO2(i)=pHfromCaTalmostclosed(CaT(i),0); % just OH from Ca(OH)2 dissolution. ignore initial alkalinity
    pHclosed(i)=pHfromCaTclosed(CaT(i),CT,NaT);  % NaT for initial alkalinity. otherwise start acidic (H2CO3)
end

%pH(1) %same as rainwater
%pHlowCO2(1) % neutral

figure(2); plot(time,pH,time,pHlowCO2,time,pHclosed)

%%%%%%%%%%%%%%%% end of main scripts %%%%%%%%%%%%%%%%%%%%%%%%%%%

end

function CaT=CaTvtime(p,t)

Cainf=p(1); k=p(2); CaT(1)=0;

for i=2:length(t)
CaT(i)=Cainf*(1-exp(-k*t(i)));
end

CaT=Cainf*(1-exp(-k*t));

end

function pH=pHfromCaT(CaT,NaT)

% equilibrium constants
Kw=10^-14; Ka1=10^-6.3; Ka2=10^-10.3; PCO2=10^-3.5; KH=10^-1.47;  

a=1;
b=2*CaT+NaT; %if there is initial alkalinity
c=-KH*Ka1*PCO2-Kw;
d=-2*KH*Ka1*Ka2*PCO2;
t=roots([a b c d]);
t=t(imag(t)==0); %sets any imaginary roots to zero
% display the positive real root
t=t(t>0);
pH=-log10(t);


end

function pH=pHfromCaTalmostclosed(CaT,NaT)

% equilibrium constants
Kw=10^-14; Ka1=10^-6.3; Ka2=10^-10.3; PCO2=10^-26.5; KH=10^-1.47;  

a=1;
b=2*CaT+NaT; %if there is initial alkalinity
c=-KH*Ka1*PCO2-Kw;
d=-2*KH*Ka1*Ka2*PCO2;
t=roots([a b c d]);
t=t(imag(t)==0); %sets any imaginary roots to zero
% display the positive real root
t=t(t>0);
pH=-log10(t);


end

function pH=pHfromCaTclosed(CaT,CT,NaT)

% equilibrium constants
Kw=10^-14; Ka1=10^-6.3; Ka2=10^-10.3; PCO2=10^-23.5; KH=10^-1.47;  

a=1;
b=2*CaT+Ka1+NaT;
c=2*CaT*Ka1+NaT*Ka1+Ka1*Ka2-CT*Ka1-Kw;
d=2*CaT*Ka1*Ka2+NaT*Ka1*Ka2-CT*Ka1*Ka2-Kw*Ka1;
e=-Kw*Ka1*Ka2;
t=roots([a b c d e]);


t=t(imag(t)==0); %sets any imaginary roots to zero
% display the positive real root
t=t(t>0);
pH=-log10(t);


end