% rough model for portlandite dissolution and calcium carbonate ppte
% pH as a function of dissolved Ca2+ open system

% equilibrium constants
Kw=10^-14; Ka1=10^-6.3; Ka2=10^-10.3; PCO2=10^-3.5; KH=10^-1.47;  

%vector of Ca concs (to be replaced by Ca(OH)2 dissolution versus time
logCaT=-6:0.1:-1; CaT=10.^logCaT; NaT=0; % initial alkalinity by adding sodium to charge bal.

for i=1:length(CaT)
    a=1;
    b=2*CaT(i)+NaT;
    c=-KH*Ka1*PCO2-Kw;
    d=-2*KH*Ka1*Ka2*PCO2;
    t=roots([a b c d]);
    t=t(imag(t)==0); %sets any imaginary roots to zero
    % display the positive real root
    t=t(t>0);
    pH(i)=-log10(t);
end

plot(logCaT,pH,'ko')