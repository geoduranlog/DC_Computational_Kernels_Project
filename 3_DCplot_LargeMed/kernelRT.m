function [KRad]= kernelRT(time,dt,dist_s_r, dist_s_rprime, dist_r_rprime,c, l)


%% Kernel BASED ON RADIATIVE TRANSFER
t=time;


%--Intensities------
eps=1e-200; %Epsilon, in case R=0, make R=eps

R=dist_s_rprime;     if R==0  R=eps; end    
Isrprime=(1./(2*pi*R)).*dirac(c*t-R).*exp(-c*t./l) + (1./(2*pi*l*c*t)).*( (1- R*R./(c*c*t.^2) ).^(-0.5) ) .* exp( (sqrt( c*c*t.^2 - R*R) -c.*t )/l )  .* heaviside(c.*t - R);

R=dist_r_rprime;    if R==0  R=eps; end 
Irrprime=(1./(2*pi*R)).*dirac(c*t-R).*exp(-c*t./l) + (1./(2*pi*l*c*t)).*( (1- R*R./(c*c*t.^2) ).^(-0.5) ) .* exp( (sqrt( c*c*t.^2 - R*R) -c.*t )/l )  .* heaviside(c.*t - R); % In Ktheo they asume same source type (as far as I know) 


R=dist_s_r;     if R==0  R=eps; end 
Itotal=(1./(2*pi*R)).*dirac(c*t-R).*exp(-c*t./l) + (1./(2*pi*l*c*t)).*( (1- R*R./(c*c*t.^2) ).^(-0.5) ) .* exp( (sqrt( c*c*t.^2 - R*R) -c.*t )/l )  .* heaviside(c.*t - R);
Itotal_RT=Itotal;

  
%   %--Kernel---
NumRad=conv(Isrprime,Irrprime).*dt;  %conv(real(Isrprime),real(Irrprime)).*dt;
KRad=NumRad(1:length(t))./real(Itotal); %I think this should be, but Num=0 in this interval [0 time(end)]


end