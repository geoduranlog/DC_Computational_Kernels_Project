%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%          DC components for Plotting - 2D ELASTIC MEDIUM   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute all DCs and save for future plotting
clear all
close all


%% --- Theoretical RT Kernel - effective parameters ---

%Large Media
stations=importdata(['/Users/alejandro/KernelsComparison_paper/Structural_Change/Decorrelation/Large_Medium_33p6kmsqr/DCnum/Sim1/M1/STATIONS']); 

%To acces the variables inside the struct: structName(indices).fieldName
stations=stations(1,1).data;   % I take only numerical values;
stations=stations(:, [1 2]);    %Take only the (x,y) values
xr=stations(:,1); zr=stations(:,2);


    % Select rcv - which corresponds to the offset used in the theoretical
     % computation of the Intensity - OJO!
     rcv=1;

     % Source is at fixed position (x,z)=(8400.24,11005.87) m   
     xs=16800.48; zs=21506.17 ; %Src of Sim1 Large Media
     
     
     %-------------------------------------------------
     %Define rprime position - where Kernel is computed 
     rp=3 %3 8 9      %r_prime is at this rcv#
     xr_prime=xr(rp);   zr_prime=zr(rp);
     %-------------------------------------------------

%PARAMETERS-- Elastic Medium--
vp=6500; %m/s
vs=vp/sqrt(3);
c=(vs*vp)/(0.75*vp+0.25*vs); %Energy Vel  0.75/0.25=Is/Ip=3 - Anne used 0.77/0.23=3.35 she got a diff Equipart due to BC and free surface 
L=1183; %940AC  %Transport mean free path (m)  ->Diff
l=1189; %900AC  %Scattering mean free path (m) ->RT
nt=60000; %Nstep seismograms Large media simulaitons
 
Norm=2.97e9; %EL force;   %Norm=3.5377e+05; %Explosion at R Sim2

mft=L/c;  %Mean free time


 %----Lame Parameters-----
rho=3750;
mu=rho*(vs)^2;
lambda=rho*(vp^2-2*vs^2);
 Cte_p=lambda+2*mu;
  Cte_s=mu;


%Distances
dist_s_r=sqrt( (xs-xr(rcv)).^2 + (zs-zr(rcv)).^2 );   %rcv-src distance (m)
dist_s_rprime=sqrt( (xs-xr_prime).^2 + (zs-zr_prime).^2 );   %rprime-src distance (m) Only for Kernel computation
dist_r_rprime=sqrt( (xr(rcv)-xr_prime).^2 + (zr(rcv)-zr_prime).^2 );   %r-rprime distance (m) Only for Kernel computation


D=c*L/2; %Diffusion coeff in 2D
dt=3e-4;
time=0.0001:dt:0.0001+(nt-1)*dt;  %time vector (s)  -> Make t>0 to avoid complex numbers in the computation
                                    %Satisfies that length(time)=nt

time_seismo=-0.06:dt:20.9397;  %nt=70mil  %17.9397; %20.9397; %This is exactly how the seismo time is.


clear Isrprime
clear Irrprime
clear Itotal

%===============================
%--Function  Theo Kernel RT ----
 [KRad]=kernelRT(time,dt,dist_s_r, dist_s_rprime, dist_r_rprime,c, l);
 
 
 %% ============= DECORRELATION - BULK ==============  

%Scatt Cross Section -Symetric GLL points %->data for dvp/vp~7%   
   sigmaPP=2.79e-4;    
   sigmaPS=6e-9;  
   sigmaSS=4.8789e-4; 
   sigmaSP=3.1864e-5;   %Margerin: sgimaPS=2*( (vp^2)/(vs^2))*sigmaSP Hyp:that's in 3D.  Then Hyp in 2D?:      (vp/vs) *sigmaSP

   
%======================
%--DC Bulk--
%=======================
%OJO:Cte used during DCnum computation 
cte_eff=0.5*vp*(sigmaPP+sigmaPS)+ 0.5* vs*(sigmaSS+sigmaSP) 

N=1; % Nb scatterers in dV'
DCbulk=N*cte_eff.*KRad; 


    
 %% ============= DCexp and DCnum ==============  
 
   %Case of Large Media
   
    %DCexp
   load(['/Users/alejandro/KernelsComparison_paper/Structural_Change/Decorrelation/Large_Medium_33p6kmsqr/DC_exp/Both_vel/dv_7pro/DCexp_BothVel_ux_Elastic_LargeMed_rp',num2str(rp,'%01.0f'),'_twm1666.mat']) 
    load(['/Users/alejandro/KernelsComparison_paper/Structural_Change/Decorrelation/Large_Medium_33p6kmsqr/DC_exp/dvp_only/dv_7pro/DCexp_P_ux_Elastic_LargeMed_rp',num2str(rp,'%01.0f'),'_twm1666.mat']) 
     load(['/Users/alejandro/KernelsComparison_paper/Structural_Change/Decorrelation/Large_Medium_33p6kmsqr/DC_exp/dvs_only/dv_7pro/DCexp_S_ux_Elastic_LargeMed_rp',num2str(rp,'%01.0f'),'_twm1666.mat']) 
    
   twm=1666;
     
     %DCnum
   load(['/Users/alejandro/KernelsComparison_paper/Structural_Change/Decorrelation/Large_Medium_33p6kmsqr/DCnum/Knum/Kernels_20Models.mat']) 
   
   
   
  DC_num_P=0.5*vp*(sigmaPP.*Kp(:,rp)./Norm+sigmaPS.*Kps(:,rp)./Norm);  
  DC_num_S=0.5* vs*(sigmaSS.*Ks(:,rp)./Norm+sigmaSP.*Ksp(:,rp)./Norm);
   
   %DCexp time - based on half time Window twm 
   time_exp=time_seismo(1+twm:length(time)-twm); 
   
   % Considering both contributions
   DC_num=DC_num_P+DC_num_S;
  
   
   
   % Equipartition Ratio
   load(['/Users/alejandro/KernelsComparison_paper/Structural_Change/Decorrelation/Large_Medium_33p6kmsqr/DCnum/Knum/Mean_Intensities_20Models.mat']) 
        
     
   
   %% DCexp Error
nmodels=20;
skip_models=[0] ; %[2 17] ->for rp3;


  %=========== Vs and Vp Placed Together======
   m=nmodels-length(skip_models);
        % median
    DCexp_both_Me=median(DC_both(:,1:m),2);
    
    %{
    %Avg but skipping outlier models 
    DC_expBoth_skip=zeros(length(DC_both),m);
    
    count=0;
    for i=1:nmodels
        
        if ismember(i,skip_models)==0
            count=count+1;
            
        DC_expBoth_skip(:,count)=DC_both(:,i);
        end
    
    end
    %}
    
    %==================================================
    %====Error DCexp - IQR  (Interquartile range)======
    IQR_both=iqr(DC_both(:,1:m),2); 
    %IQR_both=iqr(DC_expBoth_skip(:,1:m),2); 
    
    IQR_P=iqr(DC_P(:,1:m),2); 
    IQR_S=iqr(DC_S(:,1:m),2); 
    
    %==================================================
     
    if rp==8
 DC_expS=DC_expBoth; % Temporary - There is an issue  with the simulation (for rp8 dvs only). Re make...
   Ss=S_both;
   SE_s=SE_both;
end

   
   %% FIGUREs
   
%    % Equipartition Ratio Kernels
%    figure(200)
%    plot(time./mft,Ks(:,rp)./Kp(:,rp),'k')
%   %xlim([1.5 15]./mft); 
%   %xlim([9 40]); 
%    xlabel('t/t^{*} ')
%    ylabel('K_{SS}/K_{PP}')
%    set(gca,'fontsize',28)
%    ylim([0 10]); 
%       xlim([0 60]); 
   
  % Equipartition Ratio - Intensities
   figure(201)
    plot(time./mft,(Cte_s.*ISmean(:,rp))./(Cte_p.*IPmean(:,rp)),'k')
    hold on
   xlabel('t/t^{*} ')
   ylabel('I_{S}/I_{P}')
   set(gca,'fontsize',28)
    %  ylim([0 10]); 
   
   
   
     figure(3)
   plot(time,100*DC_num,'g')
     hold on
    plot(time_exp,100.*DC_expBoth,'LineWidth',3)  %./time_exp'
   hold on
   plot(time,100*DCbulk,'k--','LineWidth',5)
   legend({'DCnum_P+DCnum_S','DCexpBoth','DC_{Bulk}'},'Location','northWest','FontSize',16)
  %legend({'DCexpBoth','DC_{Bulk}'},'Location','northWest','FontSize',16)
   title(['DCs - rp=',num2str(rp)],'FontSize',18)
   xlim([1.5 15]);
   xlim([1.5 5]); 
   xlabel('time (s)')
   ylabel('DC (%)')
   set(gca,'fontsize',28)
   

   %close all
   figure(5)
   plot(time_exp,100*DC_both(:,1:20))
     hold on
    % %     %DCexp_BothVel
 plot(time_exp,DC_expBoth*100,'ko','LineWidth',4)  %plot(time(1+2*twm:length(time)),DC_expBoth*100,'k','LineWidth',4) 
xlim([1.5 15]);
title(['DC_{exp} Both - All models - rp=',num2str(rp)],'FontSize',18)
   



% ====================================================================
%% ===================  SAVE COMPONENTS FOR PLOTTING===================

folder_DC=(['/Users/alejandro/KernelsComparison_paper/Structural_Change/Decorrelation/Large_Medium_33p6kmsqr/DCplot_LargeMed/DCplot_LargeMed']);
%save([folder_DC,'DC_EL_plot_LargeMed_rp',num2str(rp,'%01.0f'),'.mat'],'time','time_exp','DCbulk','DC_num_P','DC_num_S','DC_expP','DC_expS','DC_expBoth','Sp','Ss','S_both','SE_p','SE_s','SE_both','IQR_P','IQR_S','IQR_both','DC_num_S','twm');  
   
% =================== =========================================



%%   TEST PLOT

     t=time;   dt=time(2)-time(1);
    tf=15; %seconds -> Final time for plots
         
     %DCexp time due to timewindow integration
    temp=time(1+twm:length(time)-twm);   %time(1+2*twm:length(time));  %Real time is time(1+twm:length(time)-twm);
   
    %Smooth DCexp
   smth_nb=1e4
nto=1510; %1300;  %Choose a starting time after singularities
% %--Smooth DCexp=Median----
%     Smooth_DCp_expMe=smooth(DCexpP_Me(nto:end),smth_nb);
%     IQR_P_sth=smooth(IQR_P(nto:end),smth_nb);
%      Smooth_DCs_expMe=smooth(DCexpS_Me(nto:end),smth_nb);
%     IQR_S_sth=smooth(IQR_S(nto:end),smth_nb);
%     

    %--Smooth DCexp=Mean----
     Smooth_DCp_exp=smooth(DC_expP(nto:end),smth_nb);  
   % IQR_P_sth=smooth(IQR_P(nto:end),smth_nb);
     Smooth_DCs_exp=smooth(DC_expS(nto:end),smth_nb);  
    %IQR_S_sth=smooth(IQR_S(nto:end),smth_nb);
    Smooth_DCBoth_exp=smooth(DC_expBoth(nto:end),smth_nb);  
    
    temp=temp(nto:end);  
    
    
         k=time(1+twm:length(time)-twm);
    
    for j=1:length(temp)       
    [c(j) index(j)] = min(abs(k-temp(j)));
end
    IQR_P2=IQR_P(index);  %Keep errors as they are but place them for the smoothed DCexp points
    IQR_S2=IQR_S(index); 
    
    
    %---Reduce points for plotting---
    delta=round(length(Smooth_DCp_exp)/40);
   % Smooth_DCp_expMe=Smooth_DCp_expMe(1:delta:end);
   % Smooth_DCs_expMe=Smooth_DCs_expMe(1:delta:end);
    
     Smooth_DCBoth_exp=Smooth_DCBoth_exp(1:delta:end);
    Smooth_DCp_exp=Smooth_DCp_exp(1:delta:end);
    Smooth_DCs_exp=Smooth_DCs_exp(1:delta:end);
    
   %IQR_P_sth=IQR_P_sth(1:delta:end);  %better don't use smoothed errors
   % IQR_S_sth=IQR_S_sth(1:delta:end);
    
    IQR_P2=IQR_P2(1:delta:end); 
    IQR_S2=IQR_S2(1:delta:end);
   
     temp=temp(1:delta:end); 

    
    
%-----------DCnum Smoothed--------------
          smth_nb=1e4; %6e4% 3e3; 
    nto=5300;
     Smooth_DC_num_P=smooth(DC_num_P(nto:end),smth_nb);
     Smooth_DC_num_S=smooth(DC_num_S(nto:end),smth_nb);
     time_smth_EL=time(nto:end);


%Reduce points
    delta=round(length(Smooth_DC_num_P)/40);
    Smooth_DC_num_P=Smooth_DC_num_P(1:delta:end);
    Smooth_DC_num_S=Smooth_DC_num_S(1:delta:end);
    time_smth_EL=time_smth_EL(1:delta:end); 


    
      %--Elastic DC--
      rgb=252;
     figure(10)
    alpha=0.5
% %     P-waves
%%{
H1=     shadedErrorBar( temp, Smooth_DCp_exp*100, IQR_P2*100 ,'ko', 0.9 ); %  % IQR_sth  IQR2
  set(H1.mainLine,'MarkerSize',8)
     set(H1.mainLine,'MarkerFaceColor',[128/rgb,128/rgb,128/rgb]) %Gray Color

     hold on
     H2=plot(time_smth_EL,Smooth_DC_num_P*100,'kd','MarkerSize',8);
     set(H2,'MarkerFaceColor','b') 
   hold on
 %} 
  % % S-waves
  H3=     shadedErrorBar( temp, Smooth_DCs_exp*100, IQR_S2*100 ,'ko', 0.9 ); %  % IQR_sth  IQR2
  set(H3.mainLine,'MarkerSize',8)
     set(H3.mainLine,'MarkerFaceColor',[128/rgb,128/rgb,128/rgb]) %Gray Color
     hold on
     H4=plot(time_smth_EL,Smooth_DC_num_S*100,'kd','MarkerSize',8);
  set(H4,'MarkerFaceColor','r') 
  hold on
  %plot(t(5035:length(time)),DC_num_S(5035:end)*100,'r') %Raw DCnum
     grid on
  title(['Elastic Decorrelations - r' num2str(rp),],'fontsize',14);

 legend([H1.mainLine, H2, H3.mainLine, H4], ...
    'DC_{P exp}','DC_{P num}','DC_{S exp} ','DC_{S num}', ...
    'Location', 'Northwest');
  set(gca,'fontsize',23)   
    xlabel('time (s)')
ylabel('DC (%)')
  xlim([1.4 tf]); 
 % ylim([0 2e-5]); %For P-waves 

 