%--DC plots---
clear all
close all

%Medium Type
medium='Elastic';  % Acoustic  Elastic

%Folder to save figures
folder_save=(['/Users/alejandro/KernelsComparison_paper/Structural_Change/Decorrelation/Large_Medium_33p6kmsqr/DCplot_LargeMed/']); 


 %% =======================SUBPLOTS SEVERAL rprime LOCATIONS===============================
 
 %Load Good DCnum errors
  load([folder_save,'/DCnumErrors_40sets.mat']);   %This, is from 16.8x16.8 km media
  
  %Mean free time (s)  - using effective parameters
  mft=0.2819; 
 
   num=0;
  %%{
 locations=[3 8 9]; %[5 18 20];
 
 for rp=locations;  
num=num+1;
     
  %Load DCs
   folder_DC=(['/Users/alejandro/KernelsComparison_paper/Structural_Change/Decorrelation/Large_Medium_33p6kmsqr/DCplot_LargeMed/DCplot_LargeMed']);
   load([folder_DC,'DC_EL_plot_LargeMed_rp',num2str(rp,'%01.0f'),'.mat'])
      
    t=time;   dt=time(2)-time(1);
  
    %=====================================
  % Time Limit in x_axis for plot 
   % to_mft=5; % in mft
    tf=7.0475; %seconds
    %=====================================
    
    %DCexp time due to timewindow integration
    temp=time(1+twm:length(time)-twm);   %time(1+2*twm:length(time));  %Real time is time(1+twm:length(time)-twm);  
    smth_nb=1e4%6e4% 3e3; 
    nto=1510;   %Choose a starting time after singularities
     
    %--Smooth DCexp=Mean----
     Smooth_DCp_exp=smooth(DC_expP(nto:end),smth_nb);
     Smooth_DCs_exp=smooth(DC_expS(nto:end),smth_nb);
     Smooth_DCBoth_exp=smooth(DC_expBoth(nto:end),smth_nb);
    
    temp=temp(nto:end);  
    
    
         k=time(1+twm:length(time)-twm);
        for j=1:length(temp)       
    [c(j) index(j)] = min(abs(k-temp(j)));
end
    IQR_P2=IQR_P(index);  %Keep errors as they are but place them for the smoothed DCexp points
    IQR_S2=IQR_S(index); 
    IQR_both=IQR_both(index); 
    
    
    %---Reduce points for plotting---
    delta=round(length(Smooth_DCp_exp)/40);
    Smooth_DCp_exp=Smooth_DCp_exp(1:delta:end);
    Smooth_DCs_exp=Smooth_DCs_exp(1:delta:end);
    Smooth_DCBoth_exp=Smooth_DCBoth_exp(1:delta:end);
    
    IQR_P2=IQR_P2(1:delta:end); 
    IQR_S2=IQR_S2(1:delta:end);
    IQR_both=IQR_both(1:delta:end);
   
     temp=temp(1:delta:end); 

    
     %-----------DCnum Smoothed--------------
    % Don't load smooth DCnum. Simply take the raw DCnum and smooth them
    % here:
          smth_nb=1e4; %6e4% 3e3; 
    nto=5300;
     Smooth_DC_num_P=smooth(DC_num_P(nto:end),smth_nb);
     Smooth_DC_num_S=smooth(DC_num_S(nto:end),smth_nb);
     time_smth_EL=time(nto:end);


%Reduce points
    delta=round(length(Smooth_DC_num_P)/40);
    Smooth_DC_num_P=Smooth_DC_num_P(1:delta:end);
    Smooth_DC_num_S=Smooth_DC_num_S(1:delta:end);
    
    Smooth_DC_num_both=Smooth_DC_num_P+Smooth_DC_num_S;  % DCnum=DCnumP+DCnumS
    time_smth_EL=time_smth_EL(1:delta:end); 

    % --- Locatios where DC_EL (num) was computed ----
     r_all=[2 5 12 18 20];    %Easy Model
      nr=find(rp==r_all); %Find proper location

 %%%Load Good DCnum errors -- already loaded ...
  %%%   load(['/Users/alejandro/KernelsComparison_paper/Structural_Change/Decorrelation/DCnumErrors_40sets.mat']);
        
   
    
    %=============================================================  
  %==== DCnum Error at same tsteps as SmoothDCnum ==========
% DCnumEr has more values than Smooth_DCnum (because the error was computed in another code)
%Thus, for the same time put Smooth DC with its proper error
clear c index
for j=1:length(time_smth_EL)  %SmoothDC time      
    [c(j) index(j)] = min(abs(t_sth-time_smth_EL(j)));
end

DCnumP_Er2=DCnumP_Er(index,:);
DCnumS_Er2=DCnumS_Er(index,:);
DCnumBoth_Er2=DCnumP_Er2+DCnumS_Er2;

%Temp - New labels in Large Media
DCnumP_Er2(:,3)=DCnumP_Er2(:,5); DCnumS_Er2(:,3)=DCnumS_Er2(:,5);  DCnumBoth_Er2(:,3)=DCnumBoth_Er2(:,5); 
DCnumP_Er2(:,8)=DCnumP_Er2(:,18); DCnumS_Er2(:,8)=DCnumS_Er2(:,18);  DCnumBoth_Er2(:,8)=DCnumBoth_Er2(:,18); 
DCnumP_Er2(:,9)=DCnumP_Er2(:,8); DCnumS_Er2(:,9)=DCnumS_Er2(:,8);  DCnumBoth_Er2(:,9)=DCnumBoth_Er2(:,8);


%===============================
    
     

    
    
       % ==================BIG FIGURE EASY MODEL================
 
       % Times in terms of mft
       temp=temp/mft;
       time_smth_EL=time_smth_EL/mft;
       time=time/mft;
       
       %% P-WAVE ANOMALY
       rgb=252;
    figure(5)
     subplot(3,3,num)
    
       alpha=0.5
% %     P-waves
H1=     shadedErrorBar( temp, Smooth_DCp_exp*100, IQR_P2*100 ,'ko', 0.9 ); %  % IQR_sth  IQR2
  set(H1.mainLine,'MarkerSize',8)
     set(H1.mainLine,'MarkerFaceColor',[128/rgb,128/rgb,128/rgb]) %Gray Color

     hold on
      H2=errorbar(time_smth_EL,Smooth_DC_num_P*100,DCnumP_Er2(:,rp)*100,'kd','MarkerSize',8);
  %   H2=errorbar(time_smth_EL,Smooth_DC_num_all_EL(:,nr,1)*100,DCnumP_Er2(:,rp)*100,'kd','MarkerSize',8);
  set(H2,'MarkerFaceColor','b') 
   hold on 
 
   legend([H1.mainLine, H2], ...
    'DC_{P exp}','DC_{P num}', ...
    'Location', 'Northwest');
title(['r' num2str(rp),],'fontsize',14);
  set(gca,'fontsize',23)   
   xlabel('t^{*} ') % xlabel('time (s)')
ylabel('DC (%)')
   xlim([1.4095  tf]/mft); 
   
   %--Change axes only for rp=r12--
   if rp==3
   ylim([0 6e-6]); %1e-4] For P-waves 
   else
   ylim([0 6e-6]);   %For P-waves
   end


     
 %% S-WAVE ANOMALY
%%{
     subplot(3,3,num+3) 
 hold on
  H3=     shadedErrorBar( temp, Smooth_DCs_exp*100, IQR_S2*100 ,'ko', 0.9 ); %  % IQR_sth  IQR2
  set(H3.mainLine,'MarkerSize',8)
     set(H3.mainLine,'MarkerFaceColor',[128/rgb,128/rgb,128/rgb]) %Gray Color
     hold on
     H4=errorbar(time_smth_EL,Smooth_DC_num_S*100,DCnumS_Er2(:,rp)*100,'kd','MarkerSize',8);
  set(H4,'MarkerFaceColor','r') 
  legend([H3.mainLine, H4], ...
    'DC_{S exp}','DC_{S num}', ...
    'Location', 'Northwest');
set(gca,'fontsize',23)   
  hold on
   xlabel('t^{*} ') % xlabel('time (s)')
ylabel('DC (%)')

    %--Change axes only for rp=r12--
  xlim([1.4095 tf]/mft); 
  
    if rp==3
   ylim([0 5e-5]); %1e-4] For S-waves 
   else
   ylim([0 2e-5]);   %For S-waves
   end

 %%}
 
  %% BOTH- dVP & dVS TOGETHER
%%{
     subplot(3,3,num+6) 
 hold on
  H3=     shadedErrorBar( temp, Smooth_DCBoth_exp*100, IQR_both*100 ,'ko', 0.9 ); %  % IQR_sth  IQR2
  set(H3.mainLine,'MarkerSize',8)
     set(H3.mainLine,'MarkerFaceColor',[128/rgb,128/rgb,128/rgb]) %Gray Color
     hold on
     H4=errorbar(time_smth_EL,Smooth_DC_num_both*100,DCnumBoth_Er2(:,rp)*100,'kd','MarkerSize',8);
  set(H4,'MarkerFaceColor','g') 
 %    H4=errorbar(time_smth_EL,Smooth_DC_num_both*100,DCnumBoth_Er2(:,rp)*100,'kd','MarkerSize',8);
  %set(H5,'MarkerFaceColor','k--') 
  H5=plot(time,DCbulk*100,'k--','LineWidth',3)
 
  legend([H3.mainLine, H4,H5], ...
    'DC_{exp}','DC_{num}','DC_{Bulk}', ...
    'Location', 'Northwest');
set(gca,'fontsize',23)   
  hold on
   xlabel('t^{*} ') % xlabel('time (s)')
ylabel('DC (%)')

    %--Change axes only for rp=r12--
   xlim([1.4095 tf]/mft); 
    if rp==3
   ylim([0 5e-5]); %1e-4] For S-waves 
   else
   ylim([0 2e-5]);   %For S-waves
   end
 %%}
 
 
  %---Print at the end of the loop----
   if num==length(locations) %   7 is the total number of figures to show
   
         fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 16 13];  %[0 0 15 12];
print(gcf,[folder_save,'DC_all'],'-dpng','-r800');  %dpi set at 300 to increase image resolution 
 
% 
%          fig = gcf;
% %fig.PaperUnits = 'inches';
% fig.PaperPosition =[0 0 15 12];  %[0 0 12 6];
% print(gcf,[folder_save,'DC_allTEST'],'-dpdf','-r800');  %dpi set at 300 to increase image resolution 
%  


   end
    
   
   
   
   
   
    
 
 
 end
 %}
 

