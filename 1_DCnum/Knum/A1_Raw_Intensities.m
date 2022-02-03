%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       RAW INTENSITIES FOR  KERNEL COMPUTATION - At all points 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tic
clear all; close all;

%% ---------
%Get info from simulations   
rec(:,:,1)=load(['/cluster/scratch/javierd/LargeMedium/Sim1/M1/OUTPUT_FILES/S',num2str(1,'%04.0f'),'.AA.div.semc']);  
%rec(:,:,1)=load(['/Users/alejandro/KernelsComparison_paper/Structural_Change/Decorrelation/Large_Medium_33p6kmsqr/DCnum/Sim1/M1/S',num2str(1,'%04.0f'),'.AA.div.semc']); 
time=rec(:,1,1); %time vector (seconds)
dt=abs( rec(1,1,1)-rec(2,1,1) ); %dt from the simulation
 
 %---Ctes---
L=16800*2; %m/s


%----Lame Parameters-----
vp=6500; %m/s
vs=vp/sqrt(3); %m/s
rho=3750;
mu=rho*(vs)^2;
lambda=rho*(vp^2-2*vs^2);

%lame_ratio=mu^2 / (lambda+2*mu)^2;
 Cte_p=lambda+2*mu;
  Cte_s=mu;

%%
%%{
folder_save=(['/cluster/scratch/javierd/LargeMedium']); 
%folder_save=(['/Users/alejandro/KernelsComparison_paper/Structural_Change/Decorrelation/Large_Medium_33p6kmsqr/DCnum']);


%Chose Model
for M=1:2 % 10  
    %M
   
    
% % FOLDERS - 
 folder_sim1=(['/cluster/scratch/javierd/LargeMedium/Sim1/M',num2str(M,'%1.0f'),'/OUTPUT_FILES/']);   %Sim1 used Exp source
 folder_sim2=(['/cluster/scratch/javierd/LargeMedium/Sim2/M',num2str(M,'%1.0f'),'/OUTPUT_FILES/']);  %Sim2 used elastic force fx ->
 
 %folder_sim1=([folder_save,'/Sim1/M',num2str(M,'%1.0f'),'/']);   %Sim1 used Exp source
 %folder_sim2=([folder_save,'/Sim2/M',num2str(M,'%1.0f'),'/']);   %Sim2 used elastic force fx ->
 

  
%%
%====================SIMULATION I==============================
%================SOURCE AT LOCATION "S"========================  
%data recorded at r is seismogram# 1.  CTE
rcv=1;
     
   
%     %------Total Intensity from s to r -> u^2= ux^2+uz^2--------
   ux=load([folder_sim1,'S',num2str(rcv,'%04.0f'),'.AA.BHX.semc']); 
                
              %---Intensities by Components----
              Ix=abs( hilbert(ux(:,2)) ).^2;
           
 
% Declarations
nr=10; %4903; %Total Number of rcv in OUTPUT_FILES
Ip_save=ones(nr,length(time)); Gp_save=Ip_save;
Is_save=Ip_save;  Gs_save=Ip_save;



%% Loop over space (computing Kernel at several rprime positions)---
for i=1:nr

% Data recorded at r' is seismogram: arbitraty# wherever you want the
% kernel value (at all station locations)

    rprime_div=load([folder_sim1,'S',num2str(i,'%04.0f'),'.AA.div.semc']);
     
    rprime_curl=load([folder_sim1,'S',num2str(i,'%04.0f'),'.AA.cur.semc']); 
      
    
%====================SIMULATION II ==============================
%================SOURCE AT LOCATION "r"========================
%data recorded at rprime is seismogram: arbitraty# wherever you want the kernel value

    rprime_div_simII=load([folder_sim2,'S',num2str(i,'%04.0f'),'.AA.div.semc']); 
          
    rprime_curl_simII=load([folder_sim2,'S',num2str(i,'%04.0f'),'.AA.cur.semc']); 
        
   

%------P-w Intensities-------
   Ip_save(i,:)=abs( hilbert(rprime_div(:,2)) ).^2;
   Gp_save(i,:)=abs( hilbert(rprime_div_simII(:,2) )).^2;
 

%------S-w Intensities-------
  Is_save(i,:)=abs( hilbert(rprime_curl(:,2)) ).^2;
  Gs_save(i,:)=abs( hilbert(rprime_curl_simII(:,2) )).^2;


   
end  % Loop over rcv



%% =======================
Is_save=Is_save';
Ip_save=Ip_save';
Gp_save=Gp_save';
Gs_save=Gs_save';

 indicator_RawInt=(['already run Raw Intensities_M',num2str(M,'%1.0f'),])
      

     
      %=====Save Raw Intensities=====      
 save([folder_save,'/Raw_Intensities_M',num2str(M,'%1.0f'),'.mat'],'Ip_save','Gp_save','Is_save','Gs_save','Ix');%'IP','IS','ux','uz');  %
 
 
end

%toc
