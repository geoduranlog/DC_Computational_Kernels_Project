%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    AVG INTENSITIES AND  KERNEL COMPUTATION - ALSO IN SPACE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all


%Get info from simulations 
rec(:,:,1)=load(['/cluster/scratch/javierd/LargeMedium/Sim1/M1/OUTPUT_FILES/S',num2str(1,'%04.0f'),'.AA.div.semc']);  
%rec(:,:,1)=load(['/Users/alejandro/KernelsComparison_paper/Structural_Change/Decorrelation/Large_Medium_33p6kmsqr/DCnum/Sim1/M1/S',num2str(1,'%04.0f'),'.AA.div.semc']); 
 time=rec(:,1); %time vector (seconds)
 dt=abs( rec(1,1)-rec(2,1) ); %dt from the simulation
 
 nr=10; %4903; %%Total nb of rcv in OUTPUT_FILES
 nt=length(time);   %total time steps (smooth_time)
 nm=20;  %-----total models to average -----
 
 
 %----Lame Parameters-----
vp=6500; %m/s
vs=vp/sqrt(3); %m/s
rho=3750;
mu=rho*(vs)^2;
lambda=rho*(vp^2-2*vs^2);

%lame_ratio=mu^2 / (lambda+2*mu)^2;
 Cte_p=lambda+2*mu;
  Cte_s=mu;
  
   
%AD
%{  
 %--------
 %%            Creating Mean Intensities -- Run only once
 
 %Raw Intensities
IP=ones(nt,nr,nm); IS=IP; GP=IP;  GS=IP;  %Isr_P=ones(nt,nm); Isr_S=Isr_P; 
Isr_X=ones(nt,nm);

%--Load Kernels in space (no norm)- All models--
for M=1:nm

%--Case: Raw Intensities ---
load(['/cluster/scratch/javierd/LargeMedium/Raw_Intensities_M',num2str(M,'%1.0f')]); 
%load(['/Users/alejandro/KernelsComparison_paper/Structural_Change/Decorrelation/Large_Medium_33p6kmsqr/DCnum/Raw_Intensities_M',num2str(M,'%1.0f')]); 
 
IP(:,:,M)=Ip_save;  
GP(:,:,M)=Gp_save;
IS(:,:,M)=Is_save;    
GS(:,:,M)=Gs_save;  

Isr_X(:,M)=Ix;     


end

clear Ix 





%=======Avg over models=========
IPmean=mean(IP,3);
GPmean=mean(GP,3);
ISmean=mean(IS,3);
GSmean=mean(GS,3);

Isr_Xmean=mean(Isr_X,2);  %Denominator of Kernel


 %==Save Imean==
 save(['Mean_Intensities_',num2str(nm,'%1.0f'),'Models'],'IPmean','GPmean','ISmean','GSmean','Isr_Xmean');  %'Isr_Pmean','Isr_Smean'

%}
 
 %%         Computing the kernel in space
        

load(['Mean_Intensities_',num2str(nm,'%1.0f'),'Models']); 
%load([sigma,'/Mean_SmoothIntensities_',num2str(nm,'%1.0f'),'Models']); 

clear IP GP  IS GS  Isr_X  Ip Is% To free some memory...



  %CONVOLUTION
for i=1:nr
Nump(:,i)=conv(Cte_p.*IPmean(:,i),Cte_p.*GPmean(:,i)).*dt;  
Nums(:,i)=conv(Cte_s.*ISmean(:,i),Cte_s.*GSmean(:,i)).*dt ;

Numps(:,i)=conv(Cte_p.*IPmean(:,i),Cte_s.*GSmean(:,i)).*dt ;
Numsp(:,i)=conv(Cte_s.*ISmean(:,i),Cte_p.*GPmean(:,i)).*dt ;
end


%%
%time_smooth=time;  % Only for Raw intensities (Comment this part in case
                    %of Smooth Intensities)


%Convention->you recorded Ux with an Geophone at R (I add the respectively constant in the denominator)
cte_x=2*rho*(20^2)*pi^2;    % To get Ix: Energy density due to a displacement ux


% Divide each row of Numerator by denominator
Kp=bsxfun(@rdivide,Nump(1:length(time),:),cte_x.*Isr_Xmean);
Ks=bsxfun(@rdivide,Nums(1:length(time),:),cte_x.*Isr_Xmean);

Kps =bsxfun(@rdivide, Numps(1:length(time),:),cte_x.*Isr_Xmean);
Ksp =bsxfun(@rdivide, Numsp(1:length(time),:),cte_x.*Isr_Xmean);




%% =====SAVE KERNEL====
%save(['Kernels_',num2str(nm,'%1.0f'),'Models'],'Kp','Ks','Kps','Ksp');  

%=======LOAD KERNEL======
load([sigma,'/Kernels_',num2str(nm,'%1.0f'),'Models'],'Kp','Ks');  %

rp=3
figure
plot(Kp(:,rp))
hold on
plot(Kps(:,rp))
hold on
plot(Ks(:,rp))
legend('Kp','Kps','Ks')
%xlim([])


figure
plot(ISmean(:,rp)./IPmean(:,rp))

