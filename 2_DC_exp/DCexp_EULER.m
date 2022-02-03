%% =========MEASURED DECORRELATION==========
clear all
close all

medium='Elastic'; %Acoustic   Elastic
independent='no'; %yes no  %Vp and Vs perturbations were placed separated or together as vs=vp/sqrt(3)
        
nt=60000; %60000 %Total NSTEPS

%Half Time window
twm=1666; %208; %1000;  % twm*dt ->in seconds.   The ENTIRE window corresponds to 2*twm*dt/(1/20) periods of the wave I sent
nmodels=20; %60 (acoustic);
dt=3e-4;
rcv=1; %RCV - fixed location.   RCV is 6 if using Anne Layer configuration
rp=3; %Perturbation location   

%skip_models=[];   %[16 20 24 34 35]; %These models didn't run, but I made new ones

%for rp= [  18 20 ]

CC=zeros(nt-2*twm ,nmodels);   %zeros(length(w(:,1))-2*twm ,1);  
for M=1:nmodels

    %--ELASTIC----
    if strcmp(medium,'Elastic')==1
        
        
        if strcmp(independent,'yes')==1 
           % vel=='vp'; %vs   %Choose perturbation type
           
              %Load Non-perturbed signal - i.e., Sim1
w=load(['/cluster/scratch/javierd/Epsilon/Sim1_New/M',num2str(M,'%01.0f'),'/S',num2str(rcv,'%04.0f'),'.AA.BHX.semc']);  
w=w(1:nt,:); %to make w and w_pert same size

%OJO Recordings depends on type of source at location r in sim2 during Kernel computation
%If I'm using expl src in sim 1 => DCexp formula is based on recording pressure (~DivU). Additionally, I use dvp only.
%If I'm using ElForce in x as src in sim 1 => DCexp formula is based on recording ux.


%Load Perturbed signal - i.e., Sim1_pert:  dv_s only or dv_p only
%dv_p only
%w_pert=load(['/cluster/scratch/javierd/DCexp_EL_dv100/Sim1_pert_r',num2str(rp,'%01.0f'),'/dvp_only/M',num2str(M,'%01.0f'),'/OUTPUT_FILES/S',num2str(rcv,'%04.0f'),'.AA.BHX.semc']);  

%dv_s only
%w_pert=load(['/cluster/scratch/javierd/DCexp_EL_dv100/Sim1_pert_r',num2str(rp,'%01.0f'),'/dvs_only/M',num2str(M,'%01.0f'),'/OUTPUT_FILES/S',num2str(rcv,'%04.0f'),'.AA.BHX.semc']);  

    
            
        else  %Vp and Vs perturbation together
            
        %Load Non-perturbed signal - i.e., Sim1
w=load(['/cluster/scratch/javierd/LargeMedium/Sim1/M',num2str(M,'%01.0f'),'/OUTPUT_FILES/S',num2str(rcv,'%04.0f'),'.AA.BHX.semc']);  
w=w(1:nt,:); %to make w and w_pert same size



%Load Perturbed signal - i.e., Sim1_pert
w_pert=load(['/cluster/scratch/javierd/LargeMedium/Sim1_pert_rp3/Both_dv/M',num2str(M,'%01.0f'),'/OUTPUT_FILES/S',num2str(rcv,'%04.0f'),'.AA.BHX.semc']);  



        end
        
        
for t=1+twm:length(w(:,1))-twm            
    
    
    Psum=0; Psum2=0; Psum3=0;
for tau=t-twm:t+twm
  
    %--numerator integral
    %Waves:
    Ptemp=w(tau,2).* w_pert(tau,2)+Psum;         
    % temp=abs ( hilbert(w(tau,2).* w_pert(tau,2)) ) ;   
   Psum=Ptemp;
   
     
   
   %--denominator integral
   %Waves:
   Ptemp_2=(w(tau,2).^2)+Psum2;
   %Ptemp_2=abs ( hilbert(w(tau,2)) ).^2 +Psum2;
   Psum2=Ptemp_2;
   
   Ptemp_3=(w_pert(tau,2).^2)+Psum3;
   %temp_3=abs ( hilbert(w_pert(tau,2)) ).^2 +sum2;
   Psum3=Ptemp_3;
   
   
end

%Norm_cc(t-twm)=(sum2.*dt).*(sum3.*dt);

CC(t-twm,M) = [Psum.*dt]./[sqrt( (Psum2.*dt).*(Psum3.*dt) )];


% figure 
% plot(1-CC,'k')




        end
   
%--ACOUSTIC----
    else
%Load Non-perturbed signal - i.e., Sim1
w=load(['/Users/alejandro/KernelsComparison_paper/Structural_Change/Decorrelation/Kernel_compt_Acoustic/Sim1/M',num2str(M,'%01.0f'),'/S',num2str(rcv,'%04.0f'),'.AA.PRE.semp']);

%Load Perturbed signal - i.e., Sim1_pert
w_pert=load(['/Users/alejandro/KernelsComparison_paper/Structural_Change/Decorrelation/Kernel_compt_Acoustic/Sim1_Perturbed_r',num2str(rp,'%01.0f'),'/M',num2str(M,'%01.0f'),'/S',num2str(rcv,'%04.0f'),'.AA.PRE.semp']);
    

for t=1+twm:length(w(:,1))-twm            
    
    sum=0; sum2=0; sum3=0;
for tau=t-twm:t+twm

    %numerator integral
    temp=w(tau,2).* w_pert(tau,2)+sum;         
    % temp=abs ( hilbert(w(tau,2).* w_pert(tau,2)) ) ;   
   sum=temp;

   
   %denominator integral
   temp_2=(w(tau,2).^2)+sum2;
   %temp_2=abs ( hilbert(w(tau,2)) ).^2 +sum2;
   sum2=temp_2;
   
   temp_3=(w_pert(tau,2).^2)+sum3;
   %temp_3=abs ( hilbert(w_pert(tau,2)) ).^2 +sum2;
   sum3=temp_3;
   


end

%Norm_cc(t-twm)=(sum2.*dt).*(sum3.*dt);

CC(t-twm,M) = [sum.*dt]./[sqrt( (sum2.*dt).*(sum3.*dt) )];

end

    end  %End If 


end

 if strcmp(medium,'Elastic')==1
           
DC=1-CC;
DC_P=DC;

%----%Just change name accordingly with what you place as perturbation: dv_P or dv_S or together----
% Average over models
%CC_mean=mean(CC,2);
DC_expP=mean(DC_P,2);    

%Standard Deviation of the Decorrelation
Sp =std(DC_P,0,2);    %std(A,n,2) if n=0 =>normalized over N-1 , if n=1 normalized over N

%Standard Error
SE_p=Sp./sqrt(nmodels);

%---------OJO!!-------------
%S-wave has same formula as P-case. It comes from having expl source at r  
%in sim2. Thus, I just try to separate the labels only for saving proposes
DC_S=DC_P;  DC_expS=DC_expP;   Ss=Sp;  SE_s=SE_p;

% In case of both pert together
DC_both=DC_P;  DC_expBoth=DC_expP;   S_both=Sp;  SE_both=SE_p; 

   
%--Save DC_mean ELASTIC-- Recording Displacement in x :  ux
%dvp_only
%folder_save=(['/cluster/scratch/javierd/Epsilon_saveESP/Epsilon_Layers/Weak_pert/Anne_Config/']);  
%save([folder_save,'DCexp_P_ux_Elastic_Layer_twm',num2str(twm,'%01.0f'),'.mat'],'DC_expP','Sp','SE_p','DC_P')


%dvs_only
%folder_save=(['/cluster/scratch/javierd/Epsilon_saveESP/Epsilon_Layers/Weak_pert/Anne_Config/']);  
%save([folder_save,'DCexp_S_ux_Elastic_Layer_twm',num2str(twm,'%01.0f'),'.mat'],'DC_expS','Ss','SE_s','DC_S')

%Both
folder_save=(['/cluster/scratch/javierd/LargeMedium/Sim1_pert_rp3/Both_dv/']);  
save([folder_save,'DCexp_BothVel_ux_Elastic_LargeMed_twm',num2str(twm,'%01.0f'),'.mat'],'DC_expBoth','S_both','SE_both','DC_both')  


     
 else %--Acoustic--

  
DC=1-CC;

% Average over models
%CC_mean=mean(CC,2);
DC_exp=mean(DC,2);

%Standard Deviation of the Decorrelation
S =std(DC,0,2);    %std(A,n,2) if n=0 =>normalized over N-1 , if n=1 normalized over N
 

%Standard Error
%SE=S./sqrt(nmodels);

%--Save DC_mean ACOUSTIC--
folder_save=(['/Users/alejandro/KernelsComparison_paper/Structural_Change/Decorrelation/Kernel_compt_',medium,'/']);
%save([folder_save,'DCmean60Models_rp',num2str(rp,'%01.0f'),'_twm',num2str(twm,'%01.0f'),'.mat'],'DC_exp','S','DC')

 
 end
 
%end %End for all rp

