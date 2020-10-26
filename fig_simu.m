%%%%% Figures for simulation data%%%%%                                            ;

set(0,'DefaultAxesFontSize',16);

FigFeatures.result_folder = result_folder;
%FigFeatures.print = 1;
% FigFeatures.title=0;
FigFeatures.mm=1;
FigFeatures.bar=1;
% test = 1;
%% Testing method
FigFeatures.nomtest = 'SVD_simu';
T = load(sprintf('%s/simulation/Results/simu_conv/Old/%s.mat',running_folder,FigFeatures.nomtest));
T = T.Mfinale;
[Nz,Nx,Nt] =size(T);
espace_xx=1:Nx;
espace_zz=1:Nz;
FigFeatures.color1='c';
FigFeatures.color2= '--c';
[~,TlogPWTD]=Dopplerplot(T,espace_xx,espace_zz,test,FigFeatures);    

%% RPCA
FigFeatures.nomtest = 'RPCA_simu';
R = load(sprintf('%s/simulation/Results/simu_conv/Old/%s.mat',running_folder,FigFeatures.nomtest));
R = R.Mfinale;
[Nz,Nx,Nt] =size(R);
espace_xx=1:Nx;
espace_zz=1:Nz;
FigFeatures.color1='m';
FigFeatures.color2= '--m';
[~,RlogPWTD]=Dopplerplot(R,espace_xx,espace_zz,test,FigFeatures);    

%% DRPCA
FigFeatures.nomtest = 'DRPCA_simu';
DR = load(sprintf('%s/simulation/Results/simu_conv/Old/%s.mat',running_folder,FigFeatures.nomtest));
DR = DR.Mfinale;
[Nz,Nx,Nt] =size(DR);
espace_xx=1:Nx;
espace_zz=1:Nz;
FigFeatures.color1='r';
FigFeatures.color2= '--r';
[~,DRlogPWTD]=Dopplerplot(DR,espace_xx,espace_zz,test,FigFeatures);     
%% BDPRCA
FigFeatures.nomtest = 'BDRPCA_simu_3';
BDR = load(sprintf('%s/simulation/Results/simu_conv/Old/%s.mat',running_folder,FigFeatures.nomtest));
BDR = BDR.Mfinale;
[Nz,Nx,Nt] =size(BDR);
espace_xx=1:Nx;
espace_zz=1:Nz;
FigFeatures.color1='g';
FigFeatures.color2= '--g';
[~,BDRlogPWTD]=Dopplerplot(BDR,espace_xx,espace_zz,test,FigFeatures);   

%% fDPRCA
FigFeatures.nomtest = 'fDRPCA';
%FigFeatures.print = 0;
fBDR = load(sprintf('%s/fast_BDRPCA-GitHub/Results/%s/%s.mat',running_folder,nomfichier,FigFeatures.nomtest));
fBDR = fBDR.Mfinale;
[Nz,Nx,Nt] =size(fBDR);
espace_xx=1:Nx;
espace_zz=1:Nz;
FigFeatures.color1='y';
FigFeatures.color2= '--y';
[~,fBDRlogPWTD]=Dopplerplot(fBDR,espace_xx,espace_zz,test,FigFeatures);   

%% ideal simulation figure
FigFeatures.nomtest = 'ground_truth';
IlogPWTD = -35*zeros(Nz,Nx,Nt);
IlogPWTD(200:270,45:58,:)=1;
IlogPWTD(145:180,90:100,:)=1;
FigFeatures.color1='b';
FigFeatures.color2= '--b';
%FigFeatures.bar=1;
[~,MIlogPWTD]=Dopplerplot(IlogPWTD,espace_xx,espace_zz,test,FigFeatures);                 

%%
%close all
XNx = 1:Nx;        
FigH(1)=figure();
plot(XNx,MIlogPWTD(228,:),'-b','LineWidth',1.0)
hold on
plot(XNx,TlogPWTD(228,:),'-*c')
plot(XNx,RlogPWTD(228,:),'-om','LineWidth',1.0)
plot(XNx,DRlogPWTD(228,:),'-xr','LineWidth',1.0)
plot(XNx,BDRlogPWTD(228,:),'-dg')
plot(XNx,fBDRlogPWTD(228,:),'-+k')
set(gca,'xlim',[30 90]);
set(gca,'XTick', 30:20:90);
set(gca,'ylim',[-35 0]);
xlabel('N_X [nb]')                                  ; 
ylabel('Amplitude [dB]')                                  ; 
legend('ground truth','SVD','RPCA','DRPCA','BD-RPCA','fDRPCA')

XNz = 1:Nz;        
FigH(2)=figure();
plot(XNz,MIlogPWTD(:,95),'--b','LineWidth',1.0)
hold on
plot(XNz,TlogPWTD(:,95),'--*c')
plot(XNz,RlogPWTD(:,95),'--om','LineWidth',1.0)
plot(XNz,DRlogPWTD(:,95),'--xr','LineWidth',1.0)
plot(XNz,BDRlogPWTD(:,95),'--dg','LineWidth',1.0)
plot(XNz,fBDRlogPWTD(:,95),'--+k','LineWidth',1.0)
set(gca,'xlim',[110 260]);
set(gca,'XTick', 110:50:260);
set(gca,'ylim',[-35 0]);
xlabel('N_Z [nb]')                              ; 
ylabel('Amplitude [dB]')                                  ; 
legend('ground truth','SVD','RPCA','DRPCA','BD-RPCA','fDRPCA','Location','best')

if FigFeatures.print
    for i=1:2
        export_fig(FigH(i), ... % figure handle
                    sprintf('%s/cutoff_%d.pdf', result_folder,i),... % name of output file without extension
                    '-painters', ...      % renderer
                    '-transparent', ...   % renderer
                    '-pdf', ...         % file format
                    '-r5000');             % resolution in dpi    
    end
end

%%

FigH(1)=figure();
set(FigH(1),'units','normalized','outerposition',[0 0 1 1]);
subplot(2,3,1)
set(gcf,'Color',[1,1,1]);
[ssimval,ssimmap] = ssim(MIlogPWTD,TlogPWTD);
h = imagesc(ssimmap);
set(gca, 'FontSize', 14,'LineWidth',1.5)    ;
xlabel('N_X [nb]') ;
ylabel('N_Z [nb]')                                  ;
hh=colorbar                                                      ;
xlabel(hh,'dB')                          ;
title(['SSIM of SVD:',num2str(ssimval)])
set(gca,'XTick', 1:40:espace_xx(end));
set(gca,'YTick', 1:100:espace_zz(end));       

subplot(2,3,2)
set(gcf,'Color',[1,1,1]);
[ssimval,ssimmap] = ssim(MIlogPWTD,RlogPWTD);
h = imagesc(ssimmap);
set(gca, 'FontSize', 14,'LineWidth',1.5)    ;
xlabel('N_X [nb]') ;
ylabel('N_Z [nb]')                                  ;
hh=colorbar                                                      ;
xlabel(hh,'dB')                          ;
title(['SSIM of RPCA:',num2str(ssimval)])
set(gca,'XTick', 1:40:espace_xx(end));
set(gca,'YTick', 1:100:espace_zz(end));       

subplot(2,3,3)
set(gcf,'Color',[1,1,1]);
[ssimval,ssimmap] = ssim(MIlogPWTD,DRlogPWTD);
h = imagesc(ssimmap);
set(gca, 'FontSize', 14,'LineWidth',1.5)    ;
xlabel('N_X [nb]') ;
ylabel('N_Z [nb]')                                  ; 
hh=colorbar                                                      ;
xlabel(hh,'dB')                          ;
title(['SSIM of DRPCA:',num2str(ssimval)])
set(gca,'XTick', 1:40:espace_xx(end));
set(gca,'YTick', 1:100:espace_zz(end));       

subplot(2,3,4)
set(gcf,'Color',[1,1,1]);
[ssimval,ssimmap] = ssim(MIlogPWTD,BDRlogPWTD);
h = imagesc(ssimmap);
set(gca, 'FontSize', 14,'LineWidth',1.5)    ;
xlabel('N_X [nb]') ;
ylabel('N_Z [nb]')                                  ; 
hh=colorbar                                                      ;
xlabel(hh,'dB')                          ;
title(['SSIM of BDRPCA:',num2str(ssimval)])
Ssimval = zeros(Nt,4);                                
set(gca,'XTick', 1:40:espace_xx(end));
set(gca,'YTick', 1:100:espace_zz(end)); 

subplot(2,3,5)
set(gcf,'Color',[1,1,1]);
[ssimval,ssimmap] = ssim(MIlogPWTD,fBDRlogPWTD);
h = imagesc(ssimmap);
set(gca, 'FontSize', 14,'LineWidth',1.5)    ;
xlabel('N_X [nb]') ;
ylabel('N_Z [nb]')                                  ; 
hh=colorbar                                                      ;
xlabel(hh,'dB')                          ;
title(['SSIM of fPCA:',num2str(ssimval)])
Ssimval = zeros(Nt,4);                               
set(gca,'XTick', 1:40:espace_xx(end));
set(gca,'YTick', 1:100:espace_zz(end)); 

%% PSNR

errT = immse(TlogPWTD,MIlogPWTD);
errR = immse(RlogPWTD,MIlogPWTD);
errDR = immse(DRlogPWTD,MIlogPWTD);
errBDR = immse(BDRlogPWTD,MIlogPWTD);
errfBDR = immse(fBDRlogPWTD,MIlogPWTD);


[peaksnrT, snrT] = psnr(abs(TlogPWTD),abs(MIlogPWTD));
[peaksnrR, snrR] = psnr(abs(RlogPWTD),abs(MIlogPWTD));
[peaksnrDR, snrDR] = psnr(abs(DRlogPWTD),abs(MIlogPWTD));
[peaksnrBDR, snrBDR] = psnr(abs(BDRlogPWTD),abs(MIlogPWTD));
[peaksnrfBDR, snrfBDR] = psnr(abs(fBDRlogPWTD),abs(MIlogPWTD));


ResultT=CalcPerf(abs(MIlogPWTD),abs(TlogPWTD));
ResultR=CalcPerf(abs(MIlogPWTD),abs(RlogPWTD));
ResultDR=CalcPerf(abs(MIlogPWTD),abs(DRlogPWTD));
ResultBDR=CalcPerf(abs(MIlogPWTD),abs(BDRlogPWTD));
ResultfBDR=CalcPerf(abs(MIlogPWTD),abs(fBDRlogPWTD));

%%
[EQMT,PSNRT , RMSET, NRMSET] = US_ADM_calc_PSNR(abs(MIlogPWTD),abs(TlogPWTD));
[EQMR,PSNRR , RMSER, NRMSER] = US_ADM_calc_PSNR(abs(MIlogPWTD),abs(RlogPWTD));
[EQMDR,PSNRDR , RMSEDR, NRMSEDR] = US_ADM_calc_PSNR(abs(MIlogPWTD),abs(DRlogPWTD));
[EQMBDR,PSNRBDR , RMSEBDR, NRMSEBDR] = US_ADM_calc_PSNR(abs(MIlogPWTD),abs(BDRlogPWTD));
[EQMfBDR,PSNRfBDR , RMSEfBDR, NRMSEfBDR] = US_ADM_calc_PSNR(abs(MIlogPWTD),abs(fBDRlogPWTD));

%img1=abs(MIlogPWTD);
%sqrt(ResultBDR.MSE/(sum((img1(:)).^2))*numel(img1))
%format long

ResultUS = [[EQMT,PSNRT , RMSET, NRMSET];[EQMR,PSNRR , RMSER, NRMSER];[EQMDR,PSNRDR , RMSEDR, NRMSEDR] ;[EQMBDR,PSNRBDR , RMSEBDR, NRMSEBDR];[EQMfBDR,PSNRfBDR , RMSEfBDR, NRMSEfBDR]]

%format short

%% for 3D

for i=1:Nt
    [Ssimval(i,1),~] = ssim(abs(IlogPWTD(:,:,i)),abs(T(:,:,i)));
    [Ssimval(i,2),~] = ssim(abs(IlogPWTD(:,:,i)),abs(R(:,:,i)));
    [Ssimval(i,3),~] = ssim(abs(IlogPWTD(:,:,i)),abs(DR(:,:,i)));
    [Ssimval(i,4),~] = ssim(abs(IlogPWTD(:,:,i)),abs(BDR(:,:,i)));
    [Ssimval(i,5),~] = ssim(abs(IlogPWTD(:,:,i)),abs(fBDR(:,:,i)));
end
FigH(2)=figure;
boxplot(Ssimval,'Labels',{'SVD','PRCA','DRPCA','BD-PRCA','fBDR'})
ylabel('SSIM index')
set(gca,'ylim',[0.935 0.965]);

[Ssimval3D_SVD,~] = ssim(abs(IlogPWTD),abs(T))
[Ssimval3D_PRCA,~] = ssim(abs(IlogPWTD),abs(R))
[Ssimval3D_DPRCA,~] = ssim(abs(IlogPWTD),abs(DR))
[Ssimval3D_BDPRCA,~] = ssim(abs(IlogPWTD),abs(BDR))
[Ssimval3D_fBDPRCA,~] = ssim(abs(IlogPWTD),abs(fBDR))

if FigFeatures.print
    for i=1:2
        export_fig(FigH(i), ... % figure handle
                    sprintf('%s/boxplot_%d.pdf', result_folder,i),... % name of output file without extension
                    '-painters', ...      % renderer
                    '-transparent', ...   % renderer
                    '-pdf', ...         % file format
                    '-r5000');             % resolution in dpi    
    end
end
