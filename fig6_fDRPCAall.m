%%%%% CODE COMPLET %%%%%                                               ;
% Change test = 1, 2 or 3 for RPCA, DRPCA or BD-BDPRCA - Checked OKAY for
% R1
%% INITIALISATION DE LA MATRICE 
%clear  all;
%clc
%close all 
cd C:\Users\dpham\Documents\GitHub\lrslibrary; 

running_folder = 'C:\Users\dpham\ownCloud\Working\Atempo';
%addpath(genpath(fullfile(running_folder, 'simulation')));
for test =0:4
    % 
    % running_folder = 'C:\Users\dpham\ownCloud\Working\Atempo\simulation';
    % addpath(genpath(running_folder));
    %result_folder = 'C:\Users\dpham\ownCloud\Working\Atempo\Latex_paper\figures\R1';
    %set(0,'DefaultAxesFontSize',16);

    %nomfichier='cerveau_sain' ;
    result_folder = fullfile(running_folder,'simulation\RAPID_test\Results_LSDe_BDeconv',sprintf('%s',nomfichier));

    [~,name,~] = fileparts(nomfichier);
    load(nomfichier)                                ; %chargement de la matrice
    M1=real(double(IQ));

    %%
    [Nz,Nx,Nt] = size(M1)                           ; %Attribution de la taille de la matrice RF

    %%
    if test ==0
        load(fullfile(running_folder,'simulation\Results',sprintf('%s',nomfichier),'SVD_simu.mat'));
        %nomfichier='SVD_cerveau'
    elseif test ==1
        load(fullfile(running_folder,'simulation\Results',sprintf('%s',nomfichier),'RPCA_simu.mat'));
        %nomfichier='RPCA_cerveau'
    elseif test==2
        load(fullfile(running_folder,'simulation\Results',sprintf('%s',nomfichier),'DRPCA_simu.mat'));
        %nomfichier='DRPCA_cerveau'        
    elseif test==3
        load(fullfile(running_folder,'simulation\Results',sprintf('%s',nomfichier),'BDRPCA_simu_2.mat'));
        %nomfichier='BDRPCA_cerveau'
    else
        load(fullfile(running_folder,'simulation\RAPID_test\Results_LSDe_BDeconv',sprintf('%s',nomfichier),'fDRPCA.mat'));
        %nomfichier='fRPCA_cerveau'
    end

    %% Doppler de puissance
    FigHandle0=figure();
    set(gcf,'Color',[1,1,1])                             ;
    %hold on
    colormap hot                                                    ;  
    Amp = 35                                                        ;
    PW = 1/Nt*sum(abs(Mfinale).^2,3)                                ;
    PWTD = PW/max(max(PW))                                          ;
    PWTD1 = max(PWTD,10^(-Amp/10));
    logPWTD = 10*log10(max(PWTD,10^(-Amp/10)))                      ;
    hh = imagesc(espace_x,espace_z,logPWTD(:,:),[-Amp,0])      ;
    axis ij equal tight                                             ; 
    set(gca, 'FontSize', 18, 'fontName','Arial','LineWidth',1.5)    ;
    xlabel('X [mm]','FontSize',18)                                  ; 
    ylabel('Z [mm]','FontSize',18)                                  ;      
    %h=colorbar                                                      ;
    %xlabel(h,'dB')                          ;
    drawnow              ;
    hold on
    difx = diff(espace_x);
    difz = diff(espace_z);
    xx = 73;
    zz = 95;
    dis = 30;
    SZ = 13;
    SX = 12;
    rectangle('Position',[espace_x(xx) espace_z(zz) SX*difx(1) SZ*difz(1)],'EdgeColor','w','LineWidth',3)
    rectangle('Position',[espace_x(1) espace_z(1) SX*difx(1) SZ*difz(1)],'EdgeColor','g','LineWidth',3)
    xt = [espace_x(xx)+SX*difx(1)*1.2 espace_x(10)];
    yt = [espace_z(zz)+SZ*difz(1)*1.5 4*espace_z(1)];
    if test==0
        %%
        % 
        % $$e^{\pi i} + 1 = 0$$
        % 
        xarrow = [0.25 0.33];    % adjust length and location of arrow 
        yarrow = [0.9 0.9];      % adjust hieght and width of arrow
        annotation('textarrow',xarrow,yarrow,'FontSize',13,'Linewidth',2,'Color','g')
        xarrow = [0.225 0.225];    % adjust length and location of arrow 
        yarrow = [0.88 0.80];      % adjust hieght and width of arrow
        annotation('textarrow',xarrow,yarrow,'FontSize',13,'Linewidth',2,'Color','g')
        str = {'R1','R2'};
        text(xt(1),yt(1),str(1),'FontSize',16,'Color','w')
        text(xt(2),yt(2),str(2),'FontSize',16,'Color','g')
    end
    
    N11 = 13*ones(1,20);
    N22 = 12*ones(1,16);
    PWTD2 = mat2cell(PWTD,N11,N22);
    PWTD2 = reshape(PWTD2,20*16,[]);

    % old
%     R11 = abs(PWTD(zz:zz+25,xx:xx+20));
%     R21 = abs(PWTD(zz:zz+25,xx+dis:xx+dis+20));
%     CR(test+1) = 10*log10(abs(mean(R11(:))-mean(R21(:))))
%     CNR(test+1) = round(abs(mean(R11(:))-mean(R21(:)))/sqrt(std2(R11)^2+std2(R21)^2),3)
    %
    R21 = PWTD(zz:zz+SZ-1,xx:xx+SX-1);  % fixed position
    for i=1:length(PWTD2)
        R11 = PWTD2{i};%PWTD1(zz:zz+SZ,xx+dis:xx+dis+SX);
        %CR(test+1) = abs(mean(R11(:))-mean(R21(:)))
        CR(i,test+1) = 20*log10(abs(mean(R11(:))/mean(R21(:))));
        %CNR(i,test+1) = abs(mean(R11(:))-mean(R21(:)))/sqrt(var(R11(:))+var(R21(:)));     
        %CNR(i,test+1) = abs(mean(R11(:))-mean(R21(:)))/sqrt(std2(R11)^2+std2(R21)^2);
        [CNR(i,test+1), CNR_max(i,test+1), mu_1(i,test+1), mu_2(i,test+1), var_1(i,test+1), var_2(i,test+1)] = calc_CNR(R11,R21);
    end
    
 
    if 1
    if test==0
    export_fig(FigHandle0, ... % figure handle
                sprintf('%s/%s_svd.pdf', result_folder,nomfichier),... % name of output file without extension
                '-painters', ...      % renderer
                '-transparent', ...   % renderer
                '-pdf', ...         % file format
                '-r5000');             % resolution in dpi 
    end

    if test==1
    export_fig(FigHandle0, ... % figure handle
                sprintf('%s/%s_rpca.pdf', result_folder,nomfichier),... % name of output file without extension
                '-painters', ...      % renderer
                '-transparent', ...   % renderer
                '-pdf', ...         % file format
                '-r5000');             % resolution in dpi 
    end
    %% Code for paper   
    if test==2
    export_fig(FigHandle0, ... % figure handle
                sprintf('%s/%s_drpca.pdf', result_folder,nomfichier),... % name of output file without extension
                '-painters', ...      % renderer
                '-transparent', ...   % renderer
                '-pdf', ...         % file format
                '-r5000');             % resolution in dpi      
    end
    if test==3
    export_fig(FigHandle0, ... % figure handle
                sprintf('%s/%s_bdrpca.pdf', result_folder,nomfichier),... % name of output file without extension
                '-painters', ...      % renderer
                '-transparent', ...   % renderer
                    '-pdf', ...         % file format
                    '-r5000');             % resolution in dpi      
    end
    if test==4
    export_fig(FigHandle0, ... % figure handle
                sprintf('%s/%s_fbdrpca.pdf', result_folder,nomfichier),... % name of output file without extension
                '-painters', ...      % renderer
                '-transparent', ...   % renderer
                    '-pdf', ...         % file format
                    '-r5000');             % resolution in dpi      
    end
    end
end
%[CR1,ICR1]= maxk(CR,100);
%close all
CR1 = CR;
FigHandle(1) = figure(); %set(FigHandle(1),'units','normalized','outerposition',[0 0 1 1]);
%set(FigHandle(1),'Position',[500 50 700 600])
set(FigHandle(1), 'DefaultTextFontSize', 16);
%boxplot([CR1(:,1),CR1(:,2),CR1(:,3),CR1(:,4),CR1(:,5)],'Notch','on','Labels',{'SVD','RPCA','DRPCA','BD-RPCA','fast BRPCA'})
boxplot([CR1(:,2),CR1(:,3),CR1(:,4),CR1(:,5)],'Notch','on','Labels',{'RPCA','DRPCA','BD-RPCA','fast BD-RPCA'})

% txt1 = text(1.35,4,'8.65');
% set(txt1,'Rotation',90);
% txt1 = text(2.35,9,'16.29');
% set(txt1,'Rotation',90);
% txt1 = text(3.35,38,'43.74');
% set(txt1,'Rotation',90);
% txt1 = text(4.35,38,'43.36');
% set(txt1,'Rotation',90);
ylabel('CR [dB]')
set(gca,'Fontsize',14);

[CNR1,ICNR1] = maxk(CNR,100);
%CNR1 = CNR;
FigHandle(2) = figure(); %set(FigHandle(1),'units','normalized','outerposition',[0 0 1 1]);
%set(FigHandle(2),'Position',[500 50 600 700])
set(FigHandle(2), 'DefaultTextFontSize', 20);
boxplot([CNR1(:,1),CNR1(:,2),CNR1(:,3),CNR1(:,4),CNR1(:,5)],'Notch','on','Labels',{'SVD','RPCA','DRPCA','BD-RPCA','fast BD-RPCA'})
ylabel('CNR [dB]')


if 1
 for i=1:2
     export_fig(FigHandle(i), ... % figure handle
                sprintf('%s/%s_CR_CNR_%d.pdf', result_folder,nomfichier,i),... % name of output file without extension
                '-painters', ...      % renderer
                '-transparent', ...   % renderer
                '-pdf', ...         % file format
                '-r5000');             % resolution in dpi    
 end
end