%%% code matlab of Fig 1d: fast BD-RPCA%%%%%                                              ;
close all 
clear all
%% Set Current Folder of MATLAB being BD-RPCA-GitHub and Add Path
addpath(genpath(fullfile(pwd)));

%% Some parameters
test = 1; % For figure 1d of the paper, keep test=1 
nomfichier='simu_conv' 
seuil_tissu = 2;
seuil_bruit = 15;
result_folder = fullfile(pwd,'Results');
mkdir(result_folder)
%% Some figure parameters
FigFeatures.title=1;
FigFeatures.result_folder = result_folder;
FigFeatures.mm=0;
FigFeatures.bar=1;
FigFeatures.print=0;
%% Loading data
load_data_US;
[M,m,n,p] = convert_video3d_to_2d(M1);
%% Rank Guess
fprintf(1,'Rang not specified. Trying to guess ...\n');
rang0 = guessRank(M) ;
fprintf(1,'Using Rank : %d\n',rang0);

% %% SSGoDec 
% tau = 0.025;
% power = 1;
% tGoDecStart = tic;   
% [T1,~,~,~]=SSGoDec(M,rang0,tau,power);
% tGoDecEnd = toc(tGoDecStart)      % pair 2: toc

%% Initialization using PRCA or load directly T0 from Data folder
%tRPCAStart = tic;           % pair 2: tic
%fprintf('Initialization RPCA....\n')
%[T0, ~] = RobustPCA_Doppler(M,Lambda); %
%tRPCAEnd = toc(tRPCAStart)      % pair 2: toc
%load(fullfile(pwd,'Data','T0.mat')) ; 
%save(sprintf('%s/T0.mat', result_folder),'T0')   

%% SSGoDec 
tau = 0.025;
power = 1;
tGoDecStart = tic;   
[T0,X0,~,~]=SSGoDec(M,rang0,tau,power);
tGoDecEnd = toc(tGoDecStart)      % pair 2: toc

%%
tfBDRPCAStart = tic;  
fprintf('Running estimated initial PSF ....\n')
max_iter = 5;
Mt = reshape(M-T0,Nz,Nx,Nt);
M11 = squeeze(mean(Mt,3));
[H,psf0] = Hestimate(M11,Nz,Nx,Nt);
fprintf('Initialized PSF size: %d-%d\n',size(psf0,1),size(psf0,2))
clear Mt M11 
% Print H figures
FigHandlem=figure();
imagesc(abs(psf0)); colorbar; 
%% Stop condition
tol  = 1e-3;
xtmp = M;
Ttmp = T0;
err = zeros(1,max_iter);
normM = norm(M, 'fro');

loops=20;
lambda=0.05;
iter = 1
for iter = 1:max_iter    
    fprintf('Running BDRPCA for iteration %d....\n',iter)
    [T,x] =fastDRPCA(M, H, lambda, loops, rang0, tol,[],[]);
    %% AFFICHAGE DE L'IMAGE DEROULANTE SELON Nt APRES SEUILLAGE/FILTRAGE
    Mfinale=reshape(x,Nz,Nx,Nt);
    FigFeatures.nomtest = sprintf('B_image-Iter_%d',iter);
    Dopplerplot(Mfinale,espace_xx,espace_zz,test,FigFeatures); 
    
    % Stop Condition
    Z1 = x-xtmp;    
    err(1,iter) = norm(Z1, 'fro') / normM  
    xtmp=x;       
    if (err(1,iter) > tol)    
        Mt = reshape(M-T,Nz,Nx,Nt);
        M11 = squeeze(mean(Mt,3));
        fprintf('Running estimated PSF for iteration %d....\n',iter+1)
        [H,psf1] = Hestimate(M11,Nz,Nx,Nt); 
        fprintf('PSF size for iteration %d: %d-%d\n',iter+1,size(psf1,1),size(psf1,2))   
         % Print H figures
        FigHandlem=figure();
        imagesc(abs(psf1)); colorbar; 
    else 
        break;
    end    
    clear Mt M11 psf1
    pause
    close all
end
tBDRPCAEnd = toc(tfBDRPCAStart)      % pair 2: toc
%% AFFICHAGE DE L'IMAGE DEROULANTE SELON Nt APRES SEUILLAGE/FILTRAGE
Mfinale=reshape(x,Nz,Nx,Nt);
%save(sprintf('%s/fBDRPCA_%s.mat', result_folder,nomfichier),'Mfinale')

%% Doppler de puissance
FigFeatures.nomtest = sprintf('BDRPCA_%s',nomfichier); % Name 
Dopplerplot(Mfinale,espace_xx,espace_zz,test,FigFeatures); 
clear Mfinale 
