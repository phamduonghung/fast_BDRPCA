% Code for fast DRPCA
clear  all;
close all 
%% Set Current Folder of MATLAB being BD-RPCA-GitHub and Add Path
addpath(genpath(fullfile(pwd)));

%% Some parameters
test = 1; % For figure 2a of the paper, keep test=1 
nomfichier='simu_conv' 
result_folder = fullfile(pwd,'Results');
mkdir(result_folder)
%% Loading data
load_data_US;
%% For matrix-based algorithms
[M,m,n,p] = convert_video3d_to_2d(M1);
%M = M/abs(max(M(:)));

%% Rank Guess
fprintf(1,'Rang not specified. Trying to guess ...\n');
rang0 = guessRank(M) ;
fprintf(1,'Using Rank : %d\n',rang0);

%% Fast DRPCA
lambda=0.1;
loops=20;
%% BDRPCA running
tfDRPCAStart = tic;           % pair 2: tic
fprintf('Running fast DRPCA....\n')
[L,S] =fastDRPCA(M, H, lambda, loops, rang0, [],[],[]);
tfDRPCAEnd = toc(tfDRPCAStart)     

%% AFFICHAGE DE L'IMAGE DEROULANTE SELON Nt APRES SEUILLAGE/FILTRAGE
Mfinale=reshape(S,Nz,Nx,Nt);
%save(sprintf('%s/fDRPCA.mat', result_folder),'Mfinale')   

%% Figures Parameters 
FigFeatures.title=1; % Figure title 0 ou 1
FigFeatures.result_folder = result_folder;
FigFeatures.mm=0; 
FigFeatures.bar=1; % Colorbar 0 or 1 
FigFeatures.print=0; % Pdf Figure Print 0 or 1 through export_fig 
FigFeatures.nomtest = sprintf('fast_DRPCA_%s',nomfichier); % Name 
Dopplerplot(Mfinale,espace_xx,espace_zz,test,FigFeatures); 
clear Mfinale