% Code for fast DRPCA
close all;
clear all;
clc;

%% Add Path
running_folder = 'C:\Users\dpham\ownCloud\Working\Atempo\';
addpath(genpath(fullfile(running_folder,'ISBI_2020\fast_BDRPCA-GitHub')));
%% Parameters
test = 3;
mm=0;
loops=20;
%% Testing Case
if test ==1
    nomfichier='simu_conv' 
    nomfichierS='simu';
    mm=1;
    lambda=0.1;
    tolS = 1e-3;
    rang0=7;
elseif test ==2
    nomfichier='cerveau_sain'
    nomfichierS='cerveau'; 
    lambda=5e-5;%*1./sqrt(max(Nz*Nx,Nt));
    tolS = 1e-6;  
    rang0=50;
elseif test ==3
    nomfichier='peri' 
    nomfichierS='peri';
    lambda=8e-5;%1.3/sqrt(max(Nz*Nx,Nt));
    tolS = 1e-6;
    rang0=80;
else
    nomfichier='tumeur'
    nomfichierS='tumeur';
    lambda=1.2*1e-4;
    tolS = 1e-6;
    rang0=80;
end
%%
timeOverall = tic;    
result_folder = fullfile(running_folder,'ISBI_2020\fast_BDRPCA-GitHub','Results',sprintf('%s',nomfichier));
mkdir(result_folder)
%% Loading data
iHS=0; % 
load_data_US;

%%% For matrix-based algorithms
[M,m,n,p] = convert_video3d_to_2d(M1);
M = M/abs(max(M(:)));
%%
fprintf(sprintf('performing DRPCA...\n'));
%% FPCP-working

%%
fprintf(1,'Rang not specified. Trying to guess ...\n');
rang0 = guessRank(M) ;
fprintf(1,'Using Rank : %d\n',rang0);

%% New-
[L,S] =fastDRPCA(M, H, lambda, loops, rang0, tolS,[],[]);

%% AFFICHAGE DE L'IMAGE DEROULANTE SELON Nt APRES SEUILLAGE/FILTRAGE
Mfinale=reshape(S,Nz,Nx,Nt);
save(sprintf('%s/fDRPCA.mat', result_folder),'Mfinale')   
%%
%Doppler de puissance
FigFeatures.title=1;
FigFeatures.result_folder = result_folder;
FigFeatures.mm=0;
FigFeatures.bar=1;
FigFeatures.print=0;
FigFeatures.nomtest = 'fast_DRPCA';
Dopplerplot(Mfinale,espace_xx,espace_zz,test,FigFeatures); 
clear Mfinale 

if test ==1
    fig_simu
else
    fig_reel;
end
