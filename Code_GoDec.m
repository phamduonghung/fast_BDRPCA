% Code for GoDec
close all;
clear all;
%clc;

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
    %rang0=7;
elseif test ==2
    nomfichier='cerveau_sain'
    nomfichierS='cerveau'; 
    %rang0=50;
elseif test ==3
    nomfichier='peri' 
    nomfichierS='peri';
    %rang0=80;
else
    nomfichier='tumeur'
    nomfichierS='tumeur'
 
    %rang0=80;
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
%M = M/abs(max(M(:)));
%%
fprintf(sprintf('performing GoDec...\n'));

fprintf(1,'Rang not specified. Trying to guess ...\n');
rang0 = guessRank(M) ;
fprintf(1,'Using Rank : %d\n',rang0);

%% SSGoDec 
tau = 5e3;
power = 1;
tGoDecStart = tic;   
[L,S,RMSE,error]=SSGoDec(M,rang0,tau,power);
tGoDecEnd = toc(tGoDecStart)      % pair 2: toc
%% AFFICHAGE DE L'IMAGE DEROULANTE SELON Nt APRES SEUILLAGE/FILTRAGE
Mfinale=reshape(S,Nz,Nx,Nt);
%save(sprintf('%s/GoDec.mat', result_folder),'Mfinale')   

%% Doppler de puissance
% Figures Parameters 
FigFeatures.title=1;
FigFeatures.result_folder = result_folder;
FigFeatures.mm=0;
FigFeatures.bar=1;
FigFeatures.print=0;
FigFeatures.nomtest = 'GoDec';
Dopplerplot(Mfinale,espace_xx,espace_zz,test,FigFeatures); 
clear Mfinale 
