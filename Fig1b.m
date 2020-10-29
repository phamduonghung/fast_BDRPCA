%%% code matlab of Fig1a: GoDec%%%%%                                              ;
clear  all;
close all 
%% Set Current Folder of MATLAB being BD-RPCA-GitHub and Add Path
addpath(genpath(fullfile(pwd)));

%% A modifier
test = 3;
metLS =2;

nomfichier='peri' 
seuil_tissu = 100;
seuil_bruit = 150;

result_folder = fullfile(pwd,'Results');
mkdir(result_folder)
%% Loading data
iHS=0; % Not run Oleg
load_data_US;
%%% For matrix-based algorithms
[M,m,n,p] = convert_video3d_to_2d(M1);
%M = M/abs(max(M(:)));

fprintf(sprintf('performing GoDec...\n'));

%% Rank Guess
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
FigFeatures.print=1;
FigFeatures.nomtest = 'GoDec';
Dopplerplot(Mfinale,espace_xx,espace_zz,test,FigFeatures); 
clear Mfinale 
