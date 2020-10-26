%% LRSLibrary: A library of low-rank and sparse tools for background/foreground separation in videos
close all
clear all
% restoredefaultpath;

%% First run the setup script
%lrs_setup; % or run('C:/GitHub/lrslibrary/lrs_setup')
%% LRS GUI (graphical user interface)
%lrs_gui;
%% Load configuration (for demos)
%lrs_load_conf;
%% Auxiliary functions
%%% Load video
% input_avi = fullfile(lrs_conf.lrs_dir,'dataset','demo.avi');
% video = load_video_file(input_avi);
% show_video(video);

%%% Convert video to MAT (MATLAB format)
% output_mat = fullfile(lrs_conf.lrs_dir,'dataset','demo.mat');
% video2mat(input_avi, output_mat);

%%% Convert video to matrix representation (2D array)
% M = im2double(convert_video_to_2d(video));
% show_2dvideo(M,video.height,video.width);

%%% Convert video to 3D array
% V = im2double(convert_video_to_3d(video));
% show_3dvideo(V);

%%% Convert video to 3D tensor
% T = convert_video_to_3dtensor(video);
% tensorlab.slice3(double(T)), colormap('gray'); % Visualize a third-order tensor with slices.

%%% Convert video to 4D array/tensor
% video4d = convert_video_to_4d(video);
% video4d = crop_4dvideo(video4d,1,10);
% video4d = resize_4dvideo(video4d,2);
% show_4dvideo(video4d);

%% DEMO 01 (process matrix/tensor data)
running_folder = 'C:\Users\dpham\ownCloud\Working\Atempo';
addpath(genpath(fullfile(running_folder,'simulation')));
%addpath(genpath(fullfile(running_folder,'simulation','RAPID_test')));
%addpath(genpath('C:\Users\dpham\Documents\RAPID'));
%% A modifier
test = 3;
metLS =2;
mm=0;

%%
timeOverall = tic;
if test ==1
    nomfichier='simu_conv' 
    nomfichierS='simu';
    %lambdaL = 1/sqrt(max(size(M)));
    %lambdaS=1e-1;
    mm=1;
elseif test ==2
    nomfichier='cerveau_sain'
    nomfichierS='cerveau';
    %lambdaL = 1e-5; %1/sqrt(max(size(M))); 
    %lambdaS = 6e-5; %1/sqrt(max(size(M))); 
elseif test ==3
    nomfichier='peri' 
    nomfichierS='peri';
    %lambdaL = 0;
    %lambdaS=6e3;
else 
    nomfichier='tumeur'
    nomfichierS='tumeur';
    %lambdaL = 0;
    %lambdaS=5e-5;
end
result_folder = fullfile(running_folder,'simulation\RAPID_test\Results_LSDe_BDeconv',sprintf('%s',nomfichier));
mkdir(result_folder)
%% Loading data
iHS=0; % Not run Oleg
load_data_US;
%test_ClutterPP_pcpDeconv('simu_conv',0.5, [], 0, 1,[],[],[1 10]);
% H=psf;
%test_ClutterPP_pcpDeconv('cerveau_sain',1e4, psf, 0, 1,[],[],[50, 300]);
%load(fullfile(lrs_conf.lrs_dir,'dataset','trafficdb','traffic_patches.mat'));
%V = im2double(imgdb{100});
%show_3dvideo(M1);

%%% For matrix-based algorithms
[M,m,n,p] = convert_video3d_to_2d(M1);
M = M/abs(max(M(:)));
%[DR,M]=dynarange(M); 
% show_2dvideo(M,m,n);
%% LS DECOMPOSITION
disp('LS Decomposition starts...');
timeLS = tic;

fprintf(sprintf('performing Fastpcp...\n'));
%%%%%% FPCP-working
loops=20;
if test ==1
    lambda=0.1;%25./sqrt(max(Nz*Nx,Nt));
    tolS = 1e-3;
    rang0=7;
elseif test ==2
    lambda=5e-5;%*1./sqrt(max(Nz*Nx,Nt));
    tolS = 1e-6;  
    rang0=50;
elseif test ==3
    lambda=8e-5;%1.3/sqrt(max(Nz*Nx,Nt));
    tolS = 1e-6;
    rang0=80;
else
    lambda=1.2*1e-4
    %25./sqrt(max(Nz*Nx,Nt));
    tolS = 1e-6;
    rang0=80;
end

fprintf(1,'Rang not specified. Trying to guess ...\n');
rang0 = guessRank(M) ;
fprintf(1,'Using Rank : %d\n',rang0);

tfDRPCAStart = tic;           % pair 2: tic
[L,S] =fastpcp_new2(M, H, lambda, loops, rang0, tolS,[],[],Nz,Nx,Nt,espace_xx,espace_zz);
tBDRPCAEnd = toc(tfDRPCAStart) 
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
FigFeatures.nomtest = 'B_image';
Dopplerplot(Mfinale,espace_xx,espace_zz,test,FigFeatures); 
clear Mfinale 

if test ==1
    fig2_3_fDRPCA
else
    fig6_fDRPCAall;
end
