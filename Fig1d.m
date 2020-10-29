%%% code matlab of Fig 1d: fast BD-RPCA%%%%%                                              ;
close all 
clear all
%% Set Current Folder of MATLAB being BD-RPCA-GitHub and Add Path
addpath(genpath(fullfile(pwd)));

%% Some parameters
test = 1; % For figure 2a of the paper, keep test=1 
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

%% fast BDRPCA
%% SSGoDec 
tau = 0.01;
power = 1;
tGoDecStart = tic;   
[T0,~,~,~]=SSGoDec(M,rang0,tau,power);
tGoDecEnd = toc(tGoDecStart)      % pair 2: toc

if 0
    tfBDRPCAStart = tic;           % pair 2: tic
    % Initialization SVD
    tSVDStart = tic;        
    fprintf('Initialization SVD....\n')
    Mnew = M'*M                 ; %Matrice carr?e
    [V,D2,Vt] = svd(Mnew)       ; %Application de la SVD
    D = sqrt(D2)                ; %Matrice des valeurs singuli?res
    U = M*V/D                   ; %Calcul de la matrice spatiale des vecteurs singuliers
    fprintf('Number of singular values: %d\n', length(diag(D)))

    f=ones(1,Nt)                    ; %cr?ation d'un vecteur ones
    f(seuil_tissu+1:end)=[0]            ; %Application du seuil tissu sur le vecteur 
    If=diag(f)                      ; %Matrice diagonale identit? filtr?e par les seuils
    T0=M*V*If*V'                    ; %Calcul de la matrice finale    
    tSVDEnd = toc(tSVDStart)      % pair 2: toc
end
%% Lambda1 Parameters
%  if test==1
%     Lambda = 3./sqrt(max(Nz*Nx,Nt));
%     Lambda1 = 1./sqrt(max(Nz*Nx,Nt));
% elseif test==2
%     Lambda = 1.3*1./sqrt(max(Nz*Nx,Nt));
%     Lambda1 = 1*1./sqrt(max(Nz*Nx,Nt));                       
% elseif test==3
%     %Lambda = 0.9*1./sqrt(max(Nz*Nx,Nt));
%     Lambda = 0.9*1./sqrt(max(Nz*Nx,Nt));
%     Lambda1 = 1.15./sqrt(max(Nz*Nx,Nt));
%     lambda = 4*1e3;% 1.15./sqrt(max(Nz*Nx,Nt));
% else         
%     Lambda = 1.2*1./sqrt(max(Nz*Nx,Nt));
%     Lambda1 = 1.1*1./sqrt(max(Nz*Nx,Nt)); 
% end

% tRPCAStart = tic;           % pair 2: tic
% fprintf('Initialization RPCA....\n')
% [T0, B0] = RobustPCA_Doppler(M,Lambda); %
% tRPCAEnd = toc(tRPCAStart)      % pair 2: toc
%T0 = T0/abs(max(T0(:)));

% if 0
%     Mfinale=reshape(T0,Nz,Nx,Nt);
%     FigFeatures.title=1;
%     FigFeatures.result_folder = result_folder;
%     FigFeatures.mm=0;
%     FigFeatures.bar=1;
%     FigFeatures.print=0;
%     FigFeatures.nomtest = 'RPCA';
%     Dopplerplot(Mfinale,espace_xx,espace_zz,test,FigFeatures); 
%     clear Mfinale 
% end
%%
fprintf('Running estimated initial PSF ....\n')
max_iter = 5;
Mt = reshape(M-T0,Nz,Nx,Nt);
M11 = squeeze(mean(Mt,3));
[H,psf0] = Hestimate(M11,Nz,Nx,Nt);
fprintf('Initialized PSF size: %d-%d\n',size(psf0,1),size(psf0,2))
clear Mt M11 

%% Stop condition
tol  = 1e-3;
xtmp = M;
Ttmp = T0;
err = zeros(1,max_iter);
normM = norm(M, 'fro');

loops=20;
lambda=0.1;
%tolS = 1e-3;

for iter = 1:max_iter
    iter
    fprintf('Running estimated DRPCA for iteration %d....\n',iter)
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
    else 
        break;
    end    
    clear Mt M11 psf1
end
tBDRPCAEnd = toc(tfBDRPCAStart)      % pair 2: toc

%% AFFICHAGE DE L'IMAGE DEROULANTE SELON Nt APRES SEUILLAGE/FILTRAGE
Mfinale=reshape(x,Nz,Nx,Nt);
%save(sprintf('%s/fBDRPCA.mat', result_folder),'Mfinale')   
%% Doppler de puissance
FigFeatures.nomtest = 'fBDRPCA';
Dopplerplot(Mfinale,espace_xx,espace_zz,test,FigFeatures); 
clear Mfinale 

