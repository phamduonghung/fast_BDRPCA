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
%% Figure parameters
FigFeatures.title=1; % Figure title 0 ou 1
FigFeatures.result_folder = result_folder;
FigFeatures.mm=0; 
FigFeatures.bar=1; % Colorbar 0 or 1 
FigFeatures.print=0; % Pdf Figure Print: 0 or 1 through export_fig 
%% Loading data
load_data_US;
[M,m,n,p] = convert_video3d_to_2d(M1);
%% Performing fast BD-RPCA
tfBDRPCAStart = tic;           % pair 2: tic
% Initialization using SVD
fprintf(sprintf('performing SVD...\n'))
tSVDStart = tic;           % pair 2: tic
Mnew = M'*M                 ; %Matrice carr?e
[V,D2,Vt] = svd(Mnew)       ; %Application de la SVD
D = sqrt(D2)                ; %Matrice des valeurs singuli?res
U = M*V/D                   ; %Calcul de la matrice spatiale des vecteurs singuliers
fprintf('Number of singular values: %d\n', length(diag(D)))

f=ones(1,Nt)                    ; %creation d'un vecteur ones
f(seuil_tissu+1:Nt)=[0]            ; %Application du seuil tissu sur le vecteur 
If=diag(f)                      ; %Matrice diagonale identit? filtree par les seuils
T0=M*V*If*V'                    ; %Calcul de la matrice finale    

% Estimated initial PSF
fprintf('Running estimated initial PSF ....\n')
Mt = reshape(M-T0,Nz,Nx,Nt);
M11 = squeeze(mean(Mt,3));
[H,psf0] = Hestimate(M11,Nz,Nx,Nt);
fprintf('Initialized PSF size: %d-%d\n',size(psf0,1),size(psf0,2))
clear Mt M11 

%% Rank Guess rG
fprintf(1,'Rang not specified. Trying to guess ...\n');
rang0 = guessRank(M) ;
fprintf(1,'Using Rank : %d\n',rang0);

%% Stop condition
tol  = 1e-3;
xtmp = M;
normM = norm(M, 'fro');

loops=20;
lambda=0.05;
max_iter = 20;
err = zeros(1,max_iter);

for iter = 1:max_iter    
    fprintf('Running fBDRPCA for iteration %d....\n',iter)
    [T,x] =fastDRPCA(M, H, lambda, loops, rang0, tol,[],[]);   
        
    % Stop Condition
    Z1 = x-xtmp;    
    err(1,iter) = log(norm(Z1, 'fro')) / normM; 
    xtmp=x;       
    if (err(1,iter) >= tol)       
        Mt = reshape(M-T,Nz,Nx,Nt);
        M11 = squeeze(mean(Mt,3));
        fprintf('Running estimated PSF for iteration %d....\n',iter+1)
        [H,psf1] = Hestimate(M11,Nz,Nx,Nt); 
        fprintf('PSF size for iteration %d: %d-%d\n',iter+1,size(psf1,1),size(psf1,2))   
    end  
    if (err(1,iter) < tol) || (iter>=2 &&(err(1,iter)>err(1,iter-1)))
        break
    end
    clear Mt M11 psf1
end
tfBDRPCAEnd = toc(tfBDRPCAStart)      % pair 2: toc
%% AFFICHAGE DE L'IMAGE DEROULANTE SELON Nt APRES SEUILLAGE/FILTRAGE
Mfinale=reshape(x,Nz,Nx,Nt);
%save(sprintf('%s/fBDRPCA_%s.mat', result_folder,nomfichier),'Mfinale')

%% Doppler de puissance
FigFeatures.nomtest = sprintf('fBDRPCA_%s',nomfichier); % Name 
Dopplerplot(Mfinale,espace_xx,espace_zz,test,FigFeatures); 
clear Mfinale 
