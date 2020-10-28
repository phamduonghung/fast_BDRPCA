%%% code matlab of BD-RPCA%%%%%                                              ;
%% INITIALISATION DE LA MATRICE 
clear  all;
%clc
close all 
%% Add Path
running_folder = 'C:\Users\dpham\ownCloud\Working\Atempo\';
addpath(genpath(fullfile(running_folder,'ISBI_2020\fast_BDRPCA-GitHub')));
%% A modifier
test = 3;
metLS =2;
mm=0;
%% A modifier
if test ==1
    nomfichier='simu_conv' 
    seuil_tissu = 2;
    seuil_bruit = 15;
elseif test ==2
    nomfichier='cerveau_sain'    
    seuil_tissu = 100;
    seuil_bruit = 150;
elseif test ==3
    nomfichier='peri' 
    seuil_tissu = 100;
    seuil_bruit = 150;
else 
    nomfichier='tumeur'
    seuil_tissu = 100;
    seuil_bruit = 200;
end
result_folder = fullfile(running_folder, sprintf('simulation/Results/%s',nomfichier));

%% Loading data
iHS=0; % Not run Oleg
load_data_US;

%% AFFICHAGE DE L'IMAGE INITIALE
% Some figure parameters
FigFeatures.title=0;
FigFeatures.result_folder = result_folder;
FigFeatures.mm=0;
FigFeatures.bar=1;
FigFeatures.print=0;

%% APPLICATION DE LA BDRPCA
M = reshape(M1(:),Nz*Nx,Nt) ; %Construction de la matrice de Casorati
%M = real(M);

tBDRPCAStart = tic;           % pair 2: tic
if 0
    % Initialization SVD
    tSVDStart = tic;           % pair 2: tic
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
 if test==1
    Lambda = 3./sqrt(max(Nz*Nx,Nt));
    Lambda1 = 1./sqrt(max(Nz*Nx,Nt));
elseif test==2
    Lambda = 1.3*1./sqrt(max(Nz*Nx,Nt));
    Lambda1 = 1*1./sqrt(max(Nz*Nx,Nt));                       
elseif test==3
    Lambda = 0.9*1./sqrt(max(Nz*Nx,Nt));
    Lambda1 = 1.15./sqrt(max(Nz*Nx,Nt));
else         
    Lambda = 1.2*1./sqrt(max(Nz*Nx,Nt));
    Lambda1 = 1.1*1./sqrt(max(Nz*Nx,Nt)); 
end

tRPCAStart = tic;           % pair 2: tic
fprintf('Initialization RPCA....\n')
[T0, ~] = RobustPCA_Doppler(M,Lambda); %
tRPCAEnd = toc(tRPCAStart)      % pair 2: toc
%%
fprintf('Running estimated PSF initiale....\n')
max_iter = 5;
Mt = reshape(M-T0,Nz,Nx,Nt);
M11 = squeeze(mean(Mt,3));
%M11=real(M11);
[H,psf0] = Hestimate(M11,Nz,Nx,Nt);
fprintf('Initialized PSF size: %d-%d\n',size(psf0,1),size(psf0,2))
clear Mt M11 

%% Stop condition
tol  = 1e-3;
xtmp = M;
Ttmp = T0;
err = zeros(1,max_iter);
normM = norm(M, 'fro');

for iter = 1:max_iter
    iter
    fprintf('Running estimated DRPCA for iteration %d....\n',iter)
    [T, x] = DRPCA(M,H,Lambda1); % S <-> B (blood) and  L <->T (tissue) and M <-> S  and H<-> D in paper        
                   
    % Stop Condition
    Z1 = x-xtmp;    
    err(1,iter) = norm(Z1, 'fro') / normM  
    xtmp=x;   
    
    if (err(1,iter) > tol)    
        Mt = reshape(M-T,Nz,Nx,Nt);
        M11 = squeeze(mean(Mt,3));
        %M11=real(M11);       
        fprintf('Running estimated PSF for iteration %d....\n',iter+1)
        %psftmp = zeros(size(psf));
        [H,psf1] = Hestimate(M11,Nz,Nx,Nt); 
        %save(sprintf('%s/HEst_%d.mat', result_folder,iter+1),'psf1')   
        fprintf('PSF size for iteration %d: %d-%d\n',iter+1,size(psf1,1),size(psf1,2))         
    else 
        break;
    end
    
    clear Mt M11 psf1

end
tBDRPCAEnd = toc(tBDRPCAStart)      % pair 2: toc
%% AFFICHAGE DE L'IMAGE DEROULANTE SELON Nt APRES SEUILLAGE/FILTRAGE
Mfinale=reshape(x,Nz,Nx,Nt);
%save(sprintf('%s/BDRPCA.mat', result_folder),'Mfinale')   
%%
%Doppler de puissance
FigFeatures.nomtest = 'BDRPCA';
Dopplerplot(Mfinale,espace_xx,espace_zz,test,FigFeatures); 
clear Mfinale 

