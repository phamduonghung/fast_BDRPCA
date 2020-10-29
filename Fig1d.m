%%% code matlab of fast BD-RPCA%%%%%                                              ;
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

%%
fprintf(sprintf('performing Fast BDRPCA...\n'));
%%%%%% Fast BDRPCA
loops=20;
if test ==1
    nomfichier='simu_conv' 
    lambda=0.1;%25./sqrt(max(Nz*Nx,Nt));
    tolS = 1e-3;
    rang0=7;
    seuil_tissu = 2;
    seuil_bruit = 15;
    mm=1;
elseif test ==2
    nomfichier='cerveau_sain'
    seuil_tissu = 100;
    seuil_bruit = 150;
    rank0 = 20;    
elseif test ==3
    nomfichier='peri' 
    %lambda=3e-5
    tolS = 1e-6;
    rang0=80;
    seuil_tissu = 100;
    seuil_bruit = 150;
else 
    nomfichier='tumeur'
    seuil_tissu = 100;
    seuil_bruit = 200;
    rang0=80;
end
result_folder = fullfile(running_folder,'ISBI_2020\fast_BDRPCA-GitHub','Results',sprintf('%s',nomfichier));
mkdir(result_folder)

%% Loading data
iHS=0; % Not run Oleg
load_data_US;
[M,m,n,p] = convert_video3d_to_2d(M1);
%M2=M;
%M = M/abs(max(M(:)));
%[DR,M]=dynarange(M); 
% show_2dvideo(M,m,n);

%%
fprintf(1,'Rang not specified. Trying to guess ...\n');
rang0 = guessRank(M) ;
fprintf(1,'Using Rank : %d\n',rang0);

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
    %Lambda = 0.9*1./sqrt(max(Nz*Nx,Nt));
    Lambda = 0.9*1./sqrt(max(Nz*Nx,Nt));
    Lambda1 = 1.15./sqrt(max(Nz*Nx,Nt));
    lambda = 4*1e3;% 1.15./sqrt(max(Nz*Nx,Nt));
else         
    Lambda = 1.2*1./sqrt(max(Nz*Nx,Nt));
    Lambda1 = 1.1*1./sqrt(max(Nz*Nx,Nt)); 
end

tRPCAStart = tic;           % pair 2: tic
fprintf('Initialization RPCA....\n')
[T0, B0] = RobustPCA_Doppler(M,Lambda); %
tRPCAEnd = toc(tRPCAStart)      % pair 2: toc
%T0 = T0/abs(max(T0(:)));

if 0
Mfinale=reshape(B0,Nz,Nx,Nt);
FigFeatures.title=1;
FigFeatures.result_folder = result_folder;
FigFeatures.mm=0;
FigFeatures.bar=1;
FigFeatures.print=0;
FigFeatures.nomtest = 'RPCA';
Dopplerplot(Mfinale,espace_xx,espace_zz,test,FigFeatures); 
clear Mfinale 
end
%%
fprintf('Running estimated PSF initiale....\n')
max_iter = 5;
Mt = reshape(M-T0,Nz,Nx,Nt);
M11 = squeeze(mean(Mt,3));
%M11=real(M11);
[H,psf0] = Hestimate(M11,Nz,Nx,Nt);
fprintf('Initialized PSF size: %d-%d\n',size(psf0,1),size(psf0,2))
clear Mt M11 

%%%%%%%%%%%%%%% Print PSF initiale 
FigHandlem=figure();
imagesc(psf0); colorbar;

export_fig(FigHandlem, ... % figure handle
    sprintf('%s/HEst_initialisation', result_folder),... % name of output file without extension
    '-painters', ...      % renderer
    '-transparent', ...   % renderer
    '-pdf', ...         % file format
    '-r500',... 
    '-nocrop' );             % resolution in dpi    

%%
FigFeatures.title=1;
FigFeatures.result_folder = result_folder;
FigFeatures.mm=0;
FigFeatures.bar=1;
FigFeatures.print=1;

% Stop condition
tol  = 1e-3;
xtmp = M;
Ttmp = T0;
err = zeros(1,max_iter);
normM = norm(M, 'fro');
for iter = 1:max_iter
    fprintf('Running estimated DRPCA for iteration %d....\n',iter)
    [T,x] =fastDRPCA(M, H, lambda, loops, rang0, tolS,[],[]);
    
    %% AFFICHAGE DE L'IMAGE DEROULANTE SELON Nt APRES SEUILLAGE/FILTRAGE
    Mfinale=reshape(x,Nz,Nx,Nt);
    %save(sprintf('%s/BDRPCA_simu_%d.mat', result_folder,iter),'Mfinale')   

    %% Doppler de puissance
    FigFeatures.nomtest = sprintf('B_image-Iter_%d',iter);
    Dopplerplot(Mfinale,espace_xx,espace_zz,test,FigFeatures); 
    %close
              
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

%         % Print H figures
         FigHandlem=figure();
         imagesc(abs(psf1)); colorbar; 
% 
%         % Print figure
%         export_fig(FigHandlem, ... % figure handle
%             sprintf('%s/HEst_%d', result_folder,iter+1),... % name of output file without extension
%             '-painters', ...      % renderer
%             '-transparent', ...   % renderer
%             '-pdf', ...         % file format
%             '-r500',... 
%             '-nocrop' );            % resolution in dpi    

    else 
        break;
    end
    
    clear Mt M11 psf1

end
tBDRPCAEnd = toc(tBDRPCAStart)      % pair 2: toc

%% AFFICHAGE DE L'IMAGE DEROULANTE SELON Nt APRES SEUILLAGE/FILTRAGE
Mfinale=reshape(x,Nz,Nx,Nt);
%save(sprintf('%s/fBDRPCA.mat', result_folder),'Mfinale')   
%%
%Doppler de puissance
FigFeatures.nomtest = 'fBDRPCA';
Dopplerplot(Mfinale,espace_xx,espace_zz,test,FigFeatures); 
clear Mfinale 

% if test ==1
%     fig2_3_fDRPCA
% else
%     fig6_fDRPCA_ISBI;
% end
