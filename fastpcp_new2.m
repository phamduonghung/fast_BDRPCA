function [L1, x, statsPCP] = fastpcp_new2(V, H, lambda, loops, rang0,tolS, rangThreshold, lambdaFactor,Nz,Nx,Nt,espace_xx,espace_zz)
    %% Testing with new fRPCA
    
    % Data size
    %[Nrows,Ncols] = size(V);
    % V = M;
    [m, n] = size(V);
    unobserved = isnan(V);
    V(unobserved) = 0;
    normV = norm(V, 'fro');
    
    ll=7;
    if( nargin < ll)
        lambdaFactor = 1.0;
        if( nargin < ll-1)
            rangThreshold = 0.01;  
            if( nargin < ll-2)
                rang0 = 1;                % initial rang
                if( nargin < ll-3)
                    loops = 2;              % total number of outer loops
                    if( nargin < ll-4)
                        lambda = 1./sqrt(max(m,n));
                        if( nargin < ll-5)
                            H = ones(m,n);
                        end
                    end
                end
            end
        end
    end

    if isempty(rangThreshold)
        rangThreshold = 0.01;
    end

    if isempty(lambdaFactor)
        lambdaFactor = 1.0;
    end    
    
    % Set flag (increments the rang plus one at each iteration)
    inc_rang = 1;

    % ---------------------------------
    % >>>  measure time performance <<<
    % ---------------------------------
    %t = tic;
    
    % ------------------------    
    % --- First outer loop ---
    rang = rang0;                     % current rang
    statsPCP.rang(1) = rang;          % save current rang
    
    % ------------------------    
    % initial solution
    x = zeros(m, n);
    W = zeros(m, n); % W
    %Hx = real(ifft2(H.*fft2(x)));
    mu = 1e0;
    if 0
        % Partial SVD
        %[Ulan Slan Vlan] = lansvd(V, rang, 'L');
        %[Ulan,Slan,Vlan] = svds(V, rang);
        [Ulan,Slan,Vlan] = svdsecon(V, rang); % fastest

        % Current low-rang approximation
        L1 = Ulan*Slan*Vlan';

        % Shrinkage
        z = So(lambda/mu, x + (1/mu)*W);
        x1 = real(ifft2(conj(H).*fft2(V-L1))) + z + (1/mu)*(- W);
        h1 = 1./(abs(H).^2 + ones(m,n));
        x = real(ifft2(h1.*fft2(x1)));
        %Hx = real(ifft2(H.*fft2(x)));
        Z2 = x - z;
        Z2(unobserved) = 0; % skip missing values
        W = W + mu*Z2;
        %S1 = shrink(V-L1, lambda);
        % figure
        FigFeatures.title=0;
        FigFeatures.result_folder = 'C:\Users\dpham\ownCloud\Working\Atempo\simulation\RAPID_test\Results_LSDe_BDeconv\simu_conv';
        FigFeatures.mm=0;
        FigFeatures.bar=1;
        FigFeatures.print=0;
        FigFeatures.nomtest = 'Bimage_FPCP';
        Mfinale=reshape(x,Nz,Nx,Nt);
        Dopplerplot(Mfinale,espace_xx,espace_zz,1,FigFeatures);     
    end
    % ------------------------
    % ---    Outer loops   ---
sumLoop = 0;
err2 = zeros(1,loops);
for k = 1:loops
    timeLoop = tic;
    if(inc_rang == 1)
        lambda = lambda * lambdaFactor;         % modify Lambda at each iteration
        rang = rang + 0;                        % increase rang
    end

    % low rang (partial SVD)
    %[Ulan Slan Vlan] = lansvd(V-S1, rang, 'L');
    %[Ulan,Slan,Vlan] = svds(V-S1, rang);
    [Ulan,Slan,Vlan] = svdsecon(V-x, rang); % fastest

    currentEvals = diag(Slan);                                    % extract current evals
    statsPCP.rang(k) = length( currentEvals );                    % save current rang
    statsPCP.rho(k) = currentEvals(end) / sum( currentEvals(1:end-1) );       % relative contribution of the last evec

    % simple rule to keep or increase the current rang's value
    if(statsPCP.rho(k) < rangThreshold ) 
        inc_rang = 0;
    else
        inc_rang = 1;
    end

    % Current low-rang approximation
    L1 = Ulan*Slan*Vlan';
    
    if 0
        % Shrinkage
        %S1 = shrink(V-L1, lambda);
        z = So(lambda/mu, x + (1/mu)*W);
        x1 = real(ifft2(conj(H).*fft2(V-L1))) + mu*z - W;
        h1 = 1./(abs(H).^2 + mu*ones(m,n));
    else
        % Shrinkage
        z = So(lambda/mu, x + (1/mu)*W);
        x1 = real(ifft2(conj(H).*fft2(V-L1))) + z + (1/mu)*(- W);
        h1 = 1./(abs(H).^2 + ones(m,n));
    end
    x = real(ifft2(h1.*fft2(x1)));
    %Hx = real(ifft2(H.*fft2(x)));
    Z2 = x - z;
    Z2(unobserved) = 0; % skip missing values
    W = W + mu*Z2;
    
    % figure
    FigFeatures.title=0;
%    FigFeatures.result_folder = result_folder;
    FigFeatures.mm=0;
    FigFeatures.bar=1;
    FigFeatures.print=0;
    FigFeatures.nomtest = 'Bimage_FPCP';
    Mfinale=reshape(x,Nz,Nx,Nt);
    Dopplerplot(Mfinale,espace_xx,espace_zz,1,FigFeatures);   
    
    err2(k) = norm(Z2, 'fro') / normV;
    TL = toc(timeLoop); 
    sumLoop = sumLoop + TL;
    if (err2(k) > tolS) 
        fprintf(1, 'iter: %f\terr2: %f\trang(L1): %d\tcard(S): %d\t times: %.2fs\n', ...
                k, err2(k), rank(L1), nnz(x(~unobserved)),sumLoop);
    else
        break; 
    end
    if k>=2 && err2(k)>err2(k-1)
        break;
    end
    
    pause(0.1)
    close
end

% ---------------------------------
% >>>  measure time performance <<<
% ---------------------------------
%statsPCP.time = toc(t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function u = shrink(v, lambda)
%  u = sign(v).*max(0, abs(v) - lambda);
%end 

% function r = So(tau, S)
%     % shrinkage operator
%     r = sign(S) .* max(abs(S) - tau, 0);
% end