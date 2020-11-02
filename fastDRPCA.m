function [L1, x, statsPCP] = fastDRPCA(V, H, lambda, loops, rang0,tolS, rangThreshold, lambdaFactor)
    % fast RRPCA created by PHAM Duong Hung
    % Email: duong-hung.pham@irit.fr
    % Version 02/11/2020 
    %%    
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
    
    if isempty(tolS)
        tolS = 1e-3;
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
    mu = 1e0;  
    % ------------------------
    % ---    Outer loops   ---
sumLoop = 0;
err2 = zeros(1,loops);
for k = 1:loops
    %k
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
    

    %timeLoops = tic;
    % LASSO L1
    z = So(lambda/mu, x + (1/mu)*W);
    x1 = real(ifft2(conj(H).*fft2(V-L1))) + z + (1/mu)*(- W);
    h1 = 1./(abs(H).^2 + ones(m,n));
    %TLs = toc(timeLoops)
   
    x = real(ifft2(h1.*fft2(x1)));
    Z2 = x - z;
    Z2(unobserved) = 0; % skip missing values
    W = W + mu*Z2;    
    err2(k) = norm(Z2, 'fro') / normV;
    TL = toc(timeLoop); 
    sumLoop = sumLoop + TL;
    if (err2(k) > tolS) 
        fprintf(1, 'iter: %d\terr2: %f\trang(L1): %d\tcard(S): %d\t times: %.2fs\n', ...
                k, err2(k), rank(L1), nnz(x(~unobserved)),sumLoop);
    else
        break; 
    end
    if k>=2 && err2(k)>err2(k-1)
        break;
    end
end
L1 = V-real(ifft2(H.*fft2(x)));




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