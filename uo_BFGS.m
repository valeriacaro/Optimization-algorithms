function [xk,dk,alk,iWk,betak,Hk,tauk] = uo_BFGS(x,f,g,h,epsG,kmax,almax,almin,rho,c1,c2,iW)
    % Create vectors and matrices
    xk = [x];
    dk = [];
    alk = [];
    iWk = [];
    betak = [];
    n = length(x);
    H = eye(n); % First H will be the identity
    Hk(:,:,1) = H;
    tauk = [];
    
    k = 0;
    % While a minimum has not been found and have not surpass kmax
    while norm(g(x)) > epsG && k < kmax
        % Define direction search
        d = -H*g(x);
        % Compute alpha
        [al,iWout] = uo_BLS(x,d,f,g,h,almax,almin,rho,c1,c2,iW);
        % Save actual x value
        xprev = x;
        % Compute new x
        x = x + al*d;
        % Save data
        dk = [dk,d];
        xk = [xk,x];
        iWk = [iWk, iWout];
        alk = [alk,al];
        % Compute and save new H
        sk = x - xprev;
        yk = g(x) - g(xprev);
        pk = 1/(yk'*sk);
        H = (eye(n)-(pk*sk*yk'))*H*(eye(n)-(pk*yk*sk'))+pk*sk*sk';
        Hk(:,:,k+1) = H;
        % Move to next iteration
        k = k+1;
    end
end 