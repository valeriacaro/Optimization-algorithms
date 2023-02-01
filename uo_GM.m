function [xk, dk, alk, iWk, betak, Hk, tauk] = uo_GM(x,f,g,h,epsG,kmax,almax,almin,rho,c1,c2,iW)
    % Create vectors and matrices
    xk = [x];
    dk = [];
    alk = [];
    iWk = [];
    betak = [];
    Hk = [];
    tauk = [];
    % While a minimum has not been found and have not surpass kmax
    k = 0;
    while norm(g(x)) > epsG && k < kmax 
        % Direction search is -g(x) bc GM
        d = -g(x); 
        % Find alpha using BLS
        [al, iWout] = uo_BLS(x,d,f,g,h,almax,almin,rho,c1,c2,iW);
        % Update x
        x = x + al*d;
        % Let's do another iteration
        k = k+1;
        % Save values to their respective vectors
        xk = [xk, x];
        dk = [dk,d];
        alk = [alk, al]; 
        iWk = [iWk, iWout];
        betak = [betak,0];
        Hk = [Hk,0];
        tauk = [tauk,0];
    end
end

