function [xk, dk, alk, iWk, betak, Hk, tauk] = uo_NM(x,g,h,epsG,kmax)
    % Create vectors and matrices
    xk = [x];
    dk = [];
    alk = [];
    iWk = [];
    betak = [];
    Hk = [];
    tauk = [];
    tauk = [];
   
    k = 1;
    % While a minimum has not been found and have not surpass kmax
    while k < kmax && norm(g(x)) > epsG
        % Compute p = alpha*d
        p = -inv(h(x))*g(x);
        % Compute new x
        x = x + p;
        % Save data
        xk = [xk, x];
        dk = [dk,p];
        iWk = [iWk,4];
        alk = [alk,1];
        tauk = [tauk,0];
        betak = [betak,0];
        Hk = [Hk, 0];
        % Move to next iteration
        k = k + 1;
    end
end

