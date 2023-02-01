function [xk, dk, alk, iWk, betak, Hk, tauk] = uo_CGM(x,f,g,h,epsG,kmax,almax,almin,rho,c1,c2,iW,icg,irc,nu)
    % Create vectors and matrices
    xk = [x];
    d = -g(x); % First direction search will be equal to the one in GM
    dk = [d];
    alk = [];
    iWk = [];
    betak = [0]; % First beta value will be 0 bc we are using GM
    Hk = [];
    tauk = [];
    k = 1;
    n = size(x); n = n(1);
    % If using Fletcher-Reeves,
    if icg == 1
        beta = @(g,x,xprev) (g(x)'*g(x))/norm(g(xprev))^2; 
    % Else using Polak-Riviere
    elseif icg == 2
        beta = @(g,x,xprev) max(0,g(x)'*(g(x)-g(xprev))/norm(g(xprev))^2);
    end
    % While a minimum has not been found and have not surpass kmax
    while norm(g(x)) > epsG && k < kmax
        % Find alpha using BLS
        [alpha, iWout] = uo_BLS(x,d,f,g,h,almax,almin,rho,c1,c2,iW);
        % Save actual x value
        xprev = x;
        % Update x value
        x = xprev + alpha*d; 
        % Compute beta
        bet = beta(g,x,xprev); 
        % Save values in vectors
        xk = [xk,x];
        alk = [alk, alpha]; 
        iWk = [iWk, iWout];
        betak = [betak,bet];
        % Move to next iteration
        k = k + 1;
        % If Restart Condition 1 takes place
        if irc == 1 && mod(k-1,n) == 0
            % Compute d again
            d = -g(x);
            % Save its value
            dk = [dk,d];
        % If Restart Condition 2 takes place
        elseif irc == 2 && abs(g(x)'*g(xprev)) >= (norm(g(x))^2)*nu
            % Compute d again
            d = -g(x); 
            % Save its value
            dk = [dk,d];
        % If there is not a restart condition
        else
            % Compute d using beta
            d = -g(x) + bet*d; 
            % Save its value
            dk = [dk,d];
        end
        Hk = [Hk,0];
        tauk = [tauk,0];
    end
end