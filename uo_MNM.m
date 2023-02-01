function [xk, dk, alk, iWk, betak, Hk, tauk] = uo_MNM(x,f,g,h,epsG,kmax,almax,almin,rho,c1,c2,iW,delta,isd)
    % Create vectors and matrices
    xk = [x];
    dk = [];
    alk = [];
    iWk = [];
    betak = [];
    Hk = [];
    tauk = [];

    k = 1;
    % While a minimum has not been found and have not surpass kmax
    while k < kmax && norm(g(x)) > epsG
        % isd 5 --> Descomposició espectral (SD)
        if isd == 5
            % Compute h(x) vap's
            [Q,La] = eig(h(x));
            n = size(La,1);
            LaMod = La;
            for i = 1:n
                % Add to La the max between 0 and delta-vap
                LaMod(i,i) = max(0, delta-La(i,i)); 
            end
            % E is the matrice added to the hessian
            E = Q*LaMod*Q'
            % B is the modified hessian
            B = Q*(La+LaMod)*Q';
            % Save data
            Hk(:,:,k) = B;
            tauk = [tauk,0];
        % isd 6 --> Factorització de Cholesky (CMI)
        elseif isd == 6
            % Compute B and tau
            [B, tau] = uo_MNM_CMI(h(x));
            % Save data
            Hk(:,:,k) = B;
            tauk = [tauk,tau];
        end
        % Compute d
        d = -inv(B)*g(x);
        % Compute alpha using BLS function
        [al, iWout] = uo_BLS(x,d,f,g, h,almax,almin,rho,c1,c2,iW);
        % Compute new x
        x = x + al*d; 
        % Save data
        xk = [xk, x];
        dk = [dk,d];
        alk = [alk, al]; 
        iWk = [iWk, iWout];
        betak = [betak,0];
        % Move to next iteration
        k = k + 1;
    end
end

