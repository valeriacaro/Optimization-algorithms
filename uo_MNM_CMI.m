function [B, tau] = uo_MNM_CMI(h)
    % Compute lamba UB
    laUB = norm(h,'fro');
    % First tau is 0
    tau = 0; 
    % First l is 0
    l = 0;
    n = size(h,1);
    % Compute R and p
    [R,p] = chol(h+tau*eye(n));
    % If p greater than 0
    while p > 0
        % Compute tau again
        l = l + 1;
        tau = (1.01-1/(2^l))*laUB;
        % Recompute R and p
        [R,p] = chol(h+tau*eye(n));
    end
    % Compute B
    B = R'*R;
end

