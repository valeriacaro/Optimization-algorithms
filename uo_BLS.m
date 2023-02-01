function [al,iWout] = uo_BLS(x,d,f,g,h,almax,almin,rho,c1,c2,iW)
    % First alpha value is the maximum one
    al = almax;

    % If ELS method is applied,
    if iW == 0
        al = -(g(x)'*d)/(d'*h(x)*d);
        % iWout = -1: ELS applied
        iWout = -1;
    % else if WC or SWC,
    else
        % Define Wolfe Condition 1
        WC1 = @(al) f(x+al*d) <= f(x) + c1*g(x)'*d*al;
        % Define Wolfe Condition 2
        WC2 = @(al) g(x+al*d)'*d >= c2*g(x)'*d;
        % Define Strong Wolfe Condition 2
        SWC2 = @(al) abs(g(x+al*d)'*d) <= c2*abs(g(x)'*d);
        % Define Descent Condition
        DC = @(al) g(x+al*d)'*d < 10^-6;
    end
    % Reduce alpha while it does not satisfy WC and it is bigger than alpha min
    while al >= almin && iW > 0
        if WC1(al)
            % iWout = 1: alpha satisfies (WC1)
            iWout = 1;
            if iW == 1 && WC2(al)
                % iWout = 2: alpha satisfies WC
                iWout = 2;
                break;
            elseif iW == 2 && SWC2(al) && DC(al)
                % iWout = 3: alpha satisfies SWC
                iWout = 3;
                break;
            end
        else
            % iWout = 0: alpha does not satisfy any WC
            iWout = 0;
        end
        al = al*rho;
    end
end