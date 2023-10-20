clear all
q = 1.6e-19;
Qeff = 2e11;
dox = input('dox(cm) = ');
Na = 1e14;
ni = 1.5e10;
KT = 0.0259;
eps0 = 8.85e-14;
Dit = 2.5e11;


eps_sio2 =  3.9 * eps0;
eps_si = 11.9 * eps0;
Cox = eps_sio2 / dox;
phi_B = -KT * log(Na / ni);
phi_ms = 4.1 - (4.15 + 1.12 / 2 + KT * log(Na / ni));
lambda_i = (eps_si * KT / (2 * q * ni))^0.5;
ub = phi_B / KT;


for VG = -5:0.01:2;
    %% CLF
    us = @(psi_s)  (psi_s + phi_B)/KT;
    % Qit
    Qit = @(psi_s) - (psi_s + phi_B) * Dit;
    % Qs
    Qs= @(psi_s) eps_si * sign(ub - us(psi_s)) .* (2^0.5) * KT /  lambda_i .* ((ub - us(psi_s)) * sinh(ub) - (cosh(ub) - cosh(us(psi_s)))).^0.5;
    % phi_s (帶入phi_s後算出其他參數再用Vg公式使其檢調 Vg後 = 0 反向得出當下phi_s)
    psi_s= fzero( @(psi_s) phi_ms - Qeff * q / Cox - Qit(psi_s) * q / Cox - Qs(psi_s) / Cox + psi_s - VG, 0 );
    % Cs
    Cs = - sign(ub - us(psi_s)) * (eps_si / lambda_i) .* (sinh(us(psi_s)) - sinh(ub)) ./ ((2 ^ 0.5) * ((ub - us(psi_s)) * sinh(ub) - (cosh(ub) - cosh(us(psi_s))))).^0.5;
    % Cit
    Cit = q * Dit;
    % CLF
    CLF = (1 / Cox + 1 / (Cit + Cs)) .^-1;
    % CLF ratio
    CLF_over_Cox = CLF / Cox;

    figure = plot(VG, CLF_over_Cox, 'r.');
    xlabel('Vg(V)');  
    ylabel('C / C_{ox}') 
    hold on 

    %% CHF
    if psi_s < abs(phi_B)
        
        % CHF
        CHF = (Cox * Cs) / (Cox + Cs);
        % CHF ratio
        CHF_over_Cox = CHF / Cox;
        figure = plot(VG, CHF_over_Cox, 'b.');
        xlabel('Vg(V)');  
        ylabel('C / C_{ox}') 
        hold on 

    else

        Wd = (2 *  psi_s * eps_si / (q * Na)) .^ 0.5;
        Cd = eps_si / Wd;
        CHF = (Cox * Cd) / (Cox + Cd);
        CHF_over_Cox = CHF / Cox;
        figure = plot(VG, CHF_over_Cox, 'b.');
        xlabel('Vg(V)');  
        ylabel('C / C_{ox}') 
        hold on 

   
    
    hold on 
    end
end
