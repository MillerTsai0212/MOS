clear all
%% P-MOS
% VFB
Qeff = 2e11;
dox = 5e-6;
Na = 1e14;
Dit = 2.5e11;
ni = 1.5e10;
q = 1.6e-19;
KT = 0.0259;
eps0 = 8.85e-14;
eps_sio2 =  3.9 * eps0;
eps_si = 11.9 * eps0;
cox = eps_sio2 / dox;
phi_ms = 4.1 - (4.15 + 1.12 / 2 + KT * log(Na / ni));

% Qit
phi_B =  - KT * log(Na / ni);
shi_s =  abs(phi_B);
Qit =  - q * Dit * (phi_B + shi_s); 
Qit0 = - q * Dit * phi_B;

% Qs
phi_s = phi_B + shi_s;
us = phi_s / KT;
ub = phi_B / KT;
lambda_i = (eps_si * KT / (2 * q * ni))^0.5;
Fs= sign(ub-us) .* (2^0.5) * KT /  lambda_i .* ((ub - us) * sinh(ub) - (cosh(ub) - cosh(us))).^0.5;
Qs = eps_si * Fs;

% Vg
Vfb = phi_ms - Qeff * q / cox - Qit0 / cox;
Vg = Vfb - Qs / cox  - (Qit - Qit0) / cox + shi_s;

% CHF / Cox
wd = (2 *  abs(shi_s) * eps_si / (q * Na)) .^ 0.5;
Cd = eps_si / wd;
Chf = (1 / cox + 1 / Cd) .^-1;
Chf_over_Cox = Chf / cox;

% CLF / Cox
Cs = - sign(ub - us) * (eps_si / lambda_i) .* (sinh(us) - sinh(ub)) ./ ((2 ^ 0.5) * ((ub - us) * sinh(ub) - cosh(ub) + cosh(us))).^0.5;
Cit = q * Dit;
Clf = (1 / cox + 1 / (Cit + Cs)) .^-1;
Clf_over_Cox = Clf / cox;


