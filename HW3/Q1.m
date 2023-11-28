clear all
%% P-MOS
% VFB
Qeff = 2e11;
dox = 2.5e-6;
Na = 4e14;
Dit = 2.5e11;
ni = 1.5e10;
q = 1.6e-19;
KT = 0.0259;
eps0 = 8.85e-14;
eps_sio2 =  3.9 * eps0;
eps_si = 11.9 * eps0;

phi_b =  - KT * log(Na / ni);
cox = eps_sio2 / dox;
Chf_over_Cox = 0.55;
Clf_over_Cox = 0.6;
Cd = cox * Chf_over_Cox / (1 - 0.55);
Wd = eps_si / Cd;
phi_s = (Wd ^ 2 * q * Na) / (2 * eps_si);
cit = (0.6 * cox - 0.4 * Cd) / 0.4;
dit = cit / q;
Eu_minus_Ei = phi_s+phi_b;