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

cox = eps_sio2 / dox;
phi_ms = 4.1 - (4.15 + 1.12 / 2 + KT * log(Na / ni));
phi_B =  - KT * log(Na / ni);
Qit0 = - Dit * phi_B;
Vfb = phi_ms - Qeff * q / cox - Qit0 * q / cox;
% CFB / COX
lambda_p = (KT * 11.9 * eps0 / (q * Na)) .^0.5;
Cfbs = 11.9 * eps0 / lambda_p;
Cfb_over_Cox = Cfbs / (cox + Cfbs);
% CLF / COX
Cit = q * Dit;
Clf = (1 / cox +1 / (Cfbs + Cit)) .^-1;
Clf_over_Cox = Clf / cox;

%% N-MOS
Nd = 1e14;
cox = eps_sio2 / dox;
phi_msn = 4.1 - (4.15 + 1.12 / 2 - KT * log(Nd / ni));
phi_Bn =   KT * log(Na / ni);
Qitn0 =  - Dit * phi_Bn;
Vfbn = phi_msn - Qeff * q / cox - Qitn0 * q / cox;
% CFB / COX
lambda_n = (KT * 11.9 * eps0 / (q * Nd)) .^0.5;
Cfbsn = 11.9 * eps0 / lambda_n;
Cfbn_over_Cox = Cfbsn / (cox + Cfbsn);
% CLF / COX
Citn = q * Dit;
Clfn = (1 / cox +1 / (Cfbsn + Citn)) .^-1;
Clfn_over_Cox = Clfn / cox;
