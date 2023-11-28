clear all
dox = 2.5e-6;
eps0 = 8.85e-14;
eps_si = 11.9 * eps0;
eps_sio2 = 3.9 * eps0;
Cm_over_Cox = 0.55;
Clf_over_Cox = 0.6;
Nw = 1e15;
q = 1.6e-19;

Cox = eps_sio2 / dox;
% dCoxCm__squre_dVg0 = -2 *  Cox^2 * (q * eps_si * 1e15)^-1;
% 
% dCoxCm__dVg0 = -2 *  Cox ^ 2 * 0.55 / (q * eps_si * 1e15);
d_dCm0Vg = -2 *  (q * eps_si * Nw)^-1;
dCoxCm0__squre_dVg0 = d_dCm0Vg * Cox^2;
dCoxCm0_dVg0 = dCoxCm0__squre_dVg0 * 0.55 * 0.5;

d_dCmVg = -2 * ((1 - Clf_over_Cox) / (1 - Cm_over_Cox)) * (q * eps_si * Nw)^-1;
dCoxCm__squre_dVg0 = d_dCmVg * Cox^2;
dCoxCm_dVg0 = dCoxCm__squre_dVg0 * 0.55 * 0.5;

