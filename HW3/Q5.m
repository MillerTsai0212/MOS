clear all
dox = 2.5e-6;
q = 1.6e-19;
Qm_over_q = 6e12;
eps0 = 8.85e-14;
eps_si = 11.9 * eps0;
eps_sio2 = 3.9 * eps0;
KT = 0.0259;
alpha = 5;

Cox = eps_sio2 / dox;
Qm = Qm_over_q * q;
delQeff = Qm
delVFB_max = delQeff / Cox;
delVFB = 0.25 * delVFB_max;
dashX = dox * 0.25;
area = alpha * q * Qm_over_q;