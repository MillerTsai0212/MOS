clear all
ET = 0.3;
sigma_n = 2e-15;
Vth = 1e7;
g = 1;
Mc = 2.8e19;
KT = 0.0259;
t = 5e-8;

en = g * Vth * sigma_n * Mc * exp(-ET / KT);
tau_e = 1 / en;
ft = exp(-abs(en * t));