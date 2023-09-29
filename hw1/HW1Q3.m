clear all
% Parameter
eps_si = 11.9 * 8.85 * 10^-14;
eps_sio2 = 3.9 * 8.85 * 10^-14;
ni = 1.5 * 10^10;
K = 1.38 * 10^-23;
T = 300;
KT = 0.0259;
q = 1.6 * 10^-19;
q_phi_m = 4.1;
q_x_si = 4.15;
eg_si = 1.12;
phi_B = -log(10^14 / ni) * 0.0259;
na = 1e14;
% Calculate Vfb
q_phi_s = q_x_si + eg_si/2 + KT * log(na / ni);
Vfb = q_phi_m - q_phi_s;
% Set x domain
x_max = 0.56 + abs(phi_B);
x_min = - (0.56 - abs(phi_B));
x = x_min : x_max / 1000 : x_max;
% Formula reference
phi_s=x+phi_B;
us=phi_s/KT;
ub=phi_B/KT;
% Calculate dox1 condition
dox1 = 5e-6;
Cox1 = eps_sio2 / dox1;
lambda_p1 = (eps_si * KT / (2 * q * ni))^0.5;
Fs1= sign(ub-us).* (2^0.5) * KT /  lambda_p1 .* ((ub - us) * sinh(ub) - (cosh(ub) - cosh(us))).^0.5;
Qs1 = eps_si * Fs1;
Vg1 = Vfb - Qs1 / Cox1 + x;
% Calculate dox2 condition
dox2 = 8e-6;
Cox2 = eps_sio2 / dox2;
lambda_p1 = (eps_si * KT / (2 * q * ni))^0.5;
Fs2 =sign(ub-us).* (2^0.5) * KT / lambda_p1 .* ((ub - us) * sinh(ub) - (cosh(ub) - cosh(us))).^0.5;
Qs2 = eps_si * Fs2;
Vg2 = Vfb - Qs2 / Cox2 + x;
% Draw plot
hold on
 plot(x , Vg1)
 plot(x , Vg2)
hold off
legend({'dox = 500 Å','dox = 800 Å'},'Location','southwest')
xlabel('\psi_s');  
ylabel('Vg');
grid on