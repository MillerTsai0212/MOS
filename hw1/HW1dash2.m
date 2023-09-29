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
phi_B1 = -log(10^14 / (1.5*10^10)) * 0.0259;
phi_B2 = -log(10^16 / (1.5*10^10)) * 0.0259;
% Set x domain for na1 condition
x1_max = 0.56 + abs(phi_B1);
x1_min = -(0.56 - abs(phi_B1));
x1 = x1_min : x1_max / 1000 : x1_max;
% Formula reference for na2 condition
x2_max = 0.56 + abs(phi_B2);
x2_min = -(0.56 - abs(phi_B2));
x2 = x2_min : x2_max / 1000 : x2_max;
% Formula reference for na1 condition
phi_s1=x1+phi_B1;
us1=phi_s1/KT;
ub1=phi_B1/KT;
% Formula reference for na2 condition
phi_s2=x2+phi_B2;
us2=phi_s2/KT;
ub2=phi_B2/KT;
% Calculate na1 condition
na1 = 1e14;
lambda_p1 = (eps_si * KT / (2 * q * ni))^0.5;
Fs1= sign(ub1-us1) .* (2^0.5) * KT /  lambda_p1 .* ((ub1 - us1) * sinh(ub1) - (cosh(ub1) - cosh(us1))).^0.5;
Qs1 = eps_si * Fs1;
% Calculate na2 condition
na2 = 1e16;
lambda_p2 = (eps_si * KT / (2 * q * ni))^0.5;
Fs2 = sign(ub2-us2) .* (2^0.5) * KT / lambda_p2 .* ((ub2 - us2) * sinh(ub2) - (cosh(ub2) - cosh(us2))).^0.5;
Qs2 = eps_si * Fs2;
% Draw plot
semilogy(x1 , abs(Qs1) , x2 , abs(Qs2))
% Plot tag
legend({'Na = 1e14','Na = 1e16'},'Location','southeast')
xlabel('\psi_s');  
ylabel('|Qs| C/cm^2');
grid on