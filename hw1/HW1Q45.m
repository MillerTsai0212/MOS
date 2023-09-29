clear all
type = input('NMOS enter 1, PMOS enter 2 ：');
if type == 1
    Nd = input('Nd = ');
    Na = -1;
end
if type == 2
    Na = input('Na = ');
    Nd = -1;
end
dox = 8e-6;
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
Cox = eps_sio2 / dox;

if Nd > 0
    phi_B = -log(Na / (1.5*10^10)) * 0.0259;
    shi_s = abs(phi_B);
    phi_s = phi_B + shi_s;
    q_phi_s = q_x_si + 0.5 * eg_si - KT * log(Nd / ni);
    Vfb = q_phi_m - q_phi_s;
    us=phi_s/KT;
    ub=phi_B/KT;
    lambda_i = (eps_si * KT / (2 * q * ni))^0.5;
    Fs= sign(ub-us).* (2^0.5) * KT /  lambda_i.* ((ub - us) * sinh(ub) - (cosh(ub) - cosh(us))).^0.5;
    Qs = eps_si * Fs;
    Vg = Vfb - Qs / Cox - shi_s;
    fprintf('When Nd = %e, Vg =  %f \n', Nd, Vg); 
end

if Na > 0
    phi_B = log(Na / (1.5*10^10)) * 0.0259;
    shi_s = - abs(phi_B);
    phi_s = phi_B + shi_s;
    q_phi_s = q_x_si + 0.5 * eg_si + KT * log(Na / ni);
    Vfb = q_phi_m - q_phi_s;
    us=phi_s/KT;
    ub=phi_B/KT;
    lambda_i = (eps_si * KT / (2 * q * ni))^0.5;
    Fs= sign(ub-us).* (2^0.5) * KT /  lambda_i.* ((ub - us) * sinh(ub) - (cosh(ub) - cosh(us))).^0.5;
    Qs = eps_si * Fs;
    Vg = Vfb - Qs / Cox - shi_s;
    fprintf('When Na = %e, Vg =  %f \n', Na, Vg); 
end






