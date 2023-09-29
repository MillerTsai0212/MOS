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

q_phi_m = 4.1;
q_x_si = 4.15;
eg_si = 1.12;
KT = 0.0259;
ni = 1.5 * 10^10;

if Nd > 0
    q_phi_s = q_x_si + 0.5 * eg_si - KT * log(Nd / ni);
    Vfb = q_phi_m - q_phi_s;
    fprintf('When Nd = %e, Vfb =  %f \n', Nd, Vfb); 
end

if Na > 0
    q_phi_s = q_x_si + 0.5 * eg_si + KT * log(Na / ni);
    Vfb = q_phi_m - q_phi_s;
    fprintf('When Na = %e, Vfb =  %f \n', Na, Vfb); 
end

