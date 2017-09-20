function rad_hex(T)

% this solves the problem of internal radiation for a cube

F = (ones(6,6) - eye(6,6)) * 0.2;
A = 5.67e-8 * (eye(6,6) - F);
B = 1/0.8 * eye(6,6) - .2/.8 * eye(6,6) * F;
ABinv = A * inv(B);

coupling_mat = .25 * ...
[ 0  1  1  0  0  1  1  0;...
  0  0  1  1  0  0  1  1;...
  1  0  0  1  1  0  0  1;...
  1  1  0  0  1  1  0  0;...
  1  1  1  1  0  0  0  0;...
  0  0  0  0  1  1  1  1]';


q_fe = coupling_mat * ABinv * (coupling_mat' * (T+273.16)).^4;
T_jac = (coupling_mat' * (T+273.16)).^3;
T_jac_mat = zeros(6,6);
for i=1:6,
    T_jac_mat(i,i) = T_jac(i,1);
end
jac = coupling_mat * ABinv * T_jac_mat;

T_jac_mat
coupling_mat
CABinv = coupling_mat * ABinv
q_fe
jac