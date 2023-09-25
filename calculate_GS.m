 clear;
clc;
format long
tic;

L = 1000;
% k = -pi/2 + 2*pi/L:2*pi/L:pi/2;
k = -1/2 + 2/L:2/L:1/2; % *pi
% k = 2/L:2/L:1;
U = 4;
VV = 0.5;
E_k = -2*cospi(k');
nk = length(E_k);
step = 1000;
tol = 1e-5;

phi0_2 = zeros(3,nk);
phi0_2_4 = zeros(1,nk);
phi0_2(1,:) = 1/sqrt(2);
phi0_2(2,:) = 1/sqrt(2);
e2_4 = -U/2;

m0 = zeros(step,1);
m0(1) = sqrt(2)*((phi0_2(2,:)+phi0_2(3,:))*phi0_2(1,:)')/L;

for i = 2:step
    b = -2*sqrt(2)*m0(i-1)*VV;
    for j = 1:nk
        a = 2*E_k(j);
        % spin number = 2
        H2 = [-U/2 b b;
            b a+U/2 0;
            b 0 -a+U/2];
        [V2,D2] = eig(H2);
        e2 = diag(D2);
        if e2(1) <= e2_4
            phi0_2(:,j) = V2(:,1);
            phi0_2_4(j) = 0;
        else
            phi0_2(:,j) = zeros(3,1);
            phi0_2_4(j) = 1;
        end
    end
    m0(i) = sqrt(2)*((phi0_2(2,:)+phi0_2(3,:))*(phi0_2(1,:)+phi0_2_4)')/L;
    if abs(m0(i) - m0(i-1)) < tol
        break;
    end
end
m_final = m0(end);

phi_21 = phi0_2(1,:)';
phi_22 = phi0_2(2,:)';
phi_23 = phi0_2(3,:)';

k_space_singlon = abs(phi_21.^2);
k_space_doublonk = abs(phi_22.^2);
k_space_doublonkp = abs(phi_23.^2);

k_space0(1,:) = k_space_singlon;
k_space0(2,:) = k_space_doublonk;
k_space0(3,:) = k_space_doublonkp;

filename = strcat('L = ',num2str(L), ', U = ', num2str(U), ', V = ', num2str(VV));
figure('Name',filename);
plot(k,k_space0)
legend('singlon','doublon_k','doublon_k''')

toc;