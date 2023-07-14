clear;
clc;
format long
tic;

%% parameter

L = 200;
% k = -pi/2 + 2*pi/L:2*pi/L:pi/2;
k = -1/2 + 2/L:2/L:1/2; % *pi
% k = 2/L:2/L:1;
U = 0;
VV0 = 10;
VV = 3;
delta = 0.8;
E_k = -2*cospi(k');
nk = length(E_k);
step = 1000;
tol = 1e-5;
omega = 1;

dt = 1e-3;
M = 1;
T_max = 10;
T = 0:M*dt:T_max;
nt = length(T);
nt_real = round(T_max / dt) + 1;

%% calculate GS

phi0_2 = zeros(4,nk);
phi0_2(1,:) = 1/sqrt(2);
phi0_2(2,:) = 1/sqrt(2);

m0 = zeros(step,1);
m0(1) = phi0_2(1,:)*phi0_2(2,:)'/L;


for i = 2:step
    b = -2*m0(i-1)*VV0;
    for j = 1:nk
        a = E_k(j);
        % spin number = 2
        H2 = [2*a+U/2 b b 0;
            b -U/2 0 b;
            b 0 -U/2 b;
            0 b b -2*a+U/2];
        [V2,D2] = eig(H2);
        e2 = diag(D2);
        phi0_2(:,j) = V2(:,1);
    end
    m0(i) = ((phi0_2(1,:)+phi0_2(4,:))*(phi0_2(2,:)+phi0_2(3,:))')/L;
    if abs(m0(i) - m0(i-1)) < tol
        break;
    end
end

%% quench time evolution

phi_2 = phi0_2;
m = zeros(nt,1);
m(1) = m0(i);
m_it = m(1);

count = 2;
t_it = 0;
for i = 2:nt_real
    VV_it = VV + delta*cos(2*pi*omega*t_it);
    t_it = t_it + dt;
    b = -2*m_it*VV_it;
    for j = 1:nk
        a = E_k(j);
        % spin number = 2
        H2 = [2*a+U/2 b b 0;
            b -U/2 0 b;
            b 0 -U/2 b;
            0 b b -2*a+U/2];
        [V2,D2] = eig(H2);
        e2 = diag(D2);
        trans = V2'*phi_2(:,j);
        phi_2(:,j) = V2*(exp(-1i*e2*dt).*trans);
    end
    m_it = real((phi_2(1,:)+phi_2(4,:))*(phi_2(2,:)+phi_2(3,:))')/L;
    if mod(i-1,M) == 0
        m(count) = m_it;
        count = count + 1;
    end
end

figure;
plot(T,m)

toc;