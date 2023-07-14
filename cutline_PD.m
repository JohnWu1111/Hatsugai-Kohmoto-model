clear;
clc;
format long
tic;

L = 1000;
% k = -pi/2 + 2*pi/L:2*pi/L:pi/2;
k = -1/2 + 2/L:2/L:1/2; % *pi
% k = 2/L:2/L:1;
U = 0;
VV_all = 0:0.01:5;
E_k = -2*cospi(k');
nk = length(E_k);
step = 1000;
tol = 1e-7;
nV = length(VV_all);

m_final = zeros(nV,1);
single_fill_final = zeros(nV,1);
step_store = zeros(nV,1);

parfor t = 1:nV

    VV = VV_all(t);

    phi0_2 = zeros(3,nk);
    % phi0_2_4 = zeros(1,nk);
    phi0_2(1,:) = 1/sqrt(2);
    phi0_2(2,:) = 1/sqrt(2);
    % e2_4 = -U/2;

    m0 = zeros(step,1);
    m0(1) = sqrt(2)*((phi0_2(2,:)+phi0_2(3,:))*phi0_2(1,:)')/L;

    for i = 2:step
        b = -2*sqrt(2)*m0(i-1)*VV;
        for j = 1:nk
            a = 2*E_k(j);
            % spin number = 2
            H2 = [-U/2 b b;
                b -a+U/2 0;
                b 0 a+U/2];
            [V2,D2] = eig(H2);
            e2 = diag(D2);
            phi0_2(:,j) = V2(:,1);
        end
        m0(i) = sqrt(2)*((phi0_2(2,:)+phi0_2(3,:))*phi0_2(1,:)')/L;
        if abs(m0(i) - m0(i-1)) < tol
            break;
        end
    end
    step_store(t) = i;
    m_final(t) = m0(i);
    single_fill_final(t) = 2*sum(abs(phi0_2(2,:)).^2 + abs(phi0_2(3,:)).^2)/L;
end

filename = strcat('L = ',num2str(L), ', U = ', num2str(U));
figure('Name',filename);
plot(VV_all,m_final);

toc;