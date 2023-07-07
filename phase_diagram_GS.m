clear;
clc;
format long
tic;

L = 500;
% k = -pi/2 + 2*pi/L:2*pi/L:pi/2;
k = -1/2 + 2/L:2/L:1/2; % *pi
% k = 2/L:2/L:1;
U_all = -10:0.1:20;
VV_all = 0:0.02:5;
E_k = -2*cospi(k');
nk = length(E_k);
step = 200;
tol = 1e-6;
nU = length(U_all);
nV = length(VV_all);

m_final = zeros(nU,nV);
step_store = zeros(nU,nV);

parfor n = 1:nU
    U = U_all(n);
    for t = 1:nV

        VV = VV_all(t);

        phi0_2 = zeros(4,nk);
        phi0_2(1,:) = 1/sqrt(2);
        phi0_2(2,:) = 1/sqrt(2);

        m0 = zeros(step,1);
        m0(1) = phi0_2(1,:)*phi0_2(2,:)'/L;

        for i = 2:step
            b = -2*m0(i-1)*VV;
            for j = 1:nk
                a = E_k(j);
                % spin number = 2
                H2 = [2*a+U/4 b b 0;
                    b -U/2 0 b;
                    b 0 -U/2 b;
                    0 b b -2*a+U/4];
                [V2,D2] = eig(H2);
                e2 = diag(D2);
                phi0_2(:,j) = V2(:,1);
            end
            m0(i) = ((phi0_2(1,:)+phi0_2(4,:))*(phi0_2(2,:)+phi0_2(3,:))')/L;
            if abs(m0(i) - m0(i-1)) < tol
                break;
            end
        end
        step_store(n,t) = i;
        m_final(n,t) = m0(i);
    end  
end

figure;
image(VV_all,U_all,m_final,'CDataMapping','scaled')
colorbar

toc;