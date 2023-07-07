clear;
clc;
format long
tic;

L = 200;
% k = -pi/2 + 2*pi/L:2*pi/L:pi/2;
k = -1/2 + 2/L:2/L:1/2; % *pi
% k = 2/L:2/L:1;
U_all = -3:-2;
VV = 1;
E_k = -2*cospi(k');
nk = length(E_k);
step = 100;
nU = length(U_all);

m_final = zeros(nU,1);

for n = 1:nU

    U = U_all(n);

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
    end
    m_final(n) = m0(end);
end

figure;
plot(U_all,m_final)
xlabel('U')
ylabel('CDW order parameter')

toc;