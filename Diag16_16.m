clear;
clc;
format long
tic;

L = 120;
% k = -pi/2 + 2*pi/L:2*pi/L:pi/2;
k = -1/2 + 2/L:2/L:1/2; % *pi
% k = 2/L:2/L:1;
U = 5;
VV = 0;
E_k = -2*cospi(k');
nk = length(E_k);

phi0_0 = zeros(1,nk);
phi0_1 = zeros(2,nk);
phi0_2 = zeros(4,nk);
phi0_3 = zeros(2,nk);
phi0_4 = zeros(1,nk);
phi0_2(1,:) = 1;

% m0 = 2*phi0_2(1,:)*phi0_2(2,:)'/L;
m0 = 0.5;

E0 = 0;
E1 = zeros(2,nk);
E2 = zeros(4,nk);
E2a = -U/2;
E3 = zeros(2,nk);
E4 = U/2;

b = -2*m0*VV;
for j = 1:nk
    a = E_k(j);

    % spin number = 1
    H1 = [a-U/4 b;
        b -a-U/4];
    [V1,D1] = eig(H1);
    e1 = diag(D1);
    E1(:,j) = e1;

    % spin number = 2
    H2 = [2*a+U/4 b b 0;
        b -U/2 0 b;
        b 0 -U/2 b;
        0 b b -2*a+U/4];
    [V2,D2] = eig(H2);
    e2 = diag(D2);
    E2(:,j) = e2;

    % spin number = 3
    H3 = [a b;
        b -a];
    [V3,D3] = eig(H3);
    e3 = diag(D3);
    E3(:,j) = e3;
end
%m0(i) = 2/L*(phi0_1(1,:)'*phi0_1(2,:));

E1_GS = E1(1,:);
E2_GS = min([E2(1,:);E2a*ones(1,nk)]);
E3_GS = E3(1,:);
E2_GS_sum = sum(E2_GS);
E13_GS_sum = sum(E1_GS(1:nk/4)) + sum(E1_GS(3*nk/4+1:nk)) + sum(E3_GS(nk/4+1:nk/2));
E31_GS_sum = sum(E3_GS(1:nk/4)) + sum(E3_GS(3*nk/4+1:nk)) + sum(E1_GS(nk/4+1:nk/2));

filename = strcat('L = ',num2str(L), ', U = ', num2str(U), ', V = ', num2str(VV), ', m0 = ', num2str(m0));
figure('Name',filename);
% set(gcf, 'position', [100 70 1700 900]);
hold on
plot(k,E0*ones(1,nk),'Color','black')
plot(k,E1_GS,'Color','b')
plot(k,E2_GS,'Color','r')
plot(k,E3_GS,'Color','g')
plot(k,E4*ones(1,nk),'Color','cyan')
legend('n=0','n=1','n=2','n=3','n=4')

toc;