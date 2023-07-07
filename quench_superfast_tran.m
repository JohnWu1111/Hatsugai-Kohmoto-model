clear;
clc;
format long
tic;

%% parameter

L = 200;
% k = -pi/2 + 2*pi/L:2*pi/L:pi/2;
k = -1/2 + 2/L:2/L:1/2; % *pi
% k = 2/L:2/L:1;
U = 8;
VV0 = 4;
VV = 1;
E_k = -2*cospi(k');
nk = length(E_k);
step = 1000;
tol = 1e-5;

dt = 1e-3;
M = 1;
T_max = 100;
T = 0:M*dt:T_max;
nt = length(T);
nt_real = round(T_max / dt) + 1;

%% calculate GS

% only spin number = 2 need to be considered
% 4D-basis can be divided into 3+1D-basis
% 3D:(u'd+ud')/sqrt(2),ud,u'd'; 1D:(u'd-ud')/sqrt(2)

phi0_2 = zeros(3,nk);
phi0_2_4 = zeros(1,nk);
phi0_2(1,:) = 1/sqrt(2);
phi0_2(2,:) = 1/sqrt(2);
e2_4 = -U/2;

m0 = zeros(step,1);
m0(1) = sqrt(2)*((phi0_2(2,:)+phi0_2(3,:))*phi0_2(1,:)')/L;

for i = 2:step
    b = -2*sqrt(2)*m0(i-1)*VV0;
    for j = 1:nk
        a = 2*E_k(j);
        % spin number = 2
        H2 = [-U/2 b b;
            b -a+U/4 0;
            b 0 a+U/4];
        [V2,D2] = eig(H2);
        e2 = diag(D2);
        if e2(1) < e2_4
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

%% quench time evolution

phi_2 = phi0_2;
m = zeros(nt,1);
m(1) = m0(end);
m_it = m(1);

phi_21 = phi_2(1,:)';
phi_22 = phi_2(2,:)';
phi_23 = phi_2(3,:)';

phi_norm = zeros(nt,1);
phi_norm(1) = 2*sum(abs(phi_21).^2 + abs(phi_22).^2 + abs(phi_23).^2)/L;

count = 2;
for i = 2:nt_real
    d = -2*sqrt(2)*m_it*VV;
    % spin number = 2

    % solve the eigenvalue using pol3 root formula
    a = -U/2;
    b = -2*E_k + U/4;
    c = 2*E_k + U/4;
    ab = a*b;
    ac = a*c;
    bc = b.*c;
    bpc = b+c;
    alpha = a+b+c;
    beta = ab+bc+ac-2*d^2;
    gamma = a*bc - bpc*d^2;
    p = beta - alpha.^2/3;
    q = gamma - alpha.*beta/3 + 2*alpha.^2/27;
    temp = sqrt(-p/3);
    u = 2*temp;
    temp = 1./temp.^3;
    v = acos(-q/2.*temp)/3;
    angle = [0 2*pi/3 4*pi/3];
    root3 = u.*cos(v+angle);
    root3 = sort(root3,2,'descend');

    % constructing expH
    root32 = root3.^2;
    exproot3 = exp(1i*dt*root3);
    fact = beta + 2*alpha.*root3 + 3*root32;
    exproot3_fact = exproot3./fact;
    
    A1 = sum((bc + bpc.*root3 + root32).*exproot3_fact,2);
    B = -d*sum((c + root3).*exproot3_fact,2);
    C = -d*sum((b + root3).*exproot3_fact,2);
    A2 = sum((ac - d^2 + (a+c).*root3 + root32).*exproot3_fact,2);
    D = d^2*sum(exproot3_fact,2);
    A3 = sum((ab - d^2 + (a+b).*root3 + root32).*exproot3_fact,2);
    
    phi_21t = A1.*phi_21 + B.*phi_22 + C.*phi_23;
    phi_22t = B.*phi_21 + A2.*phi_22 + D.*phi_23;
    phi_23t = C.*phi_21 + D.*phi_22 + A3.*phi_23;

    % k = pi/2 is a special point
    phik_1 = phi_21(end);
    phik_2 = phi_22(end);
    phik_3 = phi_23(end);

    % rotate basis
    phik_2n = (phik_2 + phik_3)/sqrt(2);
    phik_3 = (phik_2 - phik_3)/sqrt(2);
    phik_2 = phik_2n;

    b = sqrt(2)*d;
    a = -3*U/8;
    % H = [a b;b -a];
    fact = sqrt(a^2+b^2);
    ft = fact*dt;
    ss = sin(ft);
    ss = ss/fact;
    cc = cos(ft);
    Es = a*ss;
    bs = b*ss;
    phik_1n = (cc-1i*Es)*phik_1 +1i*bs*phik_2;
    phik_2 = (cc+1i*Es)*phik_2 +1i*bs*phik_1;
    phik_1 = phik_1n;
    
    phik_2n = (phik_2 + phik_3)/sqrt(2);
    phik_3 = (phik_2 - phik_3)/sqrt(2);
    phik_2 = phik_2n;

    phi_21 = phi_21t;
    phi_22 = phi_22t;
    phi_23 = phi_23t;
    phi_21(end) = phik_1;
    phi_22(end) = phik_2;
    phi_23(end) = phik_3;

    m_it = sqrt(2)*real((phi_22+phi_23)'*(phi_21+phi0_2_4'))/L;
    if mod(i-1,M) == 0
        m(count) = m_it;
        phi_norm(count) = 2*sum(abs(phi_21).^2 + abs(phi_22).^2 + abs(phi_23).^2)/L;
        count = count + 1;
    end
end

figure;
plot(T,m)

toc;

function [phi_21t,phi_22t,phi_23t] = time_evo(phi_21,phi_22,phi_23,E_k,U,d,dt)
    a = -U/2;
    b = -2*E_k + U/4;
    c = 2*E_k + U/4;
    ab = a*b;
    ac = a*c;
    bc = b.*c;
    bpc = b+c;
    alpha = a+b+c;
    beta = ab+bc+ac-2*d^2;
    gamma = a*bc - bpc*d^2;
    p = beta - alpha.^2/3;
    q = gamma - alpha.*beta/3 + 2*alpha.^2/27;
    temp = sqrt(-p/3);
    u = 2*temp;
    temp = 1./temp.^3;
    v = acos(-q/2.*temp)/3;
    angle = [0 2*pi/3 4*pi/3];
    root3 = u.*cos(v+angle);
    root3 = sort(root3,2,'descend');

    % constructing expH
    root32 = root3.^2;
    exproot3 = exp(1i*dt*root3);
    fact = beta + 2*alpha.*root3 + 3*root32;
    exproot3_fact = exproot3./fact;
    
    A1 = sum((bc + bpc.*root3 + root32).*exproot3_fact,2);
    B = -d*sum((c + root3).*exproot3_fact,2);
    C = -d*sum((b + root3).*exproot3_fact,2);
    A2 = sum((ac - d^2 + (a+c).*root3 + root32).*exproot3_fact,2);
    D = d^2*sum(exproot3_fact,2);
    A3 = sum((ab - d^2 + (a+b).*root3 + root32).*exproot3_fact,2);
    
    phi_21t = A1.*phi_21 + B.*phi_22 + C.*phi_23;
    phi_22t = B.*phi_21 + A2.*phi_22 + D.*phi_23;
    phi_23t = C.*phi_21 + D.*phi_22 + A3.*phi_23;

%     i = 50;
%     H1 = [A1(i) B(i) C(i); B(i) A2(i) D(i); C(i) D(i) A3(i)];
%     H2 = [a d d; d b(i) 0; d 0 c(i)];
%     expH = expm(-1i*H2*dt);
% %     expH - H1
% 
%     phi_2 = [phi_21(i),phi_22(i),phi_23(i)]';
%     [V2,D2] = eig(H2);
%     e2 = diag(D2);
%     trans = V2'*phi_2;
%     phi_2it = V2*(exp(-1i*e2*dt).*trans);
%     phi_2it - conj([phi_21t(i),phi_22t(i),phi_23t(i)])'
%     exp(-1i*e2*dt) - exp(1i*dt*root3(i,:)')
%     e2 + root3(i,:)'
end

function [phik_1,phik_2,phik_3] = time_evo_k(phik_1,phik_2,phik_3,U,b,dt)
    % rotate basis
    phik_2n = (phik_2 + phik_3)/sqrt(2);
    phik_3 = (phik_2 - phik_3)/sqrt(2);
    phik_2 = phik_2n;

    a = -3*U/8;
    % H = [a b;b -a];
    fact = sqrt(a^2+b^2);
    ft = fact*dt;
    ss = sin(ft);
    ss = ss/fact;
    cc = cos(ft);
    Es = a*ss;
    bs = b*ss;
    phik_1n = (cc-1i*Es)*phik_1 +1i*bs*phik_2;
    phik_2 = (cc+1i*Es)*phik_2 +1i*bs*phik_1;
    phik_1 = phik_1n;
    
    phik_2n = (phik_2 + phik_3)/sqrt(2);
    phik_3 = (phik_2 - phik_3)/sqrt(2);
    phik_2 = phik_2n;
end