clear;
clc;
format long
tic;

%% parameter

L = 200;
% k = -pi/2 + 2*pi/L:2*pi/L:pi/2;
k = -1/2 + 2/L:2/L:1/2; % *pi
% k = 2/L:2/L:1;
U = 2;
VV0 = 10;
VV = 3;
delta = 0.8;
E_k = -2*cospi(k');
nk = length(E_k);
step = 1000;
tol = 1e-5;
omega = 1;

dt = 1e-4;
M = 100;
T_max = 100;
T = 0:M*dt:T_max;
nt = length(T);
nt_real = round(T_max / dt) + 1;

%% calculate GS

% only spin number = 2 need to be considered
% 4D-basis can be reduced to 3+1D-basis
% 3D:(u'd+ud')/sqrt(2),ud,u'd'; 1D:(u'd-ud')/sqrt(2)

phi0_2 = zeros(3,nk);
% phi0_2_4 = zeros(1,nk);
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
            b a+U/2 0;
            b 0 -a+U/2];
        [V2,D2] = eig(H2);
        e2 = diag(D2);
        phi0_2(:,j) = V2(:,1);
        % phi0_2_4(j) = 0;
    end
    m0(i) = sqrt(2)*((phi0_2(2,:)+phi0_2(3,:))*phi0_2(1,:)')/L;
    if abs(m0(i) - m0(i-1)) < tol
        break;
    end
end

%% quench time evolution

phi_2 = phi0_2;
m = zeros(nt,1);
m(1) = m0(i);
m_it = m(1);

phi_21 = phi_2(1,:)';
phi_22 = phi_2(2,:)';
phi_23 = phi_2(3,:)';

% phi_norm = zeros(nt,1);
% phi_norm(1) = 2*sum(abs(phi_21).^2 + abs(phi_22).^2 + abs(phi_23).^2)/L;

k_space_singlon = zeros(nt,nk);
k_space_doublonk = zeros(nt,nk);
k_space_doublonkp = zeros(nt,nk);
k_space_singlon(1,:) = abs(phi_21.^2);
k_space_doublonk(1,:) = abs(phi_22.^2);
k_space_doublonkp(1,:) = abs(phi_23.^2);

count = 2;
t_it = 0;
for i = 2:nt_real
    VV_it = VV + delta*cos(2*pi*omega*t_it);
    t_it = t_it + dt;
    d = -2*sqrt(2)*m_it*VV_it;
    % spin number = 2

    % solve the eigenvalue using pol3 root formula
    a = -U/2;
    b = 2*E_k + U/2;
    c = -2*E_k + U/2;
    ab = a*b;
    ac = a*c;
    bc = b.*c;
    bpc = b+c;
    alpha = a+b+c;
    beta = ab+bc+ac-2*d^2;
    gamma = a*bc - bpc*d^2;
    AA = alpha.^2 - 3*beta;
    BB = alpha.*beta - 9*gamma;
    % CC = beta.^2 - 3*alpha.*gamma;
    temp = sqrt(AA);
    TT = (AA.*alpha - 3*BB/2)./temp.^3;
    theta = acos(TT)/3;
    angle = [0 2*pi/3 4*pi/3];
    root3 = (-alpha-2*temp.*cos(theta + angle))/3;
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
    % a = -U/2;
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

    m_it = sqrt(2)*real((phi_22+phi_23)'*phi_21)/L;
    if mod(i-1,M) == 0
        m(count) = m_it;
        k_space_singlon(count,:) = abs(phi_21.^2);
        k_space_doublonk(count,:) = abs(phi_22.^2);
        k_space_doublonkp(count,:) = abs(phi_23.^2);
%         phi_norm(count) = 2*sum(abs(phi_21).^2 + abs(phi_22).^2 + abs(phi_23).^2)/L;
        count = count + 1;
    end
end

cut = 800;
phi_f = abs(fft(m(floor(nt/2)+1:end),(nt+1)/2));
phi_f(1) = 0;
phi_f_main = phi_f(1:cut);
dw = 1/(T_max/2);
w_max = 1/(M*dt);
w = 0:dw:w_max;
w_main = w(1:cut);
[maxtab, ~]=peakdet(phi_f_main, max(phi_f_main)/3);
if ~isempty(maxtab)
    [~, peak] = max(maxtab(2));
    peak_f = w(maxtab(peak,1));
    peak_per = maxtab(peak,2)/sum(maxtab(2));
else
    peak_f = 0;
    peak_per = 0;
    warning('no peak in the inteval!')
end

k_space_all(1,:) = mean(k_space_singlon(floor(nt*0.9):end,:));
k_space_all(2,:) = mean(k_space_doublonk(floor(nt*0.9):end,:));
k_space_all(3,:) = mean(k_space_doublonkp(floor(nt*0.9):end,:));

k_space0(1,:) = k_space_singlon(1,:);
k_space0(2,:) = k_space_doublonk(1,:);
k_space0(3,:) = k_space_doublonkp(1,:);


filename = strcat('L = ',num2str(L), ', U = ', num2str(U), ', Vi = ', num2str(VV0), ', Vf = ', num2str(VV), ', delta = ', num2str(delta));
figure('Name',filename);
% plot(T(floor(nt*0.95):end),m(floor(nt*0.95):end))
% plot(w_main,phi_f_main)
% plot(T,m)
% image(k,T,k_space_doublonkp,'CDataMapping','scaled')
% imagesc(k,T,k_space_doublonkp)
plot(k,k_space0)

toc;