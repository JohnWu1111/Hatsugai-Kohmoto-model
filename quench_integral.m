clear;
clc;
format long
tic;

%% parameter

L = 1000;
% k = -pi/2 + 2*pi/L:2*pi/L:pi/2;
k = -1/2 + 2/L:2/L:1/2; % *pi
% k = 2/L:2/L:1;
nk = length(k);

U = 4;
M0 = 10000;
dt = pi/M0;
M = M0;
T_max = 100*pi;
T = 0:dt:T_max-dt;
nt = length(T);

% integral
% m = sum(cos(4*cospi(k')*T))/L;
m = besselj(0,4*T)/2;
m = m.*cos(U*T);

dw = 2*pi/T_max;
w_max = pi/dt;
w = 0:dw:w_max;
m_f = abs(fft(m));
m_f = m_f(1:floor((nt+1)/2));

T_reshape = reshape(T,M,nt/M);
m_reshape = reshape(m,M,nt/M);
T_re = T_reshape(1,:);
% m_re = sum(m_reshape)*dt;
m_re = m_reshape(1,:);
nt_re = nt/M;

% phi0_1 = 1/sqrt(2);
% phi0_2 = 1/2;
% phi0_3 = 1/2;
% E1 = -2;
% E2 = 2*(1-2*cospi(k'));
% E3 = 2*(1+2*cospi(k'));
% 
% phi1t = phi0_1.*exp(-1i*E1.*T);
% phi2t = phi0_2.*exp(-1i*E2.*T);
% phi3t = phi0_3.*exp(-1i*E3.*T);
% 
% m = sqrt(2)*real(sum((phi2t+phi3t).*conj(phi1t)))/L;

% U = 4
% m_int_t_pi_4 = -T/2.*(besselj(0,4*T).*cos(4*T) + besselj(1,4*T).*sin(4*T))...
%     - (pi+4*T)/8.*(besselj(0,4*T+pi).*cos(4*T) + besselj(1,4*T+pi).*sin(4*T));

figure;
% plot(T(1000:end),m_int_t_pi_4(1000:end))
% plot(w,m_f)
% plot(log(T_re),log(m_re))
plot(T_re,m_re);

fitx = log(T_re(floor(nt_re/10):end));
fity = log(m_re(floor(nt_re/10):end));

toc;