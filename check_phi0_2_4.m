clear;
clc;
format long
tic;

a = 0:0.01:0.5;
na = length(a);
U = 1;
E_k = -1;
e2_store = zeros(3,na);
for i = 1:na
    b = a(i);
    H2 = [-U/2 b b;
        b E_k+U/2 0;
        b 0 -E_k+U/2];
    e2 = eig(H2);
    e2_store(:,i) = e2;
end

figure;
plot(a,e2_store)

toc;