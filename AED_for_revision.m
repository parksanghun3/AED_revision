%%%%%%%%% AED ( CDR MASK �� Revision ����) %%%%%%%%%%%

clear all; close all ; clc;

M = 512; % length of impulse response vectors
% mu = 0.003;
mu = 1;
% s = audioread('m60_68_FKFB0_con.wav');
% s = audioread('m45_68_MDAC0_con.wav');
s = audioread('x_1x2.wav');
x1 = s(:,1);
x2 = s(:,2);


% �ʱ�ȭ
g2 = [1 zeros(1,M/2-1) zeros(1,M/2)];
g1 = [0 zeros(1,M/2-1) zeros(1,M/2)];

% u ������Ʈ

len = floor((length(x1)-M)/M);

for n = 1 : 5000000

    if (mod(n,len))~=0
        i = mod(n,len);
    else
        i = len;
    end
    
    x_temp(1,:) = x1((i*M) : -1 : (i-1)*M+1);
    x_temp(2,:) = x2((i*M) : -1 : (i-1)*M+1);
    x = [x_temp(1,:) x_temp(2,:)]'; % ������ x1*g2 = x2*g1 �� �����ϱ� ����, g�� �״�� �ΰ� x�� �¿츦 change

    u = [g2 -g1]';
    
    e = u' * x;
    u_temp = u - mu*e*x;
    u = u_temp / norm(u_temp,2);
    
    g2 = u(1 : M)';
    g1 = u(M+1 : 2*M)';
    
end


 u = [g2 -g1]';
 
[a b]=max(u(1:M));
[c d]=max(u(M+1:2*M));
d-b

figure(1)
subplot(211)
stem(u(1:M));
xlabel('k');ylabel('g_1(k)');
subplot(212)
stem(u(M+1:2*M));
xlabel('k');ylabel('g_2(k)');