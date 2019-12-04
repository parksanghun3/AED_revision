%%%%%%%%% AED ( CDR MASK 논문 Revision 목적) %%%%%%%%%%%

clear all; close all ; clc;

C= 343;
micdist = 0.12;
M = 512; % length of impulse response vectors
% mu = 0.003;
mu = 0.001;
% s = audioread('m60_68_FKFB0_con.wav');
% s = audioread('m45_68_MDAC0_con.wav');

RT60 = 0.6; %0.2 : 0.2 : 0.6;
azimuth =-60; %-60:30:60;
target = [{'female1'},{'female2'},{'female3'},{'male1'},{'male2'},{'male3'}];
k = 4; %1:6;
% load_dir = './1source/RT60_';
[s,fs] = audioread(['./1source/RT60_' , num2str(RT60),'/',target{k},'/mixture_1x2_SNR_0_',num2str(azimuth),'도/x_1x2.wav']);


x1 = s(:,1);
x2 = s(:,2);


% 초기화
g2 = [ zeros(1,M/2-1) 1 zeros(1,M/2)];
g1 = [0 zeros(1,M/2-1) zeros(1,M/2)];

% u 업데이트
u = [g1 g2]';
len = floor((length(x1)-M)/M);

for n = 1 : 1000000

    if (mod(n,len))~=0
        i = mod(n,len);
    else
        i = len;
    end
    
    x_temp(1,:) = x1((i*M) : -1 : (i-1)*M+1);
    x_temp(2,:) = x2((i*M) : -1 : (i-1)*M+1);
    x = [x_temp(2,:) -x_temp(1,:)]'; % 논문에서 x1*g2 = x2*g1 을 수행하기 위해, g는 그대로 두고 x의 좌우를 change

    e = u' * x;
    u_temp = u - mu*e*x;
    u = u_temp / norm(u_temp,2);
   
    if (mod(n,10000)==0)
    disp(['iteration = ' , num2str(n)])
    end
    
end

% for n = 1 : 5000000
% 
%     if (mod(n,len))~=0
%         i = mod(n,len);
%     else
%         i = len;
%     end
%     
%     x_temp(1,:) = x1((i*M) : -1 : (i-1)*M+1);
%     x_temp(2,:) = x2((i*M) : -1 : (i-1)*M+1);
%     x = [x_temp(1,:) x_temp(2,:)]'; % 논문에서 x1*g2 = x2*g1 을 수행하기 위해, g는 그대로 두고 x의 좌우를 change
% 
%     u = [g2 -g1]';
%     
%     e = u' * x;
%     u_temp = u - mu*e*x;
%     u = u_temp / norm(u_temp,2);
%     
%     g2 = u(1 : M)';
%     g1 = u(M+1 : 2*M)';
%     
% end


%  u = [g2 -g1]';
 
[a, b]=max(u(1:M));
[c, d]=max(u(M+1:2*M));

max_delay = micdist*fs/C;
est_sample_delay = -(d-b);
true_sample_delay = micdist*sin(azimuth*pi/180)*fs / C;
est_azimuth = (asin(est_sample_delay * C /fs / micdist))*180/pi;
true_azimuth = (asin(true_sample_delay * C /fs / micdist))*180/pi;


disp(['Max delay = ' ,num2str(max_delay)]);
disp(['True azimuth = ' ,num2str(true_azimuth), '˚']);
disp(['Est azimuth = ' ,num2str(est_azimuth), '˚']);
disp(['True Sample delay = ' ,num2str(true_sample_delay)]);
disp(['Est Sample delay = ' ,num2str(est_sample_delay)]);


figure(1)
subplot(211)
stem(u(1:M));
xlabel('k');ylabel('g_1(k)');
subplot(212)
stem(u(M+1:2*M));
xlabel('k');ylabel('g_2(k)');