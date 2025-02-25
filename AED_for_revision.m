%%%%%%%%% AED ( CDR MASK 논문 Revision 목적) %%%%%%%%%%%
%%
clear all; close all ; clc;

C= 343;
micdist = 1;
M = 512; % length of impulse response vectors
mu = 0.001;
% s = audioread('m60_68_FKFB0_con.wav');
% s = audioread('m45_68_MDAC0_con.wav');
fs = 16000;
IntpRatio = 16;
SNR = 20;
RT60 = 0.6; %0.2 : 0.2 : 0.6;
azimuth = -60 : 30 : 60 ;
target = [{'female1'},{'female2'},{'female3'},{'male1'},{'male2'},{'male3'}];
Nu = 0; % utterence 개수
Nu_out = 0; % est_sample_delay가 max_delay보다 큰 경우는 Nu에서 제외함.

for target_index = 1 : 6
    for azimuth_index = 1 : 5
        
        %%% 마이크 간격 0.12m
%         [s,fs_s] = audioread(['./Simulated Target/1source/RT60_' , num2str(RT60),'/',target{target_index},'/mixture_1x2_SNR_0_',num2str(azimuth(azimuth_index)),'도/x_1x2.wav']);
%         [n,fs_n] = audioread('./Simulated Noise/x_18x2.wav');
        
        %%% 마이크 간격 1m
        [s,fs_s] = audioread(['./1m_Target/RT60_' , num2str(RT60),'/',target{target_index},'/mixture_1x2_SNR_0_',num2str(azimuth(azimuth_index)),'도/x_1x2.wav']);
        [n,fs_n] = audioread(['./1m_Noise/RT60_' , num2str(RT60),'/diffuse_noise_4/', 'x_18x2.wav']);
        
        s_resample = resample(s,fs,fs_s);
        n_resample = resample(n,fs,fs_n);
        n_resample = [n_resample; n_resample];
        
        s_resample_vad = vad(s_resample);
        s = s_resample_vad;
        
        n_resample_resize = n_resample(1:length(s_resample_vad),:);
        n = n_resample_resize;
        
        
        s_norm = norm(s(:,1));
        n_norm = norm(n(:,1));
        energy_norm = s_norm/n_norm;
 
        if SNR ~= -1
            for l = 1:size(s,2)
                n(:,l) = n(:,l)*energy_norm/sqrt(10^(SNR/10));
                s(:,l) = s(:,l) + n(:,l);
            end
        end
        
%         s_energy = sqrt(sum(s_resample_vad(:,1).^2));
%         n_energy = sqrt(sum(n_resample_resize(:,1).^2));

        s1_interpol = interp(s(:,1),IntpRatio);
        s2_interpol = interp(s(:,2),IntpRatio);
        s_interpol = [s1_interpol, s2_interpol];
        
        if IntpRatio == 1
            x1 = s(:,1);
            x2 = s(:,2);
            
        elseif IntpRatio >= 2
            x1 = s_interpol(:,1);
            x2 = s_interpol(:,2);
            
        end
        
        
        % 초기화
        g2 = [ zeros(1,M*IntpRatio/2-1) 1 zeros(1,M*IntpRatio/2)];
        g1 = [0 zeros(1,M*IntpRatio/2-1) zeros(1,M*IntpRatio/2)];
        
        % u 업데이트
        u = [g1 g2]';
        len = floor((length(x1)-M*IntpRatio)/(M*IntpRatio));
        
        for n = 1 : 1000000
            
            if (mod(n,len))~=0
                i = mod(n,len);
            else
                i = len;
            end
            
            x_temp(1,:) = x1((i*M*IntpRatio) : -1 : (i-1)*M*IntpRatio+1);
            x_temp(2,:) = x2((i*M*IntpRatio) : -1 : (i-1)*M*IntpRatio+1);
            x = [x_temp(2,:) -x_temp(1,:)]'; % 논문에서 x1*g2 = x2*g1 을 수행하기 위해, g는 그대로 두고 x의 좌우를 change
            
            e = u' * x;
            u_temp = u - mu*e*x;
            u = u_temp / norm(u_temp,2);
            
        end
        [a, b]=max(u(1:M*IntpRatio));
        [c, d]=max(u(M*IntpRatio+1:2*M*IntpRatio));
        
        max_delay = micdist*fs/C;
        est_sample_delay = -(d-b)/IntpRatio;
        
        if abs(est_sample_delay) >= max_delay
            if est_sample_delay > 0 
             est_sample_delay = max_delay;
            else 
             est_sample_delay = -max_delay;
            end
            
       end
        
        true_sample_delay = micdist*sin(azimuth(azimuth_index)*pi/180)*fs / C;
        est_azimuth = (asin(est_sample_delay * C /fs / micdist))*180/pi;
        true_azimuth = (asin(true_sample_delay * C /fs / micdist))*180/pi;
        Nu = Nu + 1;
        diff_sample_delay(Nu,:) = abs(abs(est_sample_delay) - abs(true_sample_delay));
        
        
        disp(' ')
        disp(['현재 음원 = RT60(' , num2str(RT60), ')_', target{target_index} ,'_azimuth(', num2str(azimuth(azimuth_index)), ')_' , 'SNR(',num2str(SNR),')'])
%         disp(['현재 음원 = azimuth(', num2str(azimuth(azimuth_index)), ')_', target{target_index}, '_RT60(' , num2str(RT60),')_SNR(',num2str(SNR),')'])
        disp(['True azimuth = ' ,num2str(true_azimuth), '˚']);
        disp(['Est azimuth = ' ,num2str(est_azimuth), '˚']);
        disp(['True Sample delay = ' ,num2str(true_sample_delay)]);
        disp(['Est Sample delay = ' ,num2str(est_sample_delay)]);
            
    end
end

% Nu = Nu - Nu_out;
MAE = sum(diff_sample_delay) / Nu;

figure(1)
subplot(211)
stem(u(1:M*IntpRatio));
xlabel('k');ylabel('g_1(k)');
subplot(212)
stem(u(M*IntpRatio+1:2*M*IntpRatio));
xlabel('k');ylabel('g_2(k)');