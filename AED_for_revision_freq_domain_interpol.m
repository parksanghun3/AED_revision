%%%%%%%%% AED ( CDR MASK 논문 Revision 목적) %%%%%%%%%%%
%%%%%% 주파수 도메인으로 변경%%%%%%%%%%%
clear all;
close all; clc;

%Parameters
M = 2; L = 256; % length of impulse response vectors
C= 343;
micdist = 1;
fs = 16000;
max_delay = micdist*fs/C;
IntpRatio = 64;
SNR = -1;
RT60 = 0.4; %0.2 : 0.2 : 0.6;
azimuth = -60 : 30 : 60 ;
target = [{'female1'},{'female2'},{'female3'},{'male1'},{'male2'},{'male3'}];
Nu = 0; % utterence 개수
Nu_out = 0; % est_sample_delay가 max_delay보다 큰 경우는 Nu에서 제외함.

mu=3.3E-1;
mu2=7.1E-3;
lambda = (1-1/3/L)^L; 
m_frame = 100; scale_delta = 0.2; 

for target_index = 1 : 6
    for azimuth_index = 1 : 5

        err=zeros(10000000,1);
        n1=zeros(1000,L);
        n2=zeros(1000,L);
        flag=0;
        dError1=zeros(10000000,1);
        dError2=zeros(10000000,1);
        dRatio=zeros(10000000,1);
        
        % Factors
        index = 0:(L-1);
        FL = exp(-1i*2*pi/L*(index'*index));
        index = 0:2*L-1;
        F2L = exp(-1i*2*pi/(2*L)*(index'*index));
        
        W01 = FL*[zeros(L) eye(L)]*inv(F2L);
        W10 = F2L*[eye(L); zeros(L)]*inv(FL);
        
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
        
        x1=s(:,1);x2=s(:,2);
        
      
        % Initialization
        h1_init=[zeros(1,L/2-1) 1/sqrt(M) zeros(1,L/2)]';
        h2_init=[zeros(1,L/2-1) 1/sqrt(M) zeros(1,L/2)]';
        
        % Calculation
        h1 = fft(h1_init);
        h2 = fft(h2_init);
        h110_ = fft([h1_init; zeros(L,1)]);
        h210_ = fft([h2_init; zeros(L,1)]);
        
        m_Dx = zeros(2*L,1);
        for m=1:m_frame
            m_t_Dx = fft(x1((m-1)*L+1:(m+1)*L));
            m_Dx = m_Dx + abs(m_t_Dx).^2;
            m_t_Dx = fft(x2((m-1)*L+1:(m+1)*L));
            m_Dx = m_Dx + abs(m_t_Dx).^2;
        end
        
        delta = min(m_Dx)/M/m_frame*scale_delta;
        P1 = zeros(2*L,1);
        P2 = zeros(2*L,1);
        
        len=floor((length(x1)-L)/L);
        
        err(1)=0;
        
        for n=1:100000
            if (mod(n,len)~=0)
                m= mod(n,len);
            else
                m=len;
            end
            l = mod(n,1000);
            Dx1 = fft(x1((m-1)*L+1:(m+1)*L));
            Dx2 = fft(x2((m-1)*L+1:(m+1)*L));
            e12_ = W01*(Dx1.*(W10*h2) - Dx2.*(W10*h1));
            e21_ = -e12_;
            err(n) =  norm(e12_,2);
            
            % For Debugging
            dE1=W01*(Dx1.*(W10*h2));
            dE2=W01*(Dx2.*(W10*h1));
            dError1(n)=norm(dE1,2);
            dError2(n)=norm(dE2,2);
            dRatio(n)=err(n)/dError2(n);
            % End Debugging
            
            P1 = lambda*P1 + (1-lambda)*abs(Dx2).^2;
            P2 = lambda*P2 + (1-lambda)*abs(Dx1).^2;
            e1201 = fft([zeros(L,1); real(ifft(e12_))]);
            e2101 = -e1201;
            tmp_h1 = h110_- mu*conj(Dx2).*e2101./(P1+delta*ones(2*L,1));
            tmp_h2 = h210_ - mu*conj(Dx1).*e1201./(P2+delta*ones(2*L,1));
            th1 = real(ifft(tmp_h1));
            th2 = real(ifft(tmp_h2));
            norm1=norm((th1(1:L)),2);
            norm2=norm((th2(1:L)),2);
            h_th1 = th1(1:L) - norm1*mu2*sign(th1(1:L));
            h_th2 = th2(1:L) - norm2*mu2*sign(th2(1:L));
            
            % Normalization
            norm_h = norm([h_th1; h_th2]);
            n_th1 = h_th1/norm_h;
            n_th2 = h_th2/norm_h;
            if (l==0)
                flag=flag+1;
                n1(flag,:) = n_th1;
                n2(flag,:) = n_th2;
            end
            h1 = fft(n_th1);
            h2 = fft(n_th2);
            h110_ = fft([n_th1; zeros(L,1)]);
            h210_ = fft([n_th2; zeros(L,1)]);
            
            if n == 100000 %% h 업데이트 끝나고 마지막
                h1_padding = [h1(1:L/2) ; zeros(L * (IntpRatio-1)/2, 1) ; h1(L/2 +1) ; zeros(L * (IntpRatio-1)/2, 1) ; h1(end-L/2 +2 : end)];
                h2_padding = [h2(1:L/2) ; zeros(L * (IntpRatio-1)/2, 1) ; h2(L/2 +1) ; zeros(L * (IntpRatio-1)/2, 1) ; h2(end-L/2 +2 : end)];
                
                th1 = real(ifft(h1_padding));
                th2 = real(ifft(h2_padding));
               
                norm1=norm(th1,2);
                norm2=norm(th2,2);
                h_th1_final = th1 - norm1*mu2*sign(th1);
                h_th2_final = th2 - norm2*mu2*sign(th2);
                
                % Normalization
                norm_h = norm([h_th1_final; h_th2_final]);
                n_th1_final = h_th1_final/norm_h;
                n_th2_final = h_th2_final/norm_h;                
            end
        end
        
        [a, b]=max(n_th1_final);
        [c, d]=max(n_th2_final);

        est_sample_delay = ((b-d)/IntpRatio);
        
       if abs(est_sample_delay) >= max_delay
            if est_sample_delay > 0
                est_sample_delay = max_delay;
            else
                est_sample_delay = -max_delay;
            end
            
        end
        
        true_sample_delay = micdist*sin(azimuth(azimuth_index)*pi/180)*fs / C ; 
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

disp('=============================================================');
disp(['sum(diff_sample_delay) = ' ,num2str(sum(diff_sample_delay))]);
MAE = sum(diff_sample_delay) / Nu;
disp(['Nu = ' , num2str(Nu)]);
disp(['MAE = ' , num2str(MAE)]);

% figure(1);
% subplot(211);stem(n_th1);title('Adaptive Filter h1');
% subplot(212);stem(n_th2);title('Adaptive Filter h2');
% subplot(211);stem(n_th1_final);title('Adaptive Filter h1');
% subplot(212);stem(n_th2_final);title('Adaptive Filter h2');
