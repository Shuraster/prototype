%This task was done jointly by students of the group TKS-22 Eremin Pavel and Morozov Egor


clc;clear all;close all;

%% Inp parameter
Fm = 1e4; % Основная частота сигнала сообщения
Harm = [ 1 0.5 2 1 ]; % The frequency's component of message signal  
Ampl = [ 1 2 3 1 ]; % Amplitute the frequency component
sampling_rate = 1/(20*max(Fm*Harm)); % Sampling rate
range = 2/min(Fm*Harm); % Time range
t = 0:sampling_rate:range; % Syncronization

%% Signal message
message = zeros(size(t));
for k=1:length(Harm)
    message = message + Ampl(k)*sin(2*pi*Harm(k)*Fm*t);
end
minAmpl = min(message); % Min ampl in model
%%Sampling
n_sample = 5;
Fs = n_sample*max(Harm*Fm); % Frequency sampling
m_samp = zeros(size(message));
m_samp(1:1/(Fs*sampling_rate):length(t)) = message(1:1/(sampling_rate*Fs):length(t));% Out sampling signal
figure(1);
plot(t,message);
grid on 
hold on
stem(t,m_samp);
xlabel('Time');
ylabel('Signal message, sampling & Quantization signal');
title('Signal message, sampling & Quantization signal');
%legend('Signal message','Sampling signal','Quantization signal');

%% Quantization
kv_levels = 4; % Level quantizations                  
quantile = (max(m_samp) - min(m_samp))/(kv_levels); % Interval quantizations
resp_levels = min(m_samp):quantile:max(m_samp); % Representation lvl
Q = zeros(size(m_samp));                                                          
for k = 1:length(resp_levels)
    values = (( m_samp>(resp_levels(k)-quantile/2) & m_samp<(resp_levels(k)+quantile/2)));   
    Q(values) = round(resp_levels(k)); % Quantizations signal of message
end
clear values;
stem(t,Q,'r*');
grid on
%legend('Signal message, Sampling signal & Quantization signal');

%% Coding
if min(Q) >= 0
    minAmpl = 0;                                                        
end
Q1 = Q-round(minAmpl); % Offset negative values ??on the positive side for conversion to binary sampling
Bits = de2bi(Q1(1:1/(Fs*sampling_rate):length(Q)),4,'left-msb')';          
Bits = Bits(:)'; % Generating a binary sequence of a sampled quantized signal
Bits = 2 * Bits - ones(1,length(Bits));
figure(5);
stem(Bits,'*r');
hold on;
%legend('Bit sequence in transmitter','Bit sequence in receiver');

% Spreading 
spreadSeq = generateCAcode(13);%[1 1 1 1 1 -1 -1 1 -1 1 -1 -1 -1];%
spreadingCode = repmat(spreadSeq,1, length(Bits));
maskUpsampling = length(spreadSeq)*ones(1, length(Bits));
upsampleBits = repelem(Bits, maskUpsampling);

spreadingBits = spreadingCode .* upsampleBits;

figure(101);
stem(spreadingBits(1:5*length(spreadSeq)))

%% BPSK
Fc = 1e7; % Carrier frequency
Nsamp = 10; % System samples per carrier signal cycle 
Ncyc = 2; % The number of cycles of the carrier signal for one bit interval                                                        
Tbit = 0:1/(Nsamp*Fc):Ncyc/Fc; % Bit interval
t_per = 0:1/(Nsamp*Fc):(Ncyc*length(spreadingBits))/Fc+(length(spreadingBits)-1)/(Nsamp*Fc); % Total transmission time  
modulation_sig = zeros(size(t_per));                                               
l = 1;
for k = 1:length(Tbit):length(modulation_sig)
    if(spreadingBits(l) == 1)
        modulation_sig(k:k+length(Tbit)-1) = cos(2*pi*Fc*Tbit); % Carrier phase modulation to represent binary symbol 1
    else
        modulation_sig(k:k+length(Tbit)-1) = -cos(2*pi*Fc*Tbit); % Carrier phase modulation to represent binary symbol 0
    end
    l = l+1;
end

%% AWGN - channel
Per_sig = awgn(modulation_sig,10); % Transmission of a modulated carrier signal over the AWGN channel
figure(2);
plot(t_per,modulation_sig,'.-b',t_per,Per_sig,'r');
axis([0 3*Ncyc/Fc -2 2]);
xlabel('Time');
ylabel('Tx - signal transmitter & Tx - signal with noise transmitter');
title('Signal transmission and signal with noise transmission');
%legend('Signal transmitt ','Signal with noise transmitt');

%% Receiver
% Filter
F_freq = -(Nsamp*Fc)/2:(Nsamp*Fc)/length(t_per):(Nsamp*Fc)/2-(Nsamp*Fc)/length(t_per); % The frequency range used to visualize the FFT signal 
f_spread = fft(spreadingBits);
f_pr = fft(Per_sig); % Fast Fourier transform                              
figure(3);
plot(F_freq, fftshift(f_pr), F_freq, fftshift(fft(modulation_sig)), 'g', F_freq(1:21:end),  fftshift(f_spread), 'b');
grid on;
xlabel('Frequency');
ylabel('Amplitude signal');
%legend('Modulated signal','Received signal');
title('Modulated signal & Received signal');
F_rece = zeros(size(f_pr));               
FIR = (F_freq < -3*Fc | F_freq>3*Fc);
F_rece(FIR) = f_pr(FIR); % Noise filtering in frequency domain to remove noise
F_rece(~FIR) = 0.5*f_pr(~FIR);                                                                                                                      
t_rece = ifft(F_rece); % Delete noise from signal

figure(4);
plot(t_per,t_rece);
grid on;
xlabel('Time');
ylabel('Received signal');
title('Received signal after noise filtering');
delete f_freq f_tran f_rece;
clear f_freq f_tran  f_rece;

%% Demodulation
binSeq = zeros(size(spreadingBits));
k = 1:length(Tbit):length(t_per);

for  iter = 1:length(binSeq)% Extract binary data from media using correlation method                                           
        a = corr(cos(2*pi*Fc*Tbit),t_rece(k(iter):k(iter)+length(Tbit)-1));                                                              
        b = real(a);
        if b > 0.5
            binSeq(iter) = 1;
        else
            binSeq(iter) = -1;
        end
        if (iter == length(binSeq))
            continue
        endif
end
a = zeros(size(spreadSeq));
seq2 = [spreadSeq spreadSeq];
for delay = 1:length(spreadSeq)
    filtSeqPart = seq2(delay:length(spreadSeq)-1+delay);
    filtSeq = repmat(filtSeqPart,1, length(Bits));
    a(delay) = corr(filtSeq, binSeq);
    if(max(abs(a))>0.03)
        continue
    endif
endfor

demod1 = filtSeq.*binSeq;
data = zeros(size(Bits));
demod_data = zeros(size(Bits));
for i = 1:length(Bits) %downsample
    
    data(i) = sum(demod1(length(spreadSeq)*(i-1)+1:length(spreadSeq)*i));
    if(data(i)>0.5)
      demod_data(i) = 1;
    elseif(data(i)<-0.5)
      demod_data(i) = -1;
    else  % ?????? ???? ?? ??????
      demod_data(i) = 0;
    endif
endfor

figure(5);
stem(demod_data);
grid on;
xlabel('Bit position');
ylabel('The bit sequence in the receiver and transmitter');
title('The bit sequence in the receiver and transmitter');
legend('The bit sequence in the receiver','The bit sequence in the transmitter');
demod_data = 1/2*(demod_data + once(size(demod_data)));

%% Decoding
demod_data = reshape(demod_data,4,length(Bits)/4)';
mQ_rece = zeros(size(Q));
mQ_rece(1:1/(Fs*sampling_rate):length(Q)) = bi2de(demod_data,'left-msb')' + min(Q); % Extract sampled quantized data from a decoded binary sequence

%% Signal conversion
F_freq = -1/(2*sampling_rate):1/(sampling_rate*length(t)):1/(2*sampling_rate)-1/(sampling_rate*length(t));                                  
F_rece = fft(mQ_rece); % FFT of the extracted quantized extracted signal with discretization  
F_out = zeros(size(F_rece));                                            
figure(6);
plot(F_freq,fftshift(F_rece),F_freq,fftshift((fft(message))));
grid on;
xlabel('Frequence');
ylabel('Signal conversion and signal conversion');
title('Input signal and received signal');
legend('Received signal','Input signal');
F_out((F_freq < -17000 | F_freq > 17000))=F_rece((F_freq < -17000 | F_freq > 17000)); % Frequency domain filtering to recover a signal from a sample of quantized data
out = ifft(F_out); % Reconstructed output signal
figure(7);
plot(t,4*out,t,message,'r');
grid on;
xlabel('time');
ylabel('Message signal and output signal');
title('Input signal and received signal');
legend('Received signa','Input signal');
gain = 4; % Power gain
out = out*gain; % Signal after power gain
figure(8);
plot(t,out);
grid on;
xlabel('Time');
ylabel('Received signal');
title('Received signal');

## %% This program simulates BER of BPSK in AWGN channel%
## num_bit=1e5;                          %Signal length 
## max_run=20;                              %Maximum number of iterations for a single SNR
## Eb=1;                                    %Bit energy
## SNRdB=0:1:10;                             %Signal to Noise Ratio (in dB)
## SNR=10.^(SNRdB/10);                      
## hand=waitbar(0,'Please Wait....');
## for count=1:length(SNR)                  %Beginning of loop for different SNR
##     avgError=0;
##     No=Eb/SNR(count);                    %Calculate noise power from SNR
##     for run_time=1:max_run               %Beginning of loop for different runs
##         waitbar((((count-1)*max_run)+run_time-1)/(length(SNRdB)*max_run));
##         Error=0;
##         data=randint(1,num_bit);         %Generate binary data source
##         s=2*data-1;                      %Baseband BPSK modulation
##         N=sqrt(No/2)*randn(1,num_bit);   %Generate AWGN
##         Y=s+N;                           %Received Signal
##         for k=1:num_bit                  %Decision device taking hard decision and deciding error
##             if ((Y(k)>0 && data(k)==0)||(Y(k)<0 && data(k)==1))
##                 Error=Error+1;
##             end
##         end
##         Error=Error/num_bit;             %Calculate error/bit
##         avgError=avgError+Error;         %Calculate error/bit for different runs        
##     end                                  %Termination of loop for different runs
##     BER_sim(count)=avgError/max_run;     %Calculate BER for a particular SNR                                  
## end                                      %Termination of loop for different SNR 
## BER_th=(1/2)*erfc(sqrt(SNR));            %Calculate analytical BER
## close(hand);
## semilogy(SNRdB,BER_th,'r');              %Plot BER
## hold on
## semilogy(SNRdB,BER_sim,'b*');
## legend('Theoretical','Simulation',3);
## axis([min(SNRdB) max(SNRdB) 10^(-5) 1]);
## hold off
## 
 %% Constellation
##%create a random digital message
##M=2; %alphabet size
##x=data;
##%use bpsk modulation to produce y
##y=modulate(modem.qammod(M),x);
##%%transmit signal through an AWGN Channel
##ynoisy=awgn(y,100,'measured'); %without noise
##ynoisy1=awgn(y,10,'measured'); %with noise
##%Create scattet plot from noisy data
##scatterplot(ynoisy),grid;title('Constellation without noise');
##
##scatterplot(ynoisy1),grid; title('Constellation with noise');

%% Spectrum plotting 
fftOriginalPadding=2^nextpow2(length(Bits)); 
fftOriginal=fft(Bits, fftOriginalPadding); 
fftOriginal=fftOriginal/length(Bits); 
fftOriginalAbsolute=abs(fftOriginal); 
fftOriginalAxis=Fs*(0:fftOriginalPadding-1)/fftOriginalPadding; 

fftSpreadedPadding=2^nextpow2(length(spreadingBits)); 
fftSpreaded=fft(spreadingBits,fftSpreadedPadding); 
fftSpreaded=fftSpreaded/length(spreadingBits); 
fftSpreadedAbsolute=abs(fftSpreaded); 
fftSpreadedAxis=Fs*(0:fftSpreadedPadding-1)/fftSpreadedPadding; 

figure; hold on; grid on; 
plot(fftSpreadedAxis,fftSpreadedAbsolute,'r'); 
plot(fftOriginalAxis,fftOriginalAbsolute,'b'); 
title('Spectrum of the original and the spreaded signals'); 
xlabel('Frequency, Hz'); 
ylabel('|FFTsignal(F)|'); 
legend('Spreaded inf. sequence','Information sequence');