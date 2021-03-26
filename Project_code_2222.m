%% Project code

%Processing earthquake data to show characteristics and amplification of
%the earthquake in two sites in Nepal
%Alex vanderhoeff
%250858436

%%
%reading data
close all
fileID1=fopen('TVU001504250611.dat');
fileID2=fopen('THM001504250611.dat');
formatSpec='%f%f%f%f';
TVUdata=textscan(fileID1,formatSpec,'HeaderLines',34);
THMdata=textscan(fileID2,formatSpec,'HeaderLines',34);
fclose(fileID1);
fclose(fileID2);
t1=TVUdata{1,1};
t4=THMdata{1,1};
time_THM=THMdata{1,1};
time_TVU=TVUdata{1,1};
%%
%detrending and removing mean from data
%TVU
TVU2 = detrend(TVUdata{1,2}); %remove linear trend 
TVU22 = TVU2 - mean(TVUdata{1,2}); %remove mean offset
TVU3 = detrend(TVUdata{1,3}); %remove linear trend 
TVU33 = TVU3 - mean(TVUdata{1,3}); %remove mean offset
TVU4 = detrend(TVUdata{1,4}); %remove linear trend 
TVU44 = TVU4 - mean(TVUdata{1,4}); %remove mean offset

% %THM
THM2 = detrend(THMdata{1,2}); %remove linear trend 
THM22 = THM2 - mean(THMdata{1,2}); %remove mean offset
THM3 = detrend(THMdata{1,3}); %remove linear trend 
THM33 = THM3 - mean(THMdata{1,3}); %remove mean offset
THM4 = detrend(THMdata{1,4}); %remove linear trend 
THM44 = THM4 - mean(THMdata{1,4}); %remove mean offset

%%
%create the filter
fs = 100 ;                 % sampling frequency
dt = 1/fs;
order = 2;                 % order of the filter
low_cut = 0.05;              % low-cut cutoff frequency
high_cut = 20;                 %high-cut cutoff freq
[b1,a1]=butter(2,[low_cut high_cut]/(fs/2),'bandpass'); %bandpass filter for both stations
%%
%apply the butter filter
%THM
THM22_filt_accel=filtfilt(b1,a1,THM22);
THM33_filt_accel=filtfilt(b1,a1,THM33);
THM44_filt_accel=filtfilt(b1,a1,THM44);

%TVU
TVU22_filt_accel=filtfilt(b1,a1,TVU22);
TVU33_filt_accel=filtfilt(b1,a1,TVU33);
TVU44_filt_accel=filtfilt(b1,a1,TVU44);

%%
%compute velocities
v2_THM = cumsum(THM22_filt_accel)*dt;% (cm/s)
v3_THM = cumsum(THM33_filt_accel)*dt;% (cm/s)
v4_THM = cumsum(THM44_filt_accel)*dt;% (cm/s)

v2_TVU = cumsum(TVU22_filt_accel)*dt;% (cm/s)
v3_TVU = cumsum(TVU33_filt_accel)*dt;% (cm/s)
v4_TVU = cumsum(TVU44_filt_accel)*dt;% (cm/s)

%%
%compute displacement
d2_THM = cumsum(v2_THM)*dt;%(cm)
d3_THM = cumsum(v3_THM)*dt;%(cm)
d4_THM = cumsum(v4_THM)*dt;%(cm)

d2_TVU = cumsum(v2_TVU)*dt;%(cm)
d3_TVU = cumsum(v3_TVU)*dt;%(cm)
d4_TVU = cumsum(v4_TVU)*dt;%(cm)



%% Plots
%Acceleration
%THM22
figure(1)
subplot(3,1,1)
plot(time_THM,THM22_filt_accel,'k');
ylabel('Acceleration (cm/s^2)');
legend('N-S component');
axis([0 100 -250 250]);
title('THM station acceleration time series of Gorkha earthquake');

subplot(3,1,2)
plot(time_THM,THM33_filt_accel,'r');
ylabel('Acceleration (cm/s^2)');
legend('E-W component');
axis([0 100 -250 250]);

subplot(3,1,3)
plot(time_THM,THM44_filt_accel,'b');
ylabel('Acceleration (cm/s^2)');
legend('Verical component');
axis([0 100 -250 250]);
xlabel('Time (s)');

%TVU
figure(2)
subplot(3,1,1)
plot(time_TVU,TVU22_filt_accel,'k');
ylabel('Acceleration (cm/s^2)');
legend('N-S component');
axis([0 100 -250 250]);
title('TVU station acceleration time series of Gorkha earthquake');

subplot(3,1,2)
plot(time_TVU,TVU33_filt_accel,'r');
ylabel('Acceleration (cm/s^2)');
legend('E-W component');
axis([0 100 -250 250]);

subplot(3,1,3)
plot(time_TVU,TVU44_filt_accel,'b');
ylabel('Acceleration (cm/s^2)');
legend('Verical component');
axis([0 100 -250 250]);
xlabel('Time (s)');

%% Velocity
%THM22
figure(3)
subplot(3,1,1)
plot(time_THM,v2_THM,'k');
ylabel('Velocity (cm/s)');
legend('N-S component');
axis([0 100 -250 250]);
title('THM station velocity time series of Gorkha earthquake');

subplot(3,1,2)
plot(time_THM,v3_THM,'r');
ylabel('Velocity (cm/s)');
legend('E-W component');
axis([0 100 -250 250]);

subplot(3,1,3)
plot(time_THM,v4_THM,'b');
ylabel('Velocity (cm/s)');
legend('Verical component');
axis([0 100 -250 250]);
xlabel('Time (s)');

%TVU
figure(4)
subplot(3,1,1)
plot(time_TVU,v2_TVU,'k');
ylabel('Velocity (cm/s)');
legend('N-S component');
axis([0 100 -250 250]);
title('TVU station velocity time series of Gorkha earthquake');

subplot(3,1,2)
plot(time_TVU,v3_TVU,'r');
ylabel('Velocity (cm/s)');
legend('E-W component');
axis([0 100 -250 250]);

subplot(3,1,3)
plot(time_TVU,v4_TVU,'b');
ylabel('Velocity (cm/s)');
legend('Verical component');
axis([0 100 -250 250]);
xlabel('Time (s)');

%% Displacement
%THM22
figure(5)
subplot(3,1,1)
plot(time_THM,d2_THM,'k');
ylabel('Displacement (cm)');
legend('N-S component');
axis([0 90 -100 100]);
title('THM station displacement time series of Gorkha earthquake');

subplot(3,1,2)
plot(time_THM,d3_THM,'r');
ylabel('Displacement (cm)');
legend('E-W component');
axis([0 90 -100 100]);

subplot(3,1,3)
plot(time_THM,d4_THM,'b');
ylabel('Displacement (cm)');
legend('Verical component');
axis([0 90 -100 100]);
xlabel('Time (s)');

%TVU
figure(6)
subplot(3,1,1)
plot(time_TVU,d2_TVU,'k');
title('TVU station displacement time series of Gorkha earthquake');
ylabel('Displacement (cm)');
legend('N-S component');
axis([0 90 -100 100]);

subplot(3,1,2)
plot(time_TVU,d3_TVU,'r');
ylabel('Displacement (cm)');
legend('E-W component');
axis([0 90 -100 100]);

subplot(3,1,3)
plot(time_TVU,d4_TVU,'b');
ylabel('Displacement (cm)');
legend('Verical component');
axis([0 90 -100 100]);
xlabel('Time (s)');

%% Compute PGA, PGV, max displacement
THM22vel_max=max(abs(v2_THM));
THM33vel_max=max(abs(v3_THM));
THM44vel_max=max(abs(v4_THM));
THM22_filt_accel_max=max(abs(THM22_filt_accel)) ;
THM33_filt_accel_max=max(abs(THM33_filt_accel)) ; 
THM44_filt_accel_max=max(abs(THM44_filt_accel));
THM22_d_max=max(abs(d2_THM));
THM33_d_max=max(abs(d3_THM)); 
THM44_d_max=max(abs(d4_THM));

TVU22vel_max=max(abs(v2_TVU));
TVU33vel_max=max(abs(v3_TVU));
TVU44vel_max=max(abs(v4_TVU));
TVU22accel_max=max(abs(TVU22_filt_accel));
TVU33accel_max=max(abs(TVU33_filt_accel));
TVU44accel_max=max(abs(TVU44_filt_accel));
TVU22_d_max=max(abs(d2_TVU));
TVU33_d_max=max(abs(d3_TVU));
TVU44_d_max=max(abs(d4_TVU));

%% CPSD
figure(7)
subplot(3,1,1);
nfft_1=2^nextpow2(length(TVU22_filt_accel));
[cpsd_Spectra_22,f_spectra_22] = cpsd(TVU22_filt_accel,THM22_filt_accel,[ ],0,nfft_1,fs);
m1=abs(cpsd_Spectra_22);
semilogx(f_spectra_22,m1,'b');
title('Cross power spectral density between TVU and THM');
legend('N-S componenets');
ylabel('Power(W)');
axis([10^(-2) 15 0 2000]);

subplot(3,1,2);
nfft_2=2^nextpow2(length(TVU33_filt_accel));
[cpsd_Spectra_33,f_spectra_33] = cpsd(TVU33_filt_accel,THM33_filt_accel,[ ],0,nfft_2,fs);
m2=abs(cpsd_Spectra_33);
semilogx(f_spectra_33,m2,'r');
legend('E-W componenets');
ylabel('Power(W)');
axis([0 15 0 2500]);

subplot(3,1,3);
nfft_3=2^nextpow2(length(TVU44_filt_accel));
[cpsd_Spectra_44,f_spectra_44] = cpsd(TVU44_filt_accel,THM44_filt_accel,[ ],0,nfft_3,fs);
m3=abs(cpsd_Spectra_44);
semilogx(f_spectra_44,m3,'m');
legend('Vertical componenets');
xlabel('Frequency(Hz)');
ylabel('Power(W)');
axis([10^(-2) 15 0 2000]);
%% PSD 
%TVU
figure(8)
subplot(3,1,1)
nfft1=2^nextpow2(length(TVU22_filt_accel)); 
[p_TVU22,f_TVU22]=periodogram(TVU22_filt_accel,[],nfft1,200);
PSD_TVU22=abs(p_TVU22);
semilogx(f_TVU22,PSD_TVU22,'b');
axis([10^(-1) 50 0 6000]);
title('Power spectral density of N-S acceleration at station TVU');
xlabel('frequency (Hz)');
ylabel('Power(w)');

subplot(3,1,2)
nfft2=2^nextpow2(length(TVU33_filt_accel)); 
[p_TVU33,f_TVU33]=periodogram(TVU33_filt_accel,[],nfft2,200);
PSD_TVU33=abs(p_TVU33);
semilogx(f_TVU33,PSD_TVU33,'r');
axis([10^(-1) 50 0 6000]);
title('Power spectral density of E-W acceleration at station TVU');
xlabel('frequency (Hz)');
ylabel('Power(w)');

subplot(3,1,3)
nfft3=2^nextpow2(length(TVU44_filt_accel)); 
[p_TVU44,f_TVU44]=periodogram(TVU44_filt_accel,[],nfft3,200);
PSD_TVU44=abs(p_TVU44);
semilogx(f_TVU44,PSD_TVU44,'m');
axis([10^(-1) 50 0 6000]);
title('Power spectral density of vertical acceleration at station TVU');
xlabel('frequency (Hz)');
ylabel('Power(w)');

%THM
figure(9)
subplot(3,1,1);
nfft4=2^nextpow2(length(THM22_filt_accel)); 
[p_THM22,f_THM22]=periodogram(THM22_filt_accel,[],nfft4,200);
PSD_THM22=abs(p_THM22);
semilogx(f_THM22,PSD_THM22,'b');
axis([10^(-1) 50 0 6000]);
title('Power spectral density of N-S acceleration at station THM');
xlabel('frequency (Hz)');
ylabel('Power(w)');

subplot(3,1,2);
nfft5=2^nextpow2(length(THM33_filt_accel)); 
[p_THM33,f_THM33]=periodogram(THM33_filt_accel,[],nfft5,200);
PSD_THM33=abs(p_THM33);
semilogx(f_THM33,PSD_THM33,'r');
axis([10^(-1) 50 0 6000]);
title('Power spectral density of E-W acceleration at station THM');
xlabel('frequency (Hz)');
ylabel('Power(w)');

subplot(3,1,3);
nfft6=2^nextpow2(length(THM44_filt_accel)); 
[p_THM44,f_THM44]=periodogram(THM44_filt_accel,[],nfft6,200);
PSD_THM44=abs(p_THM44);
semilogx(f_THM44,PSD_THM44,'m');
axis([10^(-1) 50 0 6000]);
title('Power spectral density of vertical acceleration at station THM');
xlabel('frequency (Hz)');
ylabel('Power(w)');

%% Spectrogram
%THM
figure(10)
subplot(3,1,1);
spectrogram(THM22_filt_accel,[],50,nfft4,fs);
title('Spectrograms for station THM');
subplot(3,1,2);
spectrogram(THM33_filt_accel,[],50,nfft5,fs);
subplot(3,1,3);
spectrogram(THM44_filt_accel,[],50,nfft6,fs);


%TVU
figure(11)
subplot(3,1,1);
spectrogram(TVU22_filt_accel,[],50,nfft1,fs);
title('Spectrograms for station TVU');
subplot(3,1,2);
spectrogram(TVU33_filt_accel,[],50,nfft2,fs);
subplot(3,1,3);
spectrogram(TVU44_filt_accel,[],50,nfft3,fs);


%% Fourier transforms
%TVU
N_TVU = 2^nextpow2(length(TVU22_filt_accel));               
FT_TVU22 = abs(fft(TVU22_filt_accel)*dt);               %FT
FT_TVU33 = abs(fft(TVU33_filt_accel)*dt); 
FT_TVU44 = abs(fft(TVU44_filt_accel)*dt);
f_TVU22 = (0:length(FT_TVU22)-1)/length(FT_TVU22)*fs;     %Frequency
f_TVU33 = (0:length(FT_TVU33)-1)/length(FT_TVU33)*fs;  
f_TVU44 = (0:length(FT_TVU44)-1)/length(FT_TVU44)*fs;  

%THM
N_THM = 2^nextpow2(length(THM22_filt_accel));
FT_THM22 = abs(fft(THM22_filt_accel)*dt);
FT_THM33 = abs(fft(THM33_filt_accel)*dt); 
FT_THM44 = abs(fft(THM44_filt_accel)*dt);
f_THM22 = (0:length(FT_THM22)-1)/length(FT_THM22)*fs;
f_THM33 = (0:length(FT_THM33)-1)/length(FT_THM33)*fs;
f_THM44 = (0:length(FT_THM44)-1)/length(FT_THM44)*fs;



%% Smooth FT
%THM
coef=(ones(5,1))/5; 
smoothft_THM22 = filtfilt(coef,1,abs(FT_THM22));
smoothft_THM33 = filtfilt(coef,1,abs(FT_THM33));
smoothft_THM44 = filtfilt(coef,1,abs(FT_THM44));

%TVU
smoothft_TVU22 = filtfilt(coef,1,abs(FT_TVU22));
smoothft_TVU33 = filtfilt(coef,1,abs(FT_TVU33));
smoothft_TVU44 = filtfilt(coef,1,abs(FT_TVU44));

%% Plot fourier transforms
%THM smoothed and unsmoothed plot
figure(12)
%N-S
subplot(3,1,1);
semilogx(f_THM22,FT_THM22,'r');
axis([10^(-2) 25 0 800]);
hold on
semilogx(f_THM22,smoothft_THM22,'b');
title('Smoothed and unsmoothed Fourier transforms (FT) for station THM');
legend('Unsmoothed FT for N-S', 'Smoothed FT for N-S');
ylabel('Acceleration (cm/s^2)');
hold off

%E-W
subplot(3,1,2);
semilogx(f_THM33,FT_THM33,'k');
axis([10^(-2) 25 0 800]);
hold on
semilogx(f_THM33,smoothft_THM33,'b');
ylabel('Acceleration (cm/s^2)');
legend('Unsmoothed FT for E-W', 'Smoothed FT for E-W');
hold off

%Vertical
subplot(3,1,3);
semilogx(f_THM44,FT_THM44,'g');
axis([10^(-2) 25 0 800]);
hold on
semilogx(f_THM44,smoothft_THM44,'b');
xlabel('Frequency (Hz)');
ylabel('Acceleration (cm/s^2)');
legend('Unsmoothed FT for Vertical', 'Smoothed FT for Vertical');
hold off


figure(13)
subplot(3,1,1);
semilogx(f_TVU22,FT_TVU22,'r');
axis([10^(-2) 25 0 800]);
hold on
semilogx(f_TVU22,smoothft_TVU22,'b');
ylabel('Acceleration (cm/s^2)');
title('Smoothed and unsmoothed Fourier transforms (FT) for station TVU');
legend('Unsmoothed FT for N-S', 'Smoothed FT for N-S');

subplot(3,1,2);
semilogx(f_TVU33,FT_TVU33,'m');
axis([10^(-2) 25 0 800]);
hold on
semilogx(f_TVU33,smoothft_TVU33,'b');
ylabel('Acceleration (cm/s^2)');
legend('Unsmoothed FT for E-W', 'Smoothed FT for E-W');
hold off

subplot(3,1,3);
semilogx(f_TVU44,FT_TVU44,'k');
axis([10^(-2) 25 0 800]);
hold on
semilogx(f_TVU44,smoothft_TVU44,'b');
ylabel('Acceleration (cm/s^2)');
xlabel('Frequency (Hz)');
legend('Unsmoothed FT for Vertical', 'Smoothed FT for Vertical');
hold off

%% Spectral ratio
%HVSR
%compute hvsr TVU
geomean_TVU = sqrt(smoothft_TVU22.*smoothft_TVU33);
hv_av_TVU = geomean_TVU./smoothft_TVU44;

%compute hvsr THM
geomean_THM = sqrt(smoothft_THM22.*smoothft_THM33);
hv_av_THM = geomean_THM./smoothft_THM44;

figure(14)
semilogx(f_TVU22,hv_av_TVU,'g');
axis([10^(-2) 25 0 20]);
hold on
semilogx(f_THM22,hv_av_THM,'r');
xlabel('Frequency (Hz)');
ylabel('Amplification');
title('HV spectral ratios');
legend('Station TVU','Station THM');
