
function [snr_wb_dB] = sinr_fn(RxWaveform, IntRxWaveform, bandwidth, txInfo)
        
RMS = @(x) sqrt(mean(abs(x).^2));

%%            
fs = txInfo.SamplingRate;
RB = 200e3;
Nfft = txInfo.Nfft;
resolution = fs/Nfft;

%% Noise power
wbN0 = -174 + 10 * log10(bandwidth) - 30;  
wbN0 = 10^(wbN0/10);
wbPnoise = wbN0;               % wideband noise power

%% wideband signal power
wbPsignal = RMS(RxWaveform)^2;

%% Interference power
wbPint = RMS(IntRxWaveform)^2; 

%% snr calculations
snr_wb      = wbPsignal / (wbPint + wbPnoise);
snr_wb_dB   = 10 * log10(snr_wb);

end




