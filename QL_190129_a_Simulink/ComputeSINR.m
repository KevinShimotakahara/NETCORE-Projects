
function [snr_wb_dB, snr_sb_dB] = ComputeSINR(RxWaveform, IntRxWaveform,...
    nRBGs, BW, TxInfo)
           
%%            
Fs          = TxInfo.SamplingRate;
RB          = 200e3;
fsubband    = [0:49] * RB;
N           = length(RxWaveform);
Res         = Fs/N;
f           = 0 : Res : Fs-Res;
RBbins      = floor(RB/Res);

%% Noise power
No_wb       = -174 + 10 * log10(BW) - 30;  
No_wb       = 10^(No_wb/10);
Pnoise_wb   = No_wb;

No_sb       = -174 + 10 * log10(180e3) - 30;
No_sb       = 10^(No_sb/10);
Pnoise_sb   = No_sb;

%% Signal power
Sig_fft         = fft(RxWaveform);

Psig_wb = bandpower(Sig_fft)/length(RxWaveform);
Psig_sb = zeros(size(fsubband));

for i = 1:length(fsubband)
    Psig_sb(i) = bandpower(Sig_fft((i-1)*RBbins+1 : i*RBbins))/length(RxWaveform);
end

%% Interference power
Int_fft     = fft(IntRxWaveform);

Pint_wb = bandpower(Int_fft)/length(IntRxWaveform);
Pint_sb = zeros(size(fsubband));

for i = 1:length(fsubband)
    Pint_sb(i) = bandpower(Int_fft((i-1)*RBbins+1 : i*RBbins))/length(IntRxWaveform);
end

%% snr calculations
snr_wb      = Psig_wb / (Pint_wb + Pnoise_wb);
snr_wb_dB   = 10 * log10(snr_wb);

snr_sb      = Psig_sb ./ (Pint_sb + Pnoise_sb);
snr_sb_dB   = 10 * log10(snr_sb);

%% Calculate SINR for RBGs

nRBs = 50/nRBGs;

for i = 1:nRBGs
    snr_temp(i) = mean(snr_sb_dB( (i-1)*nRBs+1 : (i*nRBs) ));
end

snr_sb_dB = snr_temp';

end




