f_samp = 330e3;

%Band Edge speifications
fs1 = 39.3e3;
fp1 = 43.3e3;
fp2 = 63.3e3;
fs2 = 67.3e3;

Wc1 = fp1*2*pi/f_samp;
Wc2  = fp2*2*pi/f_samp;

%Kaiser paramters
A = -20*log10(0.15);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

w_t = 2*pi*(fp1 - fs1)/f_samp;
N_min = ceil((A-8) / (2.285*w_t));           %empirical formula for N_min

%Window length for Kaiser Window
n=N_min + 18;

%Ideal bandpass impulse response of length "n"
bp_ideal = ideal_lp(pi*(fp2 + fs2)/f_samp,n) - ideal_lp(pi*(fp1 + fs1)/f_samp,n);

%Kaiser Window of length "n" with shape paramter beta calculated above
kaiser_win = (kaiser(n,beta))';

FIR_BandPass = bp_ideal .* kaiser_win;
fvtool(FIR_BandPass);         %frequency response

%magnitude response
[H,f] = freqz(FIR_BandPass,1,1024, f_samp);
plot(f,abs(H), "DisplayName", "Magnitude Response")
hold on 
title("Magnitude Response of Bandpass FIR Filter")
xlabel("Frequency (Hz)")
ylabel("Magnitude")
plot(f, ones(size(f))*0.15, '--', "DisplayName", "0.15");
plot(f, ones(size(f))*0.85, '--', "DisplayName", "0.85");
plot(f, ones(size(f))*1.15, '--', "DisplayName", "1.15");
plot([1 1]*fs1, [0 1.15], '--', "DisplayName", "f_{s1}");
plot([1 1]*fp1, [0 1.15], '--', "DisplayName", "f_{p1}");
plot([1 1]*fp2, [0 1.15], '--', "DisplayName", "f_{p2}");
plot([1 1]*fs2, [0 1.15], '--', "DisplayName", "f_{s2}");
grid
legend
hold off