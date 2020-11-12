%% Set input parameters
m = 3;
delta = 0.15;
f_s = 330;

%% Get band edges
q = floor(0.1*m);
r = m - 10*q;
f_pl = 25 + 1.7*q + 6.1*r;
f_ph = f_pl + 20;
f_sl = f_pl - 4;
f_sh = f_ph + 4;
%% Normalize and convert to low-pass analog
f = [f_sl f_pl f_ph f_sh];
w = 2*pi*f/f_s;
sigma = tan(w/2);
sigma_0 = sqrt(sigma(2)*sigma(3));
B = sigma(3) - sigma(2);
sigma_l = (sigma.^2 - sigma_0^2)./(B*sigma);
sigma_lp = 1;
sigma_ls = min(abs([sigma_l(1), sigma_l(4)]));

%% Analog Lowpass Transfer function
D1 = 1/(1-delta)^2 - 1;
epsilon = sqrt(D1);
D2 = 1/delta^2 - 1;
N = ceil(acosh(sqrt(D2/D1))/acosh(sigma_ls/sigma_lp));
B_k = asinh(1/epsilon)/N;
angles = pi*((2*N+1)/(2*N):1/N:(4*N-1)/(2*N));
poles = sinh(B_k)*sin(angles) + cosh(B_k)*cos(angles)*1i;

% This part is specific to N=4
p1 = [1 -2*real(poles(1)) abs(poles(1)^2)];
p2 = [1 -2*real(poles(2)) abs(poles(2)^2)];
denom = conv(p1, p2);
num = denom(5)/sqrt(1+epsilon^2);

syms s z;
H_analog_lpf(s) = poly2sym(num, s)/poly2sym(denom, s);

%% Bandpass and Discrete Transfer Function
H_analog_bpf(s) = H_analog_lpf((s^2 + sigma_0^2)/(s*B));
H_discrete_bpf = H_analog_bpf((z-1)/(z+1));
[num_discrete, denom_discrete] = numden(H_discrete_bpf);
num_discrete = sym2poly(num_discrete);
denom_discrete = sym2poly(denom_discrete);
k = denom_discrete(1);
num_discrete = num_discrete/k;
denom_discrete = denom_discrete/k;
[H,f] = freqz(num_discrete, denom_discrete, 1000, 330e3);

%% Required Plots
fvtool(num_discrete, denom_discrete, 'Analysis', 'freq')
figure
hold on
plot(f,abs(H), "DisplayName", "Magnitude Response")
title("Magnitude Response of Bandpass IIR Filter")
xlabel("Frequency (Hz)")
ylabel("Magnitude")
plot(f, ones(size(f))*0.15, '--', "DisplayName", "0.15");
plot(f, ones(size(f))*0.85, '--', "DisplayName", "0.85");
plot(f, ones(size(f))*1.15, '--', "DisplayName", "1.15");
plot([1 1]*f_sl*1e3, [0 1.15], '--', "DisplayName", "f_{s1}");
plot([1 1]*f_pl*1e3, [0 1.15], '--', "DisplayName", "f_{p1}");
plot([1 1]*f_ph*1e3, [0 1.15], '--', "DisplayName", "f_{p2}");
plot([1 1]*f_sh*1e3, [0 1.15], '--', "DisplayName", "f_{s2}");
% plot(f,abs(angle(H)), "DisplayName", "Phase Response")
grid
legend
hold off