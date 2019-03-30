clc; close all; clear all;
N = 50; Cn = 1/sqrt(N); E = 1; T = 0.001; Ts = 5;
wm = 140; % beta*v; 50 in hertz; 2*pi*f/c *v=139.6
shi = -pi + 2*pi.*rand(1,N);
phi = -pi + 2*pi.*rand(1,N);
theta = -pi + 2*pi.*rand(1,N);
JakeCommat=[];
for m = 1 : T : Ts
    JakeRe = 0; JakeIm = 0; 
    for k = 1 : N
        alpha(k) = (2*pi*k-pi+theta(k))/( 4*N + 2 );
        Xc = E * Cn * (cos(wm*m*cos(alpha(k)) + phi(k)));
        Xs = E * Cn * (sin(wm*m*cos(alpha(k)) + phi(k)));
        JakeRe = JakeRe + Xc;
        JakeIm = JakeIm + Xs;
    end    
    JakeCom = JakeRe + 1i*JakeIm;
    JakeCommat = [JakeCommat JakeCom];
end

% AutoCorrelation 
rea = real(JakeCommat);
imo = imag(JakeCommat);
envelope = sqrt(rea.^2 + imo.^2);
[a1,a2] = xcorr(JakeCommat,'Coeff');
zayz = xcorr(JakeCommat,'Coeff');

% Power Spectral Density
Fres = 1000; N = length(zayz); psdd = fft(zayz);
xdft = psdd(1:N/2+1);
psdx = (1/(Fres*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fres/length(psdd):Fres/2;

% Plotting of Autocorrelation & Power Spectral Density
len = length(a2);
[b1,b2] = xcorr(rea, imo, 'Coeff');
aa = -(len-1)/2 : 1 : (len-1)/2;
bb = (len-2001)./2;
cc = bb+1 : 1 : bb+2001;
dd = -1000 : 1 : 1000;
tdd = dd*T; z1 = wm*tdd;  sigma0 = 1;
T_bessel = sigma0.^2.*besselj(0,z1);
figure;
subplot(2,2,1);
plot(tdd,real(a1(cc)),'-',tdd,real(T_bessel),'.');
title('AutoCorrelation of Jakes Model'); xlabel('x-axis not scaled')
legend('simulated','theoritical');
subplot(2,2,2);
plot(freq,psdx); grid on;
title('PSD of AutoCorrelation of Jakes Model')
xlabel('Frequency (Hz)')

%
co1 = 1 : 1000;
subplot(2,2,3);
semilogy(co1*T,envelope(1:1000)); title('Rayleigh Coefficient'); xlabel('t(second)');
length_r = length(envelope);
pdf_env = zeros(1,501);
count=0;
temp = round(100.*envelope);
for k = 1 : length_r
    if temp(k) <= 500
        count = count + 1;
        pdf_env(1,temp(k)+1) = pdf_env(1,temp(k)+1)+1;
    end
end
pdf_env = pdf_env./count./0.01;
sgma2 = 0.5;
x = [0:0.01:5];
pdf_theory = (x./sgma2).*exp(-1.*x.^2./(2.*sgma2));

subplot(2,2,4);
plot(x,pdf_env,'-',x,pdf_theory,'*');
legend('Simulated','Theory'); title('PDF of jakes envelope');
