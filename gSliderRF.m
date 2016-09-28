% Test script to design gSlider pulses.
% Start by designing beta filters.
N = 128; % # time points in pulse
G = 5; % gSlider factor
Gind = 2; % sub-slice to design for
Gpulse = 'ex';
tbG = 12; % overall tb product of encoding pulse
tbOther = 8; % tb product of non-encoding pulse
usecvx = false;
if strcmp(Gpulse,'ex')
  bsf = sqrt(1/2); % excitation pulse
  d1 = 0.01;d2 = 0.01; % passband and stopband ripples of the overall profile
  d1 = sqrt(d1/2); % Mxy passband ripple
  d2 = d2/sqrt(2); % Mxy stopband ripple
  phi = pi; % slice phase - e.g., for ex, will be pi, for se, will be pi/2
elseif strcmp(Gpulse,'se')
  bsf = 1; % spin echo pulse
  d1 = 0.01;d2 = 0.01; % passband and stopband ripples of the overall profile
  d1 = d1/4;
  d2 = sqrt(d2);
  phi = pi/2; % slice phase
end

% design the beta filter
if usecvx
  b = bsf*gSliderBeta_cvx(N,G,Gind,tbG,d1,d2,phi,8);
else % use firls (less accurate)
  b = bsf*gSliderBeta(N,G,Gind,tbG,d1,d2,phi);
end

% plot the beta profile
%figure;hold on
B = fftshift(fft(b(:).',8*N).*exp(1i*2*pi/(8*N)*N/2*(0:8*N-1)));
if strcmp(Gpulse,'se')
  B = B.^2;
end
%plot(real(B));
%plot(imag(B));
%plot(abs(B));

% Did a conventional design with same characteristics for comparison,
% while exploring islr symmetry issue. 
% [rf12,b12] = dzrf(N,tbG,Gpulse,'ls',d1,d2);
% b12 = b12.*exp(1i*2*pi/N*4*(0:N-1));
% rf12 = b2rf(b12);
% B12 = fftshift(fft(b12(:).',8*N).*exp(1i*2*pi/(8*N)*N/2*(0:8*N-1)));
% if strcmp(Gpulse,'se')
%   B12 = B12.^2;
% end
% [ap12,bp12] = abr(rf12,-N/2:1/8:N/2-1/8);
% if strcmp(Gpulse,'ex')
%   Mxy12 = 2*conj(ap12).*bp12.*exp(1i*2*pi/(8*N)*N/2*(0:8*N-1)');
% elseif strcmp(Gpulse,'se')
%   Mxy12 = bp12.^2;
% end

% scale and solve for rf - note that b2rf alone doesn't work bc 
% b is not flipped correctly wrt a for non-symmetric profiles
a = b2a(b);
rf = cabc2rf(a,fliplr(b)); % iSLR for min-power pulse

% simulate the pulse
[ap,bp] = abr(rf,-N/2:1/8:N/2-1/8);
if strcmp(Gpulse,'ex')
  Mxy = 2*conj(ap).*bp.*exp(1i*2*pi/(8*N)*N/2*(0:8*N-1)');
elseif strcmp(Gpulse,'se')
  Mxy = bp.^2;
end

figure;hold on
plot(abs(Mxy));
plot(real(Mxy));
plot(imag(Mxy));
title(sprintf('gSlider factor = %d; Slice index = %d',G,Gind));
legend('|Mxy|','Mx','My');
