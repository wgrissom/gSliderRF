% Test script to design gSlider pulses.
% Start by designing beta filters.
N = 128; % # time points in pulse
G = 5; % gSlider factor
Gind = 2; % sub-slice to design for
Gpulse = 'se';
tbG = 12; % overall tb product of encoding pulse
tbOther = 8; % tb product of non-encoding pulse
usecvx = false; % use Boyd's cvx toolbox for beta filter design
dt = 2.5e-3; % ms, final dwell time of pulses
T = 11; % ms, pulse duration of gSlider pulse; other pulse duration will be tbOther/tbG*T
slThick = 3.3; % mm, gSlider slice thickness
otherThickFactor = 1.15; % factor to increase slice thickness of non-gSlider pulse
gSlew = 150; % mT/m/ms, gradient slew rate for ramps
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

% print out some info about what we are doing
printf('--------gSlider RF Pulse Design---------');
printf('Designing an %s gSlider encoding pulse.',Gpulse);
printf('Number of bands: %d',G);
printf('Band index: %d',Gind);
printf('Time-bandwidth product: %d',tbG);
printf('Duration: %d ms',T);
printf('Slice Thickness: %d mm',slThick);
printf('Output dwell time: %d ms',dt);
printf('Time-bandwidth product of other pulse: %d',tbOther);
printf('Slice thickness ratio of other pulse: %d',otherThickFactor);

% design the beta filter
if usecvx
  printf('Designing beta filter using cvx');
  b = bsf*gSliderBeta_cvx(N,G,Gind,tbG,d1,d2,phi,8);
else % use firls (less accurate but fast)
  printf('Designing beta filter using firls');
  b = bsf*gSliderBeta(N,G,Gind,tbG,d1,d2,phi);
end

% plot the beta profile
%figure;hold on
%B = fftshift(fft(b(:).',8*N).*exp(1i*2*pi/(8*N)*N/2*(0:8*N-1)));
%if strcmp(Gpulse,'se')
%  B = B.^2;
%end
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
rfEnc = cabc2rf(a,fliplr(b)); % iSLR for min-power pulse

% simulate the pulse
[ap,bp] = abr(rfEnc,-N/2:1/8:N/2-1/8);
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

% design the other pulse
if strcmp(Gpulse,'ex')
  Gother = 'se';
else
  Gother = 'ex';
end
rfOther = dzrf(N,tbOther,Gother,'ls',d1,d2);

% interpolate pulses to target dwell time
Nout = T/dt;
rfEncOut = interp1((0:N-1)./N*T,rfEnc,(0:Nout-1)*dt,'spline',0);
rfEncOut = rfEncOut./sum(real(rfEncOut))*sum(real(rfEnc));
Tother = T*tbOther/tbG/otherThickFactor;
NoutOther = round(Tother/dt);
rfOtherOut = interp1((0:N-1)./N*Tother,rfOther,(0:NoutOther-1)*dt,'spline',0);
rfOtherOut = rfOtherOut./sum(real(rfOtherOut))*sum(real(rfOther));

% convert to uT - TODO: Convert to Volts
rfEncOut = rfEncOut./(2*pi*42.58*dt*10^-3);
rfOtherOut = rfOtherOut./(2*pi*42.58*dt*10^-3);

% define gradients at 10 us, to be same duration as pulses
dtG = 10e-3;
gAmp = tbG/(T*10^-3)/(slThick*10^-3)/42580; % mT/m
gEnc = gAmp + zeros(ceil(length(rfEncOut)*dt/dtG),1);
gOther = gAmp + zeros(ceil(length(rfOtherOut)*dt/dtG),1);
nRamp = ceil(gAmp/(gSlew*dtG));
gRamp = (1:nRamp-1)./nRamp*gAmp;
gEnc = [gRamp(:);gEnc(:);flipud(gRamp(:))];
gOther = [gRamp(:);gOther(:);flipud(gRamp(:))];

% zero-pad pulses to same duration as gradients
nRFpad = length(gRamp)*dtG/dt;
rfEncOut = [zeros(nRFpad,1);rfEncOut(:);zeros(nRFpad,1)];
rfOtherOut = [zeros(nRFpad,1);rfOtherOut(:);zeros(nRFpad,1)];

% plot the pulses
tRFEnc = (0:length(rfEncOut)-1)*dt;
tRFOther = (0:length(rfOtherOut)-1)*dt;
tGEnc = (0:length(gEnc)-1)*dtG;
tGOther = (0:length(gOther)-1)*dtG;
figure
subplot(211),hold on
plot(tRFEnc,real(rfEncOut));
plot(tRFEnc,imag(rfEncOut));
plot(tRFEnc,abs(rfEncOut));
plot(tGEnc,gEnc);
title(sprintf('%s pulse',Gpulse));
xlabel 'ms'
ylabel 'uT or mT/m'

subplot(212),hold on
plot(tRFOther,real(rfOtherOut));
plot(tRFOther,imag(rfOtherOut));
plot(tRFOther,abs(rfOtherOut));
plot(tGOther,gOther);
title(sprintf('Other pulse'));
xlabel 'ms'
ylabel 'uT or mT/m'

% write out the pulses - TODO: Should we phase shift the 90 or 180 by pi/2 for CPMG? 
