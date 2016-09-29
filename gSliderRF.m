% Test script to design gSlider pulses.
% Start by designing beta filters.
N = 128; % # time points in pulse
G = 5; % gSlider factor
Gind = 2; % sub-slice to design for
Gpulse = 'ex';
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

% design and simulate the other pulse
if strcmp(Gpulse,'ex')
  Gother = 'se';
else
  Gother = 'ex';
end
rfOther = dzrf(N,tbOther,Gother,'ls',d1,d2);
[apO,bpO] = abr(rfOther,-N/2:1/8:N/2-1/8);
if strcmp(Gpulse,'ex')
  MxyO = bpO.^2;
elseif strcmp(Gpulse,'se')
  MxyO = 2*conj(apO).*bpO.*exp(1i*2*pi/(8*N)*N/2*(0:8*N-1)');
end

% plot both the pulse profiles
if strcmp(Gpulse,'ex')
    plind = 311;
    plindO = 312;
    titleText = sprintf('EX profile; gSlider factor = %d; Slice index = %d',G,Gind);
    titleTextO = 'SE pulse';
else
    plind = 312;
    plindO = 311;
    titleText = sprintf('SE pulse; gSlider factor = %d; Slice index = %d',G,Gind);
    titleTextO = 'EX pulse';
end
% plot the encoding pulse
figure;subplot(plind),hold on
zG = (-N/2:1/8:N/2-1/8)*slThick/tbG;
plot(zG,abs(Mxy));
plot(zG,real(Mxy));
plot(zG,imag(Mxy));
title(titleText);
legend('|Mxy|','Mx','My');
xlabel 'mm'

% plot the other pulse
subplot(plindO),hold on
zO = (-N/2:1/8:N/2-1/8)*slThick*otherThickFactor/tbOther;
plot(zO,abs(MxyO));
plot(zO,real(MxyO));
plot(zO,imag(MxyO));
title(titleTextO);
legend('|Mxy|','Mx','My');
xlabel 'mm'

% plot them together, to compare thicknesses
subplot(313),hold on
plot(zG,abs(Mxy));
plot(zO,abs(MxyO));
legend('|Mxy|','|Mxy|, other');
xlabel 'mm'

% interpolate pulses to target dwell time
Nout = T/dt;
rfEncOut = interp1((0:N-1)./N*T,rfEnc,(0:Nout-1)*dt,'spline',0);
rfEncOut = rfEncOut./sum(real(rfEncOut))*sum(real(rfEnc));
Tother = T*tbOther/tbG/otherThickFactor;
NoutOther = round(Tother/dt);
rfOtherOut = interp1((0:N-1)./N*Tother,rfOther,(0:NoutOther-1)*dt,'spline',0);
rfOtherOut = rfOtherOut./sum(real(rfOtherOut))*sum(real(rfOther));

% convert to uT
rfEncOut = rfEncOut./(2*pi*42.58*dt*10^-3);
rfOtherOut = rfOtherOut./(2*pi*42.58*dt*10^-3);

gAmp = tbG/(T*10^-3)/(slThick*10^-3)/42580; % mT/m, gradient amplitude

% plot the final pulses
tRFEnc = (0:length(rfEncOut)-1)*dt;
tRFOther = (0:length(rfOtherOut)-1)*dt;
figure
subplot(211),hold on
plot(tRFEnc,real(rfEncOut));
plot(tRFEnc,imag(rfEncOut));
plot(tRFEnc,abs(rfEncOut));
title(sprintf('%s pulse',Gpulse));
xlabel 'ms'
ylabel 'uT'
legend('|RF|','Re{RF}','Im{RF}');

subplot(212),hold on
plot(tRFOther,real(rfOtherOut));
plot(tRFOther,imag(rfOtherOut));
plot(tRFOther,abs(rfOtherOut));
title(sprintf('Other pulse'));
xlabel 'ms'
ylabel 'uT'

% write out the pulses, gradient amplitude, target slice thickness, 
% expected flip at middle of slice
% TODO: Should we phase shift the 90 or 180 by pi/2 for CPMG? 
