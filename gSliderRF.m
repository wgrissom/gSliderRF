% Script to design gSlider pulses.
addpath rf_tools/ % JP's tools: gets dinf, b2a, cabc2rf, abr...
N = 128; % # time points in pulse
G = 5; % gSlider factor
Gpulse = 'ex';
tbG = 12; % overall tb product of encoding pulse
tbOther = 8; % tb product of non-encoding pulse
usecvx = false; % use Boyd's cvx toolbox for beta filter design
dt = 2.5e-3; % ms, final dwell time of pulses
T = 11; % ms, pulse duration of gSlider pulse; other pulse duration will be tbOther/tbG*T
slThick = 3.3; % mm, gSlider slice thickness
otherThickFactor = 1.15; % factor to increase slice thickness of non-gSlider pulse
gSlew = 150; % mT/m/ms, gradient slew rate for ramps
DFTphs = false; % do DFT phases
if strcmp(Gpulse,'ex')
    bsf = sqrt(1/2); % excitation pulse
    d1 = 0.001;d2 = 0.01; % passband and stopband ripples of the overall profile
    d1 = sqrt(d1/2); % Mxy passband ripple
    d2 = d2/sqrt(2); % Mxy stopband ripple
    d1O = 0.01;d2O = 0.01; % passband and stopband ripples of the se profile
    phi = pi; % slice phase - e.g., for ex, will be pi, for se, will be pi/2
elseif strcmp(Gpulse,'se')
    bsf = 1; % spin echo pulse
    d1 = 0.001;d2 = 0.01; % passband and stopband ripples of the overall profile
    d1 = d1/4;
    d2 = sqrt(d2);
    d1O = 0.01;d2O = 0.01; % passband and stopband ripples of the ex profile
    phi = pi/2; % slice phase
end

% print out some info about what we are doing
printf('--------gSlider RF Pulse Design---------');
printf('Designing %s gSlider encoding pulses.',Gpulse);
printf('Number of sub-slices: %d',G);
printf('Slab time-bandwidth product: %g',tbG);
printf('Duration: %d ms',T);
printf('Slice Thickness: %g mm',slThick);
printf('Output dwell time: %g ms',dt);
printf('Time-bandwidth product of other pulse: %g',tbOther);
printf('Slice thickness ratio of other pulse: %g',otherThickFactor);

rfEnc = zeros(N,ceil(G/2));
Mxy = zeros(N*8,ceil(G/2));
nomFlips = zeros(ceil(G/2),1);
for Gind = 1:ceil(G/2) % sub-slice to design for
    
    printf('Designing pulse for sub-slice %d of %d',Gind,G);
    
    % design the beta filter
    if ~DFTphs
        if usecvx
            printf('Designing beta filter using cvx');
            b = bsf*gSliderBeta_cvx(N,G,Gind,tbG,d1,d2,phi,8);
        else % use firls (less accurate but fast)
            printf('Designing beta filter using firls');
            b = bsf*gSliderBeta(N,G,Gind,tbG,d1,d2,phi);
        end
    else
        printf('Designing DFT beta filter using firls');
        phs = 2*pi/G*(ceil(-G/2):ceil(G/2)-1)*(Gind+ceil(-G/2)-1);
        if strcmp(Gpulse,'se'); phs = phs./2; end
        b = bsf*gSliderBetaDFT(N,phs,tbG,d1,d2);
    end
    
    % scale and solve for rf - note that b2rf alone doesn't work bc
    % b is not flipped correctly wrt a for non-symmetric profiles
    rfEnc(:,Gind) = cabc2rf(b2a(b),fliplr(b)); % iSLR for min-power pulse
    
    % simulate the pulse
    [ap,bp] = abr(rfEnc(:,Gind),-N/2:1/8:N/2-1/8);
    % calculate target flip angle of pulse in degrees
    nomFlips(Gind) = 2*asin(abs(bp(length(bp)/2+1)))*180/pi; 
    if strcmp(Gpulse,'ex')
        Mxy(:,Gind) = 2*conj(ap).*bp.*exp(1i*2*pi/(8*N)*N/2*(0:8*N-1)');
    elseif strcmp(Gpulse,'se')
        Mxy(:,Gind) = bp.^2;
    end
    
end

% plot all the profiles
zG = (-N/2:1/8:N/2-1/8)*slThick/tbG;
figure;
for ii = 1:ceil(G/2)
    if strcmp(Gpulse,'ex')
        titleText = sprintf('EX profile; gSlider factor %d; sub-slice %d',G,ii);
    else
        titleText = sprintf('SE pulse; gSlider factor %d; sub-slice %d',G,ii);
    end
    
    % plot the first encoding pulse
    subplot(ceil(G/2)*100 + 10 + ii),hold on
    plot(zG,abs(Mxy(:,ii)));
    plot(zG,real(Mxy(:,ii)));
    plot(zG,imag(Mxy(:,ii)));
    title(titleText);
    legend('|Mxy|','Mx','My');
    xlabel 'mm'
    axis([min(zG) max(zG) -1 1]);
end

% design and simulate the other pulse
if strcmp(Gpulse,'ex')
    Gother = 'se';
else
    Gother = 'ex';
end
rfOther = dzrf(N,tbOther,Gother,'ls',d1O,d2O);
[apO,bpO] = abr(rfOther,-N/2:1/8:N/2-1/8);
if strcmp(Gpulse,'ex')
    MxyO = bpO.^2;
elseif strcmp(Gpulse,'se')
    MxyO = 2*conj(apO).*bpO.*exp(1i*2*pi/(8*N)*N/2*(0:8*N-1)');
end

% plot both the pulse profiles
if strcmp(Gpulse,'ex')
    plind = 1;
    plindO = 2;
    titleText = sprintf('EX profile; gSlider factor = %d; sub-slice 1',G);
    titleTextO = 'SE profile';
else
    plind = 2;
    plindO = 1;
    titleText = sprintf('SE pulse; gSlider factor = %d; sub-slice 1',G);
    titleTextO = 'EX profile';
end
% plot the first encoding pulse
figure;subplot(310+plind),hold on
zG = (-N/2:1/8:N/2-1/8)*slThick/tbG;
plot(zG,abs(Mxy(:,1)));
plot(zG,real(Mxy(:,1)));
plot(zG,imag(Mxy(:,1)));
title(titleText);
legend('|Mxy|','Mx','My');
xlabel 'mm'
axis([min(zG) max(zG) -1 1]);

% plot the other pulse
subplot(310+plindO),hold on
zO = (-N/2:1/8:N/2-1/8)*slThick*otherThickFactor/tbOther;
plot(zO,abs(MxyO));
plot(zO,real(MxyO));
plot(zO,imag(MxyO));
title(titleTextO);
legend('|Mxy|','Mx','My');
xlabel 'mm'
axis([min(zG) max(zG) -1 1]);

% plot them together, to compare thicknesses
subplot(313),hold on
plot(zG,abs(Mxy(:,1)));
plot(zO,abs(MxyO));
legend('|Mxy|','|Mxy|, other');
xlabel 'mm'
axis([-2*slThick 2*slThick 0 1]);
title 'Both profile amplitudes'

% interpolate pulses to target dwell time
Nout = T/dt;
rfEncOut = zeros(Nout,G);
for ii = 1:ceil(G/2)
    rfEncOut(:,ii) = interp1((0:N-1)./N*T,rfEnc(:,ii),(0:Nout-1)*dt,'spline',0);
    rfEncOut(:,ii) = rfEncOut(:,ii)./sum(real(rfEncOut(:,ii)))*sum(real(rfEnc(:,ii)));
end
% time-reverse to get opposite pulses
rfEncOut(:,ceil(G/2)+1:end) = fliplr(flipud(rfEncOut(:,1:floor(G/2))));
Tother = T*tbOther/tbG/otherThickFactor;
NoutOther = round(Tother/dt);
rfOtherOut = interp1((0:N-1)./N*Tother,rfOther,(0:NoutOther-1)*dt,'spline',0);
rfOtherOut = rfOtherOut./sum(real(rfOtherOut))*sum(real(rfOther));

% convert to uT
rfEncOut = rfEncOut./(2*pi*42.58*dt*10^-3);
rfOtherOut = rfOtherOut./(2*pi*42.58*dt*10^-3);

gAmp = tbG/(T*10^-3)/(slThick*10^-3)/42580; % mT/m, gradient amplitude

% plot the final pulses
tRFEnc = (0:size(rfEncOut,1)-1)*dt;
tRFOther = (0:length(rfOtherOut)-1)*dt;
figure
subplot(210+plind),hold on
plot(tRFEnc,real(rfEncOut));
plot(tRFEnc,imag(rfEncOut));
plot(tRFEnc,abs(rfEncOut));
title(sprintf('%s gSlider Pulses',Gpulse));
xlabel 'ms'
ylabel 'uT'
legend('|RF|','Re\{RF\}','Im\{RF\}');

subplot(210+plindO),hold on
plot(tRFOther,real(rfOtherOut));
plot(tRFOther,imag(rfOtherOut));
plot(tRFOther,abs(rfOtherOut));
if strcmp(Gpulse,'ex')
    title(sprintf('SE pulse'));
else
    title(sprintf('EX pulse'));
end
xlabel 'ms'
ylabel 'uT'

% write out the pulses, gradient amplitude, target slice thickness,
% expected flip at middle of slice
% TODO: Should we phase shift the 90 or 180 by pi/2 for CPMG?
save(sprintf('gSliderRF_%dx_%s',G,Gpulse),'rfEncOut','rfOtherOut',...
    'gAmp','slThick','dt','nomFlips');
