function b = gSliderBetaHadamard(N,G,ind,tb,d1,d2)

H = hadamard(G);
encode = H(ind,:);

% Script to design gSlider beta filters.
ftw = dinf(d1,d2)/tb; % fractional transition width of the slab profile

%% Design a gSlider beta polynomial - can get negative
% band positions by flipping the waveform
if ind == 1
    % No phase, simple uniform pulse
    f = [0 (1-ftw)*(tb/2) (1+ftw)*(tb/2) (N/2)]/(N/2);
    m = [1 1 0 0];
    w = [1 d1/d2];
    b = firls(N-1,f,m,w); % the filter
else
    % Each sub-slice has its own phase;
    % design two shifted filters that we can add to kill off the left band,
    % then demodulate the result back to DC
    shift = N/4;
    f = [0 shift-(1+ftw)*(tb/2)]; % left stopband
    for ii = 1:length(encode)
        Gcent = shift+(ii-G/2-1/2)*tb/G;
        f = [f Gcent-(tb/G/2-ftw*(tb/2)) Gcent+(tb/G/2-ftw*(tb/2))];
    end
    f = [f shift+(1+ftw)*(tb/2) (N/2)]./(N/2); % right stopband
    m = [0 0 kron(encode,[1 1]) 0 0];
    w = [d1/d2 ones(1,G) d1/d2];
    b = firls(N-1,f,m,w); % the filter
    b = hilbert(b).*exp(-1i*2*pi/N*shift*(0:N-1))/2*exp(-1i*pi/N*shift);
end
