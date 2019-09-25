function b = gSliderBetaHadamard_pn(N,G,ind,tb,d1,d2)

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
    % left stopband
    f = [0 shift-(1+ftw)*(tb/2)]; 
    m = [0 0];
    w = d1/d2;
    % first sub-band
    ii = 1;
    Gcent = shift+(ii-G/2-1/2)*tb/G; % first band center
    f = [f Gcent-(tb/G/2-ftw*(tb/2))]; % first band left edge
    m = [m encode(ii)]; % first band left edge amp
    if encode(ii) ~= encode(ii+1)
        % add a weight for first band, and the right edge and its amp
        f = [f Gcent+(tb/G/2-ftw*(tb/2))];
        m = [m encode(ii)];
        w = [w 1];
    end
    % middle sub-bands
    for ii = 2:length(encode)-1
        Gcent = shift+(ii-G/2-1/2)*tb/G; % center of this band
        if encode(ii) ~= encode(ii-1)
            % add a left edge and amplitude for this band
            f = [f Gcent-(tb/G/2-ftw*(tb/2))];
            m = [m encode(ii)];
        end
        if encode(ii) ~= encode(ii+1)
            % add a weight for this band, and the right edge and its amp
            f = [f Gcent+(tb/G/2-ftw*(tb/2))];
            m = [m encode(ii)];
            w = [w 1];
        end
    end
    % last sub-band
    ii = length(encode);
    Gcent = shift+(ii-G/2-1/2)*tb/G; % last band center
    if encode(ii) ~= encode(ii-1)
        f = [f Gcent-(tb/G/2-ftw*(tb/2))]; % last band left edge
        m = [m encode(ii)]; % last band left edge amp
    end
    % add a weight for last band, and the right edge and its amp
    f = [f Gcent+(tb/G/2-ftw*(tb/2))];
    m = [m encode(ii)];
    w = [w 1];
    
    % add the right stopband
    f = [f shift+(1+ftw)*(tb/2) (N/2)]./(N/2); % right stopband
    %m = [0 0 kron(encode,[1 1]) 0 0];
    m = [m 0 0];
    mp = double(m > 0);
    mn = double(m < 0);
    w = [w d1/d2];
    bp = firls(N-1,f,mp,w); % the filter
    bn = firls(N-1,f,mn,w);
    b = hilbert(bp - bn).*exp(-1i*2*pi/N*shift*(0:N-1))/2*exp(-1i*pi/N*shift);
end
