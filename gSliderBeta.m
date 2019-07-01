function [b,b1,b2] = gSliderBeta(N,G,Gind,tb,d1,d2,phi)

% Script to design gSlider beta filters.
ftw = dinf(d1,d2)/tb; % fractional transition width of the slab profile

if ~verLessThan('matlab','9.6')
    warning('Your MATLAB version is >= 9.6, forcing real-valued gSlider profile');
    realM = true;
else
    realM = false;
end

%% Design a gSlider beta polynomial - can get negative
% band positions by flipping the waveform
if rem(G,2) && Gind == ceil(G/2) % centered sub-slice
    if G == 1 % no sub-slices, as a sanity check
        f = [0 (1-ftw)*(tb/2) (1+ftw)*(tb/2) (N/2)]/(N/2);
        m = [1 1 0 0];
        w = [1 d1/d2];
        b = firls(N-1,f,m,w); % the filter
    else
        f = [0 (1/G-ftw)*(tb/2) (1/G+ftw)*(tb/2) (1-ftw)*(tb/2) (1+ftw)*(tb/2) (N/2)]/(N/2);
        m = [exp(1i*phi) exp(1i*phi) 1 1 0 0];
        if realM
            m = real(m);
        end
        w = [1 1 d1/d2];
        b = firls(N-1,f,m,w); % the filter
    end
else
    % design two shifted filters that we can add to kill off the left band,
    % then demodulate the result back to DC
    shift = N/4;
    Gcent = shift+(Gind-G/2-1/2)*tb/G;
    if Gind > 1 && Gind < G
        % separate transition bands for slab+slice
        f = [0 shift-(1+ftw)*(tb/2) shift-(1-ftw)*(tb/2) ...
            Gcent-(tb/G/2+ftw*(tb/2)) Gcent-(tb/G/2-ftw*(tb/2)) ...
            Gcent+(tb/G/2-ftw*(tb/2)) Gcent+(tb/G/2+ftw*(tb/2)) ...
            shift+(1-ftw)*(tb/2) shift+(1+ftw)*(tb/2) (N/2)]/(N/2);
        m = [0 0 1 1 exp(1i*phi) exp(1i*phi) 1 1 0 0];
        w = [d1/d2 1 1 1 d1/d2];
    elseif Gind == 1
        % the slab and slice share a left transition band
        f = [0 shift-(1+ftw)*(tb/2) shift-(1-ftw)*(tb/2) ...
            Gcent+(tb/G/2-ftw*(tb/2)) Gcent+(tb/G/2+ftw*(tb/2)) ...
            shift+(1-ftw)*(tb/2) shift+(1+ftw)*(tb/2) (N/2)]/(N/2);
        m = [0 0 exp(1i*phi) exp(1i*phi) 1 1 0 0];
        w = [d1/d2 1 1 d1/d2];
    elseif Gind == G
        % the slab and slice share a right transition band
        f = [0 shift-(1+ftw)*(tb/2) shift-(1-ftw)*(tb/2) ...
            Gcent-(tb/G/2+ftw*(tb/2)) Gcent-(tb/G/2-ftw*(tb/2)) ...
            shift+(1-ftw)*(tb/2) shift+(1+ftw)*(tb/2) (N/2)]/(N/2);
        m = [0 0 1 1 exp(1i*phi) exp(1i*phi) 0 0];
        w = [d1/d2 1 1 d1/d2];
    end
    if realM
        m = real(m);
    end
    b = firls(N-1,f,m,w); % the filter
    % hilbert transform to suppress negative passband, and then demod to DC
    b = hilbert(b).*exp(-1i*2*pi/N*shift*(0:N-1))/2*exp(-1i*pi/N*shift);
end
