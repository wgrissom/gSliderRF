function [b,bNotch,b2] = gSliderBeta(N,G,Gind,tb,d1,d2,phi)

% Script to design gSlider beta filters.
ftw = dinf(d1,d2)/tb; % fractional transition width of the slab profile

if rem(G,2) && Gind == ceil(G/2) % centered sub-slice
    if G == 1 % no sub-slices, as a sanity check
        f = [0 (1-ftw)*(tb/2) (1+ftw)*(tb/2) (N/2)]/(N/2);
        m = [1 1 0 0];
        w = [1 d1/d2];
        b = firls(N-1,f,m,w); % the filter
    else

        % Design 2 filters, to allow arbitrary phases on the subslice
        % The first is a wider notch filter with '0's where it the
        % subslice appears, and the second is the subslice.
        % Multiply the subslice by its phase and add the filters.
        
        f = [0 (1/G-ftw)*(tb/2) (1/G+ftw)*(tb/2) (1-ftw)*(tb/2) ...
            (1+ftw)*(tb/2) (N/2)]/(N/2);
        mNotch = [0 0 1 1 0 0];
        mSub = [1 1 0 0 0 0];
        w = [1 1 d1/d2];
        
        bNotch = firls(N-1,f,mNotch,w); % the notched filter
        bSub = firls(N-1,f,mSub,w); % the subslice filter
        % add them with the subslice phase
        b = bNotch + exp(1i*phi)*bSub;
            
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
        mNotch = [0 0 1 1 0 0 1 1 0 0];
        mSub = [0 0 0 0 1 1 0 0 0 0];
        
        w = [d1/d2 1 1 1 d1/d2];
    elseif Gind == 1
        % the slab and slice share a left transition band
        f = [0 shift-(1+ftw)*(tb/2) shift-(1-ftw)*(tb/2) ...
            Gcent+(tb/G/2-ftw*(tb/2)) Gcent+(tb/G/2+ftw*(tb/2)) ...
            shift+(1-ftw)*(tb/2) shift+(1+ftw)*(tb/2) (N/2)]/(N/2);
        mNotch = [0 0 0 0 1 1 0 0];
        mSub = [0 0 1 1 0 0 0 0];
        w = [d1/d2 1 1 d1/d2];
    elseif Gind == G
        % the slab and slice share a right transition band
        f = [0 shift-(1+ftw)*(tb/2) shift-(1-ftw)*(tb/2) ...
            Gcent-(tb/G/2+ftw*(tb/2)) Gcent-(tb/G/2-ftw*(tb/2)) ...
            shift+(1-ftw)*(tb/2) shift+(1+ftw)*(tb/2) (N/2)]/(N/2);
        mNotch = [0 0 1 1 0 0 0 0];
        mSub = [0 0 0 0 1 1 0 0];
        w = [d1/d2 1 1 d1/d2];
    end

    bNotch = firls(N-1,f,mNotch,w); % the notched filter
    % hilbert transform to suppress negative passband, and then demod to DC
    bNotch = hilbert(bNotch).*exp(-1i*2*pi/N*shift*(0:N-1))/2*exp(-1i*pi/N*shift);
    bSub = firls(N-1,f,mSub,w); % the subslice filter
    % hilbert transform to suppress negative passband, and then demod to DC
    bSub = hilbert(bSub).*exp(-1i*2*pi/N*shift*(0:N-1))/2*exp(-1i*pi/N*shift);
    % add them with the subslice phase
    b = bNotch + exp(1i*phi)*bSub;

end
