function [rfw,bw] = rootFlip(b,d1,flip,tbw)

% root flip the envelope to minimize peak B1/SAR
b = b./max(abs(freqz(b)));
b = b*sin(flip/2 + atan(d1*2)/2); % scale to target flip angle
r = roots(b); % calculate roots of single-band beta polynomial
% sort the roots by phase angle
[~,idx] = sort(angle(r));
r = r(idx);r = r(:);
r = leja_fast(r);

candidates = abs(1-abs(r)) > 0.004 & abs(angle(r)) < tbw/length(b)*pi;

maxRF = Inf;
for ii = 1:2^sum(candidates)
    
    % convert this ii to binary pattern
    doFlipStr = fliplr(dec2bin(ii-1));
    doFlip = zeros(1,length(doFlipStr));
    for jj = 1:length(doFlipStr)
       doFlip(jj) = str2num(doFlipStr(jj)); 
    end
    % add zeros for the rest of the passbands
    doFlip = [doFlip(:); zeros(sum(candidates)-length(doFlip),1)];
    % embed in all-roots vector
    tmp = zeros(length(b)-1,1);
    tmp(candidates) = doFlip;
    doFlip = tmp;
    % flip those indices
    rt = r(:);
    rt(doFlip == 1) = conj(1./rt(doFlip == 1));
    % get polynomial coefficients back from flipped roots
    tmp = poly(rt);
    % important: normalize before modulating! doesn't work
    % for negative shifted slices otherwise
    tmp = tmp./max(abs(freqz(tmp)));
    % store the rf
    bt = tmp(:)*sin(flip/2 + atan(d1*2)/2);
    rft = b2rf(bt);
    % take RF with min peak amplitude
    if max(abs(rft)) < maxRF
        maxRF = max(abs(rft));
        rfw = rft;
        bw = bt;
    end
    
end
