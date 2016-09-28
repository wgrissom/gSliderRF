function b = gSliderBeta_cvx(N,G,Gind,tb,d1,d2,phi,osfact)

% Script to design a gSlider beta filter using cvx.
ftw = dinf(d1,d2)/tb; % fractional transition width of the slab profile

% set up target patterns and weights
nn = (-N/2*osfact:N/2*osfact-1)'/(N/2*osfact);  % -1 to 1 - profile indices
d = zeros(N*osfact,1);          % target pattern
s = nn <= -(1+ftw)*(tb/2)/(N/2) | nn >= (1+ftw)*(tb/2)/(N/2); % stopband mask
wts = 1./abs(nn).^0.5; % decaying ripple weights
wts = wts.*s;wts(isnan(wts)) = 0;wts = wts./max(wts);

Gcent = (Gind-G/2-1/2)*tb/G/(N/2);
if Gind > 1 && Gind < G
  % separate transition bands for slab+slice
  d = d + double(nn >= -(1-ftw)*(tb/2)/(N/2) & nn <= Gcent-(tb/G/2+ftw*(tb/2))/(N/2));
  d = d + exp(1i*phi)*double(nn >= Gcent-(tb/G/2-ftw*(tb/2))/(N/2) & nn <= Gcent+(tb/G/2-ftw*(tb/2))/(N/2));
  d = d + double(nn >= Gcent+(tb/G/2+ftw*(tb/2))/(N/2) & nn <= (1-ftw)*(tb/2)/(N/2));
elseif Gind == 1
  
elseif Gind == G
  
end

dd = d(s | logical(abs(d)));
wtsd = wts(s | logical(abs(d)));

ddri = [real(dd);imag(dd)];
ddriabs = [abs(dd);abs(dd)];
wtsdri = [wtsd;wtsd];
A = exp(-1i*2*pi/N*(-N/2:1/osfact:N/2-1/osfact)'*(-N/2:N/2-1));
Ad = A(s | logical(abs(d)),:);
Adri = [real(Ad) -imag(Ad);imag(Ad) real(Ad)]; % split to real+imaginary for cvx

% use cvx to do the constrained optimization
cvx_begin
    variable delta(1)
    variable x(2*N)
    minimize( delta )
    subject to
        -delta*ddriabs <= Adri*x - ddri <= delta*ddriabs + delta*d2/d1*wtsdri
cvx_end

% get complex filter back
b = x(1:N)+1i*x(N+1:end);
