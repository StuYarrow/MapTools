function ret = fftconv2(h, x)

fftSize = max([size(h) ; size(x)]);
outSize = abs(size(h) - size(x)) + 1;
pad = (fftSize - outSize) ./ 2;

H = fft2(h, fftSize(1), fftSize(2));
X = fft2(x, fftSize(1), fftSize(2));

HX = H .* X;

hx = real(ifft2(HX));

ret = hx(pad(1)+1 : end-pad(1), pad(2)+1 : end-pad(2));