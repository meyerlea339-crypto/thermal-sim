// sehr einfache FFT Lib
export function fft(arr) {
  const N = arr.length;
  if (N <= 1) return [{ re: arr[0], im: 0 }];

  const even = fft(arr.filter((_, i) => i % 2 === 0));
  const odd  = fft(arr.filter((_, i) => i % 2 === 1));

  const out = new Array(N);

  for (let k = 0; k < N / 2; k++) {
    const t = -2 * Math.PI * k / N;
    const exp = { re: Math.cos(t), im: Math.sin(t) };
    const oddk = {
      re: exp.re * odd[k].re - exp.im * odd[k].im,
      im: exp.re * odd[k].im + exp.im * odd[k].re
    };

    out[k] = {
      re: even[k].re + oddk.re,
      im: even[k].im + oddk.im
    };
    out[k + N/2] = {
      re: even[k].re - oddk.re,
      im: even[k].im - oddk.im
    };
  }

  return out;
}

export function ifft(arr) {
  const N = arr.length;
  const conj = arr.map(v => ({ re: v.re, im: -v.im }));
  const t = fft(conj).map(v => ({ re: v.re / N, im: -v.im / N }));
  return t;
}
