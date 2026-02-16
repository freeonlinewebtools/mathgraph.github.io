/* ============================================================
   NumberGraph — Math Engine v2.0
   Number Theory + Calculus + Complex Numbers + Expression Parser
   ============================================================ */

const MathEngine = (() => {
    'use strict';

    // ========================================================
    // PRIME SIEVE & INFRASTRUCTURE
    // ========================================================
    const SIEVE_LIMIT = 1_100_000;
    let sieveBits = null;
    let primesList = [];
    let primeSet = null;

    function initSieve() {
        if (sieveBits) return;
        sieveBits = new Uint8Array(SIEVE_LIMIT + 1);
        sieveBits[0] = sieveBits[1] = 1;
        for (let i = 2; i * i <= SIEVE_LIMIT; i++) {
            if (!sieveBits[i]) {
                for (let j = i * i; j <= SIEVE_LIMIT; j += i) {
                    sieveBits[j] = 1;
                }
            }
        }
        primesList = [];
        for (let i = 2; i <= SIEVE_LIMIT; i++) {
            if (!sieveBits[i]) primesList.push(i);
        }
        primeSet = new Set(primesList);
    }
    initSieve();

    // ========================================================
    // COMPLEX NUMBER CLASS
    // ========================================================
    class Complex {
        constructor(re, im = 0) {
            this.re = re;
            this.im = im;
        }
        static from(v) {
            if (v instanceof Complex) return v;
            if (typeof v === 'number') return new Complex(v, 0);
            return new Complex(NaN, NaN);
        }
        add(other) {
            const b = Complex.from(other);
            return new Complex(this.re + b.re, this.im + b.im);
        }
        sub(other) {
            const b = Complex.from(other);
            return new Complex(this.re - b.re, this.im - b.im);
        }
        mul(other) {
            const b = Complex.from(other);
            return new Complex(
                this.re * b.re - this.im * b.im,
                this.re * b.im + this.im * b.re
            );
        }
        div(other) {
            const b = Complex.from(other);
            const d = b.re * b.re + b.im * b.im;
            if (d === 0) return new Complex(NaN, NaN);
            return new Complex(
                (this.re * b.re + this.im * b.im) / d,
                (this.im * b.re - this.re * b.im) / d
            );
        }
        neg() { return new Complex(-this.re, -this.im); }
        conj() { return new Complex(this.re, -this.im); }
        abs() { return Math.sqrt(this.re * this.re + this.im * this.im); }
        arg() { return Math.atan2(this.im, this.re); }

        pow(other) {
            const b = Complex.from(other);
            if (this.re === 0 && this.im === 0) {
                if (b.re > 0) return new Complex(0, 0);
                return new Complex(NaN, NaN);
            }
            const r = this.abs();
            const theta = this.arg();
            const lnr = Math.log(r);
            // z^w = exp(w * ln(z)), ln(z) = ln|z| + i*arg(z)
            const newRe = b.re * lnr - b.im * theta;
            const newIm = b.im * lnr + b.re * theta;
            const expRe = Math.exp(newRe);
            return new Complex(expRe * Math.cos(newIm), expRe * Math.sin(newIm));
        }

        static exp(z) {
            z = Complex.from(z);
            const r = Math.exp(z.re);
            return new Complex(r * Math.cos(z.im), r * Math.sin(z.im));
        }
        static log(z) {
            z = Complex.from(z);
            return new Complex(Math.log(z.abs()), z.arg());
        }
        static sin(z) {
            z = Complex.from(z);
            return new Complex(
                Math.sin(z.re) * Math.cosh(z.im),
                Math.cos(z.re) * Math.sinh(z.im)
            );
        }
        static cos(z) {
            z = Complex.from(z);
            return new Complex(
                Math.cos(z.re) * Math.cosh(z.im),
                -Math.sin(z.re) * Math.sinh(z.im)
            );
        }
        static tan(z) {
            return Complex.sin(z).div(Complex.cos(z));
        }
        static sqrt(z) {
            z = Complex.from(z);
            const r = z.abs();
            const newR = Math.sqrt(r);
            const theta = z.arg() / 2;
            return new Complex(newR * Math.cos(theta), newR * Math.sin(theta));
        }
        static sinh(z) {
            z = Complex.from(z);
            const ep = Complex.exp(z);
            const en = Complex.exp(z.neg());
            return ep.sub(en).mul(new Complex(0.5));
        }
        static cosh(z) {
            z = Complex.from(z);
            const ep = Complex.exp(z);
            const en = Complex.exp(z.neg());
            return ep.add(en).mul(new Complex(0.5));
        }
        static gamma(z) {
            z = Complex.from(z);
            // Lanczos approximation for complex gamma
            if (z.re < 0.5) {
                // Reflection: Γ(z) = π / (sin(πz) Γ(1-z))
                const oneMinusZ = new Complex(1 - z.re, -z.im);
                const sinPiZ = Complex.sin(new Complex(Math.PI * z.re, Math.PI * z.im));
                const gammaOneMinusZ = Complex.gamma(oneMinusZ);
                return new Complex(Math.PI, 0).div(sinPiZ.mul(gammaOneMinusZ));
            }
            z = z.sub(new Complex(1));
            const g = 7;
            const c = [
                0.99999999999980993, 676.5203681218851, -1259.1392167224028,
                771.32342877765313, -176.61502916214059, 12.507343278686905,
                -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7
            ];
            let x = new Complex(c[0]);
            for (let i = 1; i < g + 2; i++) {
                x = x.add(new Complex(c[i]).div(z.add(new Complex(i))));
            }
            const t = z.add(new Complex(g + 0.5));
            const sqrtTwoPi = new Complex(Math.sqrt(2 * Math.PI));
            return sqrtTwoPi.mul(t.pow(z.add(new Complex(0.5)))).mul(Complex.exp(t.neg())).mul(x);
        }
        static zeta(z) {
            z = Complex.from(z);
            // Dirichlet eta method: ζ(s) = η(s) / (1 - 2^(1-s))
            // η(s) = Σ (-1)^(n+1) / n^s using Euler transform for convergence
            const N = 40;
            let sum = new Complex(0);
            // Borwein's method coefficients
            const d = [];
            d[0] = 1;
            for (let k = 1; k <= N; k++) {
                d[k] = d[k-1] + N * (N + k - 1) * (N - k + 1) / ((2 * k - 1) * (2 * k) / 4) / d[k-1] * d[k-1]; // simplified
            }
            // Simpler: use Euler-Maclaurin
            for (let n = 1; n <= N; n++) {
                const nz = new Complex(n).pow(z.neg());
                sum = sum.add(nz);
            }
            // Euler-Maclaurin correction
            const Ns = new Complex(N).pow(new Complex(1).sub(z));
            const correction = Ns.div(z.sub(new Complex(1)));
            const half = new Complex(N).pow(z.neg()).mul(new Complex(0.5));
            sum = sum.add(correction).add(half);
            // B2 correction
            const B2 = 1/6;
            const nNs = new Complex(N).pow(z.neg());
            sum = sum.add(nNs.mul(z).mul(new Complex(B2 / N)));
            return sum;
        }

        isFinite() { return isFinite(this.re) && isFinite(this.im); }
        isReal() { return Math.abs(this.im) < 1e-12; }
        toReal() { return this.re; }
    }

    // ========================================================
    // NUMBER THEORY FUNCTIONS
    // ========================================================
    function isPrime(n) {
        n = Math.floor(n);
        if (n < 2) return 0;
        if (n <= SIEVE_LIMIT) return sieveBits[n] === 0 ? 1 : 0;
        if (n % 2 === 0) return 0;
        if (n % 3 === 0) return 0;
        for (let i = 5; i * i <= n; i += 6) {
            if (n % i === 0 || n % (i + 2) === 0) return 0;
        }
        return 1;
    }

    function nthPrime(n) {
        n = Math.floor(n);
        if (n < 1) return NaN;
        if (n <= primesList.length) return primesList[n - 1];
        return NaN;
    }

    function primeCount(x) {
        x = Math.floor(x);
        if (x < 2) return 0;
        if (x > SIEVE_LIMIT) x = SIEVE_LIMIT;
        let lo = 0, hi = primesList.length - 1;
        while (lo <= hi) {
            let mid = (lo + hi) >> 1;
            if (primesList[mid] <= x) lo = mid + 1;
            else hi = mid - 1;
        }
        return lo;
    }

    function nextPrime(n) {
        n = Math.floor(n) + 1;
        if (n < 2) return 2;
        while (n <= SIEVE_LIMIT) {
            if (!sieveBits[n]) return n;
            n++;
        }
        return NaN;
    }

    function prevPrime(n) {
        n = Math.floor(n) - 1;
        if (n < 2) return NaN;
        if (n > SIEVE_LIMIT) n = SIEVE_LIMIT;
        while (n >= 2) {
            if (!sieveBits[n]) return n;
            n--;
        }
        return NaN;
    }

    function primorial(n) {
        n = Math.floor(n);
        if (n < 2) return 1;
        let result = 1;
        for (const p of primesList) {
            if (p > n) break;
            result *= p;
            if (!isFinite(result)) return Infinity;
        }
        return result;
    }

    function twinPrimeCount(n) {
        n = Math.floor(n);
        if (n < 3) return 0;
        if (n > SIEVE_LIMIT) n = SIEVE_LIMIT;
        let count = 0;
        for (const p of primesList) {
            if (p > n - 2) break;
            if (p + 2 <= SIEVE_LIMIT && !sieveBits[p + 2]) count++;
        }
        return count;
    }

    function primeGap(n) {
        n = Math.floor(n);
        if (n < 2) return NaN;
        const np = nextPrime(n);
        if (isPrime(n)) return np - n;
        return NaN;
    }

    // --- Factorization ---
    function factorize(n) {
        n = Math.floor(Math.abs(n));
        if (n < 2) return [];
        const factors = [];
        for (let p = 2; p * p <= n; p++) {
            while (n % p === 0) { factors.push(p); n /= p; }
        }
        if (n > 1) factors.push(n);
        return factors;
    }

    function omega(n) {
        n = Math.floor(Math.abs(n));
        if (n < 2) return 0;
        const seen = new Set();
        for (let p = 2; p * p <= n; p++) {
            if (n % p === 0) { seen.add(p); while (n % p === 0) n /= p; }
        }
        if (n > 1) seen.add(n);
        return seen.size;
    }

    function bigOmega(n) { return factorize(n).length; }
    function greatestPrimeFactor(n) { const f = factorize(n); return f.length > 0 ? f[f.length - 1] : NaN; }
    function leastPrimeFactor(n) { const f = factorize(n); return f.length > 0 ? f[0] : NaN; }
    function radical(n) {
        n = Math.floor(Math.abs(n));
        if (n < 1) return NaN;
        if (n === 1) return 1;
        let result = 1;
        for (let p = 2; p * p <= n; p++) {
            if (n % p === 0) { result *= p; while (n % p === 0) n /= p; }
        }
        if (n > 1) result *= n;
        return result;
    }

    // --- Arithmetic functions ---
    function totient(n) {
        n = Math.floor(Math.abs(n));
        if (n < 1) return 0;
        if (n === 1) return 1;
        let result = n, temp = n;
        for (let p = 2; p * p <= temp; p++) {
            if (temp % p === 0) { while (temp % p === 0) temp /= p; result -= result / p; }
        }
        if (temp > 1) result -= result / temp;
        return Math.round(result);
    }

    function jordanTotient(n, k) {
        n = Math.floor(Math.abs(n));
        k = Math.floor(k);
        if (n < 1) return 0;
        if (n === 1) return 1;
        let result = Math.pow(n, k), temp = n;
        for (let p = 2; p * p <= temp; p++) {
            if (temp % p === 0) {
                while (temp % p === 0) temp /= p;
                result *= (1 - Math.pow(p, -k));
            }
        }
        if (temp > 1) result *= (1 - Math.pow(temp, -k));
        return Math.round(result);
    }

    function mobius(n) {
        n = Math.floor(Math.abs(n));
        if (n < 1) return 0;
        if (n === 1) return 1;
        let numFactors = 0;
        for (let p = 2; p * p <= n; p++) {
            if (n % p === 0) { n /= p; if (n % p === 0) return 0; numFactors++; }
        }
        if (n > 1) numFactors++;
        return (numFactors % 2 === 0) ? 1 : -1;
    }

    function sigma(n, k = 1) {
        n = Math.floor(Math.abs(n));
        k = Math.floor(k);
        if (n < 1) return 0;
        let sum = 0;
        for (let d = 1; d * d <= n; d++) {
            if (n % d === 0) { sum += Math.pow(d, k); if (d !== n / d) sum += Math.pow(n / d, k); }
        }
        return sum;
    }

    function numDivisors(n) { return sigma(n, 0); }

    function liouville(n) { return (bigOmega(n) % 2 === 0) ? 1 : -1; }

    function mangoldt(n) {
        n = Math.floor(Math.abs(n));
        if (n < 2) return 0;
        for (let p = 2; p * p <= n; p++) {
            if (n % p === 0) {
                let val = n;
                while (val % p === 0) val /= p;
                return val === 1 ? Math.log(p) : 0;
            }
        }
        return Math.log(n);
    }

    function mertens(n) {
        n = Math.floor(n);
        if (n < 1) return 0;
        let sum = 0;
        for (let k = 1; k <= Math.min(n, 100000); k++) sum += mobius(k);
        return sum;
    }

    function chebyshev1(x) {
        x = Math.floor(x);
        if (x < 2) return 0;
        let sum = 0;
        for (const p of primesList) { if (p > x) break; sum += Math.log(p); }
        return sum;
    }

    function chebyshev2(x) {
        x = Math.floor(x);
        if (x < 2) return 0;
        let sum = 0;
        for (let n = 2; n <= Math.min(x, SIEVE_LIMIT); n++) sum += mangoldt(n);
        return sum;
    }

    function divisors(n) {
        n = Math.floor(Math.abs(n));
        if (n < 1) return [];
        const divs = [];
        for (let d = 1; d * d <= n; d++) {
            if (n % d === 0) { divs.push(d); if (d !== n / d) divs.push(n / d); }
        }
        return divs.sort((a, b) => a - b);
    }

    function isperfect(n) {
        n = Math.floor(Math.abs(n));
        return n > 1 && sigma(n, 1) === 2 * n ? 1 : 0;
    }

    function isAbundant(n) {
        n = Math.floor(Math.abs(n));
        return n > 1 && sigma(n, 1) > 2 * n ? 1 : 0;
    }

    function isDeficient(n) {
        n = Math.floor(Math.abs(n));
        return n > 1 && sigma(n, 1) < 2 * n ? 1 : 0;
    }

    function digitsum(n, base = 10) {
        n = Math.floor(Math.abs(n));
        base = Math.floor(base);
        if (base < 2) return NaN;
        let sum = 0;
        while (n > 0) { sum += n % base; n = Math.floor(n / base); }
        return sum;
    }

    function digitalroot(n) {
        n = Math.floor(Math.abs(n));
        if (n === 0) return 0;
        return 1 + ((n - 1) % 9);
    }

    function collatz(n) {
        n = Math.floor(Math.abs(n));
        if (n < 1) return NaN;
        let steps = 0;
        while (n !== 1 && steps < 10000) {
            n = (n % 2 === 0) ? n / 2 : 3 * n + 1;
            steps++;
        }
        return steps;
    }

    // ========================================================
    // ANALYTIC & SPECIAL FUNCTIONS
    // ========================================================
    function gammaFunc(z) {
        if (z < 0.5) return Math.PI / (Math.sin(Math.PI * z) * gammaFunc(1 - z));
        z -= 1;
        const g = 7;
        const c = [
            0.99999999999980993, 676.5203681218851, -1259.1392167224028,
            771.32342877765313, -176.61502916214059, 12.507343278686905,
            -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7
        ];
        let x = c[0];
        for (let i = 1; i < g + 2; i++) x += c[i] / (z + i);
        const t = z + g + 0.5;
        return Math.sqrt(2 * Math.PI) * Math.pow(t, z + 0.5) * Math.exp(-t) * x;
    }

    function lnGamma(x) {
        if (x <= 0) return NaN;
        return Math.log(gammaFunc(x));
    }

    function betaFunc(a, b) {
        return gammaFunc(a) * gammaFunc(b) / gammaFunc(a + b);
    }

    function digamma(x) {
        // Psi function (digamma) - numerical approximation
        if (x <= 0 && x === Math.floor(x)) return NaN;
        let result = 0;
        // Use recurrence to shift x > 6
        while (x < 6) { result -= 1 / x; x += 1; }
        // Asymptotic series
        result += Math.log(x) - 1 / (2 * x);
        const x2 = x * x;
        result -= 1 / (12 * x2);
        result += 1 / (120 * x2 * x2);
        result -= 1 / (252 * x2 * x2 * x2);
        return result;
    }

    function logIntegral(x) {
        if (x <= 2) return 0;
        const steps = Math.max(100, Math.min(2000, Math.floor(x)));
        const h = (x - 2) / steps;
        let sum = 1 / Math.log(2) + 1 / Math.log(x);
        for (let i = 1; i < steps; i++) {
            const t = 2 + i * h;
            sum += ((i % 2 === 0) ? 2 : 4) / Math.log(t);
        }
        return (h / 3) * sum;
    }

    function zetaReal(s) {
        if (Math.abs(s - 1) < 1e-10) return NaN;
        if (s < 0) {
            const reflected = zetaReal(1 - s);
            return Math.pow(2, s) * Math.pow(Math.PI, s - 1) *
                Math.sin(Math.PI * s / 2) * gammaFunc(1 - s) * reflected;
        }
        if (s === 0) return -0.5;
        const N = 50;
        let sum = 0;
        for (let n = 1; n <= N; n++) sum += Math.pow(n, -s);
        sum += Math.pow(N, 1 - s) / (s - 1) + 0.5 * Math.pow(N, -s);
        const B2 = 1 / 6, B4 = -1 / 30;
        let Nns = Math.pow(N, -s);
        sum += B2 * s * Nns / N;
        sum += B4 * s * (s + 1) * (s + 2) * Nns / (24 * N * N * N);
        return sum;
    }

    function harmonic(n) {
        n = Math.floor(n);
        if (n < 1) return 0;
        if (n > 100000) {
            // Asymptotic: H_n ≈ ln(n) + γ + 1/(2n) - 1/(12n²)
            return Math.log(n) + 0.5772156649015329 + 1 / (2 * n) - 1 / (12 * n * n);
        }
        let sum = 0;
        for (let k = 1; k <= n; k++) sum += 1 / k;
        return sum;
    }

    function generalizedHarmonic(n, m) {
        n = Math.floor(n);
        if (n < 1) return 0;
        let sum = 0;
        for (let k = 1; k <= Math.min(n, 100000); k++) sum += Math.pow(k, -m);
        return sum;
    }

    // Error function
    function erf(x) {
        const a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741;
        const a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;
        const sign = x < 0 ? -1 : 1;
        x = Math.abs(x);
        const t = 1.0 / (1.0 + p * x);
        const y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
        return sign * y;
    }

    function erfc(x) { return 1 - erf(x); }

    // Lambert W function (principal branch)
    function lambertW(x) {
        if (x < -1 / Math.E) return NaN;
        if (x === 0) return 0;
        let w = x > 1 ? Math.log(x) - Math.log(Math.log(x)) : 0;
        for (let i = 0; i < 50; i++) {
            const ew = Math.exp(w);
            const wew = w * ew;
            const delta = (wew - x) / (ew * (w + 1) - (w + 2) * (wew - x) / (2 * w + 2));
            w -= delta;
            if (Math.abs(delta) < 1e-14) break;
        }
        return w;
    }

    // Exponential integral Ei(x)
    function expIntegral(x) {
        if (x <= 0) return NaN;
        if (x < 40) {
            let sum = 0.5772156649015329 + Math.log(x);
            let term = x;
            for (let n = 1; n < 100; n++) {
                sum += term / (n * n);
                term *= x / (n + 1);
                if (Math.abs(term / (n * n)) < 1e-15) break;
            }
            return sum;
        }
        // Asymptotic for large x
        let sum = 0, term = 1;
        for (let n = 1; n <= 20; n++) {
            term *= n / x;
            sum += term;
        }
        return Math.exp(x) / x * (1 + sum);
    }

    // Sine integral Si(x)
    function sineIntegral(x) {
        if (x === 0) return 0;
        const sign = x < 0 ? -1 : 1;
        x = Math.abs(x);
        // Simpson's rule
        const n = Math.max(100, Math.floor(x * 10));
        const h = x / n;
        let sum = Math.sin(0 + 1e-15) / (0 + 1e-15) + Math.sin(x) / x;
        for (let i = 1; i < n; i++) {
            const t = i * h;
            sum += ((i % 2 === 0) ? 2 : 4) * Math.sin(t) / t;
        }
        return sign * (h / 3) * sum;
    }

    // Cosine integral Ci(x)
    function cosineIntegral(x) {
        if (x <= 0) return NaN;
        // Ci(x) = γ + ln(x) + ∫₀ˣ (cos(t)-1)/t dt
        const n = Math.max(100, Math.floor(x * 10));
        const h = x / n;
        let sum = 0;
        for (let i = 1; i <= n; i++) {
            const t = i * h;
            const f = (Math.cos(t) - 1) / t;
            if (i === 1 || i === n) sum += f;
            else sum += ((i % 2 === 0) ? 2 : 4) * f;
        }
        return 0.5772156649015329 + Math.log(x) + (h / 3) * sum;
    }

    // Fresnel integrals
    function fresnelS(x) {
        const n = Math.max(200, Math.floor(Math.abs(x) * 20));
        const h = x / n;
        let sum = Math.sin(Math.PI / 2 * 0 * 0);
        sum += Math.sin(Math.PI / 2 * x * x);
        for (let i = 1; i < n; i++) {
            const t = i * h;
            sum += ((i % 2 === 0) ? 2 : 4) * Math.sin(Math.PI / 2 * t * t);
        }
        return (h / 3) * sum;
    }

    function fresnelC(x) {
        const n = Math.max(200, Math.floor(Math.abs(x) * 20));
        const h = x / n;
        let sum = Math.cos(Math.PI / 2 * 0 * 0);
        sum += Math.cos(Math.PI / 2 * x * x);
        for (let i = 1; i < n; i++) {
            const t = i * h;
            sum += ((i % 2 === 0) ? 2 : 4) * Math.cos(Math.PI / 2 * t * t);
        }
        return (h / 3) * sum;
    }

    // Bessel functions J0, J1, Jn (first kind)
    function besselJ(n, x) {
        n = Math.floor(n);
        // Series expansion: Jn(x) = Σ (-1)^k / (k! Γ(n+k+1)) (x/2)^(n+2k)
        let sum = 0;
        for (let k = 0; k < 50; k++) {
            const sign = (k % 2 === 0) ? 1 : -1;
            const term = sign / (factorial(k) * gammaFunc(n + k + 1)) * Math.pow(x / 2, n + 2 * k);
            sum += term;
            if (Math.abs(term) < 1e-15 * Math.abs(sum)) break;
        }
        return sum;
    }

    // Bessel Y0 (second kind, order 0)
    function besselY(n, x) {
        if (x <= 0) return NaN;
        // Use cross product relation: Yn(x) ≈ (Jn(x)cos(nπ) - J(-n)(x)) / sin(nπ) for non-integer n
        // For integer n, use series
        n = Math.floor(n);
        const eps = 0.001;
        const jn_plus = besselJ(n + eps, x);
        const jn_minus = besselJ(n - eps, x);
        // Numerical approximation via finite difference of Bessel J
        // Actually, let's use: Y_n(x) ≈ (2/π) [J_n(x)(γ + ln(x/2)) - sum...]
        // Simpler numerical approach:
        const h = 0.0001;
        return (besselJ(n + h, x) * Math.cos((n + h) * Math.PI) - besselJ(-(n + h), x)) / Math.sin((n + h) * Math.PI);
    }

    // Airy function Ai(x)
    function airyAi(x) {
        // Use series for small x, asymptotic for large x
        if (Math.abs(x) < 5) {
            let sum1 = 0, sum2 = 0;
            for (let k = 0; k < 30; k++) {
                sum1 += Math.pow(x, 3 * k) / (Math.pow(3, 2 * k / 3) * factorial(k) * gammaFunc(k + 2 / 3));
                sum2 += Math.pow(x, 3 * k + 1) / (Math.pow(3, (2 * k + 1) / 3) * factorial(k) * gammaFunc(k + 4 / 3));
            }
            return sum1 / (Math.pow(3, 2/3) * gammaFunc(2/3)) - sum2 / (Math.pow(3, 1/3) * gammaFunc(1/3));
        }
        if (x > 0) {
            const xi = 2/3 * Math.pow(x, 1.5);
            return Math.exp(-xi) / (2 * Math.sqrt(Math.PI) * Math.pow(x, 0.25));
        }
        const xi = 2/3 * Math.pow(-x, 1.5);
        return Math.sin(xi + Math.PI/4) / (Math.sqrt(Math.PI) * Math.pow(-x, 0.25));
    }

    // Sinc function
    function sinc(x) {
        if (Math.abs(x) < 1e-10) return 1;
        return Math.sin(Math.PI * x) / (Math.PI * x);
    }

    // Heaviside step function
    function heaviside(x) { return x < 0 ? 0 : (x === 0 ? 0.5 : 1); }

    // --- Sequences ---
    const fibCache = new Map([[0, 0], [1, 1]]);
    function fibonacci(n) {
        n = Math.floor(n);
        if (n < 0) return NaN;
        if (fibCache.has(n)) return fibCache.get(n);
        const val = fibonacci(n - 1) + fibonacci(n - 2);
        fibCache.set(n, val);
        return val;
    }

    const lucasCache = new Map([[0, 2], [1, 1]]);
    function lucas(n) {
        n = Math.floor(n);
        if (n < 0) return NaN;
        if (lucasCache.has(n)) return lucasCache.get(n);
        const val = lucas(n - 1) + lucas(n - 2);
        lucasCache.set(n, val);
        return val;
    }

    function catalan(n) {
        n = Math.floor(n);
        if (n < 0) return NaN;
        return binomial(2 * n, n) / (n + 1);
    }

    const partitionCache = new Map([[0, 1]]);
    function partition(n) {
        n = Math.floor(n);
        if (n < 0) return 0;
        if (partitionCache.has(n)) return partitionCache.get(n);
        let sum = 0;
        for (let k = 1; ; k++) {
            const g1 = k * (3 * k - 1) / 2;
            const g2 = k * (3 * k + 1) / 2;
            if (g1 > n && g2 > n) break;
            const sign = (k % 2 === 1) ? 1 : -1;
            if (g1 <= n) sum += sign * partition(n - g1);
            if (g2 <= n) sum += sign * partition(n - g2);
        }
        partitionCache.set(n, sum);
        return sum;
    }

    function bernoulli(n) {
        n = Math.floor(n);
        if (n < 0) return NaN;
        if (n === 0) return 1;
        if (n === 1) return -0.5;
        if (n % 2 === 1) return 0; // B_n = 0 for odd n > 1
        const B = [1, -0.5];
        for (let m = 2; m <= n; m++) {
            if (m % 2 === 1 && m > 1) { B.push(0); continue; }
            let sum = 0;
            for (let k = 0; k < m; k++) {
                sum += binomial(m + 1, k) * B[k];
            }
            B.push(-sum / (m + 1));
        }
        return B[n];
    }

    function bell(n) {
        n = Math.floor(n);
        if (n < 0) return NaN;
        if (n === 0) return 1;
        // Bell triangle
        let prev = [1];
        for (let i = 1; i <= n; i++) {
            const curr = [prev[prev.length - 1]];
            for (let j = 1; j <= i; j++) {
                curr.push(curr[j - 1] + prev[j - 1]);
            }
            prev = curr;
        }
        return prev[0];
    }

    // --- Utility ---
    function gcd(a, b) {
        a = Math.floor(Math.abs(a));
        b = Math.floor(Math.abs(b));
        while (b) { [a, b] = [b, a % b]; }
        return a;
    }

    function lcm(a, b) {
        a = Math.floor(Math.abs(a));
        b = Math.floor(Math.abs(b));
        return (a === 0 || b === 0) ? 0 : (a / gcd(a, b)) * b;
    }

    function factorial(n) {
        n = Math.floor(n);
        if (n < 0) return gammaFunc(n + 1); // extend via gamma
        if (n <= 1) return 1;
        if (n > 170) return Infinity;
        let result = 1;
        for (let i = 2; i <= n; i++) result *= i;
        return result;
    }

    function binomial(n, k) {
        n = Math.floor(n);
        k = Math.floor(k);
        if (k < 0 || k > n) return 0;
        if (k === 0 || k === n) return 1;
        if (k > n - k) k = n - k;
        let result = 1;
        for (let i = 0; i < k; i++) result = result * (n - i) / (i + 1);
        return Math.round(result);
    }

    function isqrt(n) { return Math.floor(Math.sqrt(Math.floor(Math.abs(n)))); }

    function stirling1(n, k) {
        // Unsigned Stirling numbers of the first kind
        n = Math.floor(n); k = Math.floor(k);
        if (n < 0 || k < 0 || k > n) return 0;
        if (n === 0 && k === 0) return 1;
        if (n === 0 || k === 0) return 0;
        // Recurrence: |s(n,k)| = (n-1)*|s(n-1,k)| + |s(n-1,k-1)|
        // Use memoized table
        const s = Array.from({ length: n + 1 }, () => new Array(k + 1).fill(0));
        s[0][0] = 1;
        for (let i = 1; i <= n; i++) {
            for (let j = 1; j <= Math.min(i, k); j++) {
                s[i][j] = (i - 1) * s[i - 1][j] + s[i - 1][j - 1];
            }
        }
        return s[n][k];
    }

    function stirling2(n, k) {
        // Stirling numbers of the second kind
        n = Math.floor(n); k = Math.floor(k);
        if (n < 0 || k < 0 || k > n) return 0;
        if (n === 0 && k === 0) return 1;
        if (n === 0 || k === 0) return 0;
        let sum = 0;
        for (let j = 0; j <= k; j++) {
            const sign = ((k - j) % 2 === 0) ? 1 : -1;
            sum += sign * binomial(k, j) * Math.pow(j, n);
        }
        return Math.round(sum / factorial(k));
    }

    // Polygamma function
    function polygamma(m, x) {
        m = Math.floor(m);
        if (m === 0) return digamma(x);
        // ψ^(m)(x) = (-1)^(m+1) m! Σ 1/(x+k)^(m+1)
        const sign = (m % 2 === 0) ? -1 : 1;
        let sum = 0;
        for (let k = 0; k < 200; k++) {
            sum += Math.pow(x + k, -(m + 1));
            if (Math.abs(Math.pow(x + k, -(m + 1))) < 1e-15 * Math.abs(sum)) break;
        }
        return sign * factorial(m) * sum;
    }

    // Dedekind eta function (on the imaginary axis)
    function dedekindEta(tau_im) {
        // η(iτ) = e^(-π τ/12) Π (1 - e^(-2π n τ))
        if (tau_im <= 0) return NaN;
        let prod = 1;
        const q = Math.exp(-2 * Math.PI * tau_im);
        let qn = q;
        for (let n = 1; n <= 100; n++) {
            prod *= (1 - qn);
            qn *= q;
            if (qn < 1e-30) break;
        }
        return Math.exp(-Math.PI * tau_im / 12) * prod;
    }

    // ========================================================
    // CALCULUS FUNCTIONS
    // ========================================================

    // Numerical derivative: evaluate AST with x = point ± h
    function numericalDerivative(ast, point, vars, variable = 'x') {
        const h = Math.max(1e-7, Math.abs(point) * 1e-7);
        const varsPlus = { ...vars, [variable]: point + h };
        const varsMinus = { ...vars, [variable]: point - h };
        const varsPlus2 = { ...vars, [variable]: point + 2 * h };
        const varsMinus2 = { ...vars, [variable]: point - 2 * h };
        try {
            // 5-point stencil for better accuracy
            const f1 = evaluate(ast, varsMinus2);
            const f2 = evaluate(ast, varsMinus);
            const f3 = evaluate(ast, varsPlus);
            const f4 = evaluate(ast, varsPlus2);
            return (-f4 + 8 * f3 - 8 * f2 + f1) / (12 * h);
        } catch {
            return NaN;
        }
    }

    // Second derivative
    function numericalSecondDerivative(ast, point, vars, variable = 'x') {
        const h = Math.max(1e-5, Math.abs(point) * 1e-5);
        const v0 = { ...vars, [variable]: point };
        const vp = { ...vars, [variable]: point + h };
        const vm = { ...vars, [variable]: point - h };
        try {
            return (evaluate(ast, vp) - 2 * evaluate(ast, v0) + evaluate(ast, vm)) / (h * h);
        } catch {
            return NaN;
        }
    }

    // Numerical integration (Simpson's 3/8 rule, adaptive)
    function numericalIntegral(ast, a, b, vars, variable = 'x') {
        if (a === b) return 0;
        if (a > b) return -numericalIntegral(ast, b, a, vars, variable);
        const steps = 200;
        const h = (b - a) / steps;
        let sum = 0;
        const evalAt = (t) => {
            const v = { ...vars, [variable]: t };
            try { return evaluate(ast, v); } catch { return 0; }
        };
        sum += evalAt(a) + evalAt(b);
        for (let i = 1; i < steps; i++) {
            const t = a + i * h;
            sum += ((i % 2 === 0) ? 2 : 4) * evalAt(t);
        }
        return (h / 3) * sum;
    }

    // Summation: Σ f(n) from n=a to n=b
    function computeSum(ast, varName, a, b, vars) {
        a = Math.floor(a);
        b = Math.floor(b);
        if (b - a > 100000) b = a + 100000;
        let sum = 0;
        for (let n = a; n <= b; n++) {
            const v = { ...vars, [varName]: n };
            try {
                const val = evaluate(ast, v);
                if (isFinite(val)) sum += val;
            } catch { /* skip */ }
        }
        return sum;
    }

    // Product: Π f(n) from n=a to n=b
    function computeProduct(ast, varName, a, b, vars) {
        a = Math.floor(a);
        b = Math.floor(b);
        if (b - a > 100000) b = a + 100000;
        let prod = 1;
        for (let n = a; n <= b; n++) {
            const v = { ...vars, [varName]: n };
            try {
                const val = evaluate(ast, v);
                if (isFinite(val)) prod *= val;
            } catch { /* skip */ }
        }
        return prod;
    }

    // ========================================================
    // FUNCTION REGISTRY
    // ========================================================
    const FUNCTIONS = {
        // -- Number theory --
        isprime:    { fn: (n) => isPrime(n), args: 1, desc: 'Primality test', discrete: true },
        prime:      { fn: (n) => nthPrime(n), args: 1, desc: 'Nth prime', discrete: true },
        primepi:    { fn: (x) => primeCount(x), args: 1, desc: 'Prime counting π(x)', discrete: false },
        nextprime:  { fn: (n) => nextPrime(n), args: 1, desc: 'Next prime', discrete: true },
        prevprime:  { fn: (n) => prevPrime(n), args: 1, desc: 'Previous prime', discrete: true },
        primorial:  { fn: (n) => primorial(n), args: 1, desc: 'Primorial', discrete: true },
        twinprimes: { fn: (n) => twinPrimeCount(n), args: 1, desc: 'Twin prime count', discrete: false },
        primegap:   { fn: (n) => primeGap(n), args: 1, desc: 'Gap to next prime', discrete: true },
        totient:    { fn: (n) => totient(n), args: 1, desc: "Euler's totient φ(n)", discrete: true },
        jordan:     { fn: (n, k) => jordanTotient(n, k), args: 2, desc: "Jordan's totient", discrete: true },
        mobius:     { fn: (n) => mobius(n), args: 1, desc: 'Möbius function μ(n)', discrete: true },
        sigma:      { fn: (n, k) => sigma(n, k === undefined ? 1 : k), args: [1, 2], desc: 'Divisor sum σ(n)', discrete: true },
        numdiv:     { fn: (n) => numDivisors(n), args: 1, desc: 'Number of divisors d(n)', discrete: true },
        liouville:  { fn: (n) => liouville(n), args: 1, desc: 'Liouville function', discrete: true },
        mangoldt:   { fn: (n) => mangoldt(n), args: 1, desc: 'Von Mangoldt Λ(n)', discrete: true },
        omega:      { fn: (n) => omega(n), args: 1, desc: 'Distinct prime factors ω(n)', discrete: true },
        bigomega:   { fn: (n) => bigOmega(n), args: 1, desc: 'Prime factors Ω(n)', discrete: true },
        mertens:    { fn: (n) => mertens(n), args: 1, desc: 'Mertens function M(n)', discrete: true },
        chebyshev1: { fn: (x) => chebyshev1(x), args: 1, desc: 'Chebyshev θ(x)', discrete: false },
        chebyshev2: { fn: (x) => chebyshev2(x), args: 1, desc: 'Chebyshev ψ(x)', discrete: false },
        li:         { fn: (x) => logIntegral(x), args: 1, desc: 'Log integral Li(x)', discrete: false },
        zeta:       { fn: (s) => zetaReal(s), args: 1, desc: 'Riemann zeta ζ(s)', discrete: false },
        harmonic:   { fn: (n) => harmonic(n), args: 1, desc: 'Harmonic H(n)', discrete: true },
        gharmonic:  { fn: (n, m) => generalizedHarmonic(n, m), args: 2, desc: 'Generalized harmonic', discrete: true },
        radical:    { fn: (n) => radical(n), args: 1, desc: 'Radical of n', discrete: true },
        gpf:        { fn: (n) => greatestPrimeFactor(n), args: 1, desc: 'Greatest prime factor', discrete: true },
        lpf:        { fn: (n) => leastPrimeFactor(n), args: 1, desc: 'Least prime factor', discrete: true },
        isperfect:  { fn: (n) => isperfect(n), args: 1, desc: 'Perfect number test', discrete: true },
        isabundant: { fn: (n) => isAbundant(n), args: 1, desc: 'Abundant number test', discrete: true },
        isdeficient:{ fn: (n) => isDeficient(n), args: 1, desc: 'Deficient number test', discrete: true },
        digitsum:   { fn: (n, b) => digitsum(n, b === undefined ? 10 : b), args: [1, 2], desc: 'Digit sum', discrete: true },
        digitalroot:{ fn: (n) => digitalroot(n), args: 1, desc: 'Digital root', discrete: true },
        collatz:    { fn: (n) => collatz(n), args: 1, desc: 'Collatz steps', discrete: true },

        // -- Sequences --
        fibonacci:  { fn: (n) => fibonacci(n), args: 1, desc: 'Fibonacci F(n)', discrete: true },
        fib:        { fn: (n) => fibonacci(n), args: 1, desc: 'Fibonacci F(n)', discrete: true },
        lucas:      { fn: (n) => lucas(n), args: 1, desc: 'Lucas numbers', discrete: true },
        catalan:    { fn: (n) => catalan(n), args: 1, desc: 'Catalan numbers', discrete: true },
        partition:  { fn: (n) => partition(n), args: 1, desc: 'Partition p(n)', discrete: true },
        bernoulli:  { fn: (n) => bernoulli(n), args: 1, desc: 'Bernoulli numbers', discrete: true },
        bell:       { fn: (n) => bell(n), args: 1, desc: 'Bell numbers', discrete: true },
        stirling1:  { fn: (n, k) => stirling1(n, k), args: 2, desc: 'Stirling 1st kind', discrete: true },
        stirling2:  { fn: (n, k) => stirling2(n, k), args: 2, desc: 'Stirling 2nd kind', discrete: true },

        // -- Special / Analytic --
        gamma:      { fn: (z) => gammaFunc(z), args: 1, desc: 'Gamma function Γ(z)' },
        lngamma:    { fn: (x) => lnGamma(x), args: 1, desc: 'Log-gamma ln Γ(x)' },
        beta:       { fn: (a, b) => betaFunc(a, b), args: 2, desc: 'Beta function B(a,b)' },
        digamma:    { fn: (x) => digamma(x), args: 1, desc: 'Digamma ψ(x)' },
        polygamma:  { fn: (m, x) => polygamma(m, x), args: 2, desc: 'Polygamma ψ^(m)(x)' },
        erf:        { fn: (x) => erf(x), args: 1, desc: 'Error function' },
        erfc:       { fn: (x) => erfc(x), args: 1, desc: 'Complementary error' },
        lambertw:   { fn: (x) => lambertW(x), args: 1, desc: 'Lambert W function' },
        ei:         { fn: (x) => expIntegral(x), args: 1, desc: 'Exponential integral' },
        si:         { fn: (x) => sineIntegral(x), args: 1, desc: 'Sine integral Si(x)' },
        ci:         { fn: (x) => cosineIntegral(x), args: 1, desc: 'Cosine integral Ci(x)' },
        fresnels:   { fn: (x) => fresnelS(x), args: 1, desc: 'Fresnel S(x)' },
        fresnelc:   { fn: (x) => fresnelC(x), args: 1, desc: 'Fresnel C(x)' },
        besselj:    { fn: (n, x) => besselJ(n, x), args: 2, desc: 'Bessel J_n(x)' },
        bessely:    { fn: (n, x) => besselY(n, x), args: 2, desc: 'Bessel Y_n(x)' },
        airy:       { fn: (x) => airyAi(x), args: 1, desc: 'Airy Ai(x)' },
        sinc:       { fn: (x) => sinc(x), args: 1, desc: 'Sinc function' },
        heaviside:  { fn: (x) => heaviside(x), args: 1, desc: 'Heaviside step' },
        eta:        { fn: (t) => dedekindEta(t), args: 1, desc: 'Dedekind eta' },

        // -- Utility --
        gcd:        { fn: (a, b) => gcd(a, b), args: 2, desc: 'GCD', discrete: true },
        lcm:        { fn: (a, b) => lcm(a, b), args: 2, desc: 'LCM', discrete: true },
        factorial:  { fn: (n) => factorial(n), args: 1, desc: 'Factorial n!' },
        binomial:   { fn: (n, k) => binomial(n, k), args: 2, desc: 'Binomial C(n,k)' },
        mod:        { fn: (a, b) => ((a % b) + b) % b, args: 2, desc: 'Modulo', discrete: true },
        isqrt:      { fn: (n) => isqrt(n), args: 1, desc: 'Integer sqrt', discrete: true },
        clamp:      { fn: (x, lo, hi) => Math.min(Math.max(x, lo), hi), args: 3, desc: 'Clamp value' },
        lerp:       { fn: (a, b, t) => a + (b - a) * t, args: 3, desc: 'Linear interpolation' },
        piecewise:  { fn: (...args) => {
            for (let i = 0; i < args.length - 1; i += 2) {
                if (args[i]) return args[i + 1];
            }
            return args.length % 2 === 1 ? args[args.length - 1] : 0;
        }, args: [2, 10], desc: 'Piecewise(cond,val,...)' },

        // -- Standard math --
        sin:    { fn: Math.sin, args: 1, desc: 'Sine' },
        cos:    { fn: Math.cos, args: 1, desc: 'Cosine' },
        tan:    { fn: Math.tan, args: 1, desc: 'Tangent' },
        sec:    { fn: (x) => 1 / Math.cos(x), args: 1, desc: 'Secant' },
        csc:    { fn: (x) => 1 / Math.sin(x), args: 1, desc: 'Cosecant' },
        cot:    { fn: (x) => Math.cos(x) / Math.sin(x), args: 1, desc: 'Cotangent' },
        asin:   { fn: Math.asin, args: 1, desc: 'Arcsine' },
        acos:   { fn: Math.acos, args: 1, desc: 'Arccosine' },
        atan:   { fn: Math.atan, args: 1, desc: 'Arctangent' },
        atan2:  { fn: Math.atan2, args: 2, desc: 'Arctangent (2-arg)' },
        asec:   { fn: (x) => Math.acos(1 / x), args: 1, desc: 'Arcsecant' },
        acsc:   { fn: (x) => Math.asin(1 / x), args: 1, desc: 'Arccosecant' },
        acot:   { fn: (x) => Math.atan(1 / x), args: 1, desc: 'Arccotangent' },
        sinh:   { fn: Math.sinh, args: 1, desc: 'Hyperbolic sine' },
        cosh:   { fn: Math.cosh, args: 1, desc: 'Hyperbolic cosine' },
        tanh:   { fn: Math.tanh, args: 1, desc: 'Hyperbolic tangent' },
        asinh:  { fn: Math.asinh, args: 1, desc: 'Inverse hyp sine' },
        acosh:  { fn: Math.acosh, args: 1, desc: 'Inverse hyp cosine' },
        atanh:  { fn: Math.atanh, args: 1, desc: 'Inverse hyp tangent' },
        sech:   { fn: (x) => 1 / Math.cosh(x), args: 1, desc: 'Hyperbolic secant' },
        csch:   { fn: (x) => 1 / Math.sinh(x), args: 1, desc: 'Hyperbolic cosecant' },
        coth:   { fn: (x) => Math.cosh(x) / Math.sinh(x), args: 1, desc: 'Hyperbolic cotangent' },
        ln:     { fn: Math.log, args: 1, desc: 'Natural logarithm' },
        log:    { fn: Math.log10, args: 1, desc: 'Log base 10' },
        log2:   { fn: Math.log2, args: 1, desc: 'Log base 2' },
        log10:  { fn: Math.log10, args: 1, desc: 'Log base 10' },
        logb:   { fn: (b, x) => Math.log(x) / Math.log(b), args: 2, desc: 'Log base b' },
        sqrt:   { fn: Math.sqrt, args: 1, desc: 'Square root' },
        cbrt:   { fn: Math.cbrt, args: 1, desc: 'Cube root' },
        nroot:  { fn: (n, x) => Math.pow(x, 1 / n), args: 2, desc: 'Nth root' },
        abs:    { fn: Math.abs, args: 1, desc: 'Absolute value' },
        sign:   { fn: Math.sign, args: 1, desc: 'Sign function' },
        floor:  { fn: Math.floor, args: 1, desc: 'Floor' },
        ceil:   { fn: Math.ceil, args: 1, desc: 'Ceiling' },
        round:  { fn: Math.round, args: 1, desc: 'Round' },
        trunc:  { fn: Math.trunc, args: 1, desc: 'Truncate' },
        frac:   { fn: (x) => x - Math.floor(x), args: 1, desc: 'Fractional part' },
        exp:    { fn: Math.exp, args: 1, desc: 'Exponential e^x' },
        min:    { fn: Math.min, args: 2, desc: 'Minimum' },
        max:    { fn: Math.max, args: 2, desc: 'Maximum' },
        pow:    { fn: Math.pow, args: 2, desc: 'Power' },
        hypot:  { fn: Math.hypot, args: 2, desc: 'Hypotenuse' },

        // -- Complex --
        re:     { fn: (x) => x, args: 1, desc: 'Real part' },
        im:     { fn: () => 0, args: 1, desc: 'Imaginary part' },
        carg:   { fn: (x) => x >= 0 ? 0 : Math.PI, args: 1, desc: 'Complex argument' },
        cabs:   { fn: Math.abs, args: 1, desc: 'Complex modulus' },
    };

    const CONSTANTS = {
        pi: Math.PI,
        e: Math.E,
        phi: (1 + Math.sqrt(5)) / 2,
        tau: 2 * Math.PI,
        eulergamma: 0.5772156649015329,
        catalan: 0.9159655941772190,
        apery: 1.2020569031595943,
        inf: Infinity,
        infinity: Infinity,
    };

    // ========================================================
    // TOKENIZER
    // ========================================================
    const TokenType = {
        NUMBER: 'NUMBER', IDENT: 'IDENT',
        PLUS: '+', MINUS: '-', STAR: '*', SLASH: '/',
        PERCENT: '%', CARET: '^', LPAREN: '(', RPAREN: ')',
        COMMA: ',', BANG: '!', PIPE: '|', EQUALS: '=',
        LT: '<', GT: '>', LTE: '<=', GTE: '>=', EQEQ: '==', NEQ: '!=',
        AND: '&&', OR: '||', SEMICOLON: ';', LBRACKET: '[', RBRACKET: ']',
        EOF: 'EOF',
    };

    function tokenize(input) {
        const tokens = [];
        let i = 0;
        const len = input.length;

        while (i < len) {
            if (input[i] === ' ' || input[i] === '\t') { i++; continue; }

            // Number
            if (isDigit(input[i]) || (input[i] === '.' && i + 1 < len && isDigit(input[i + 1]))) {
                let num = '';
                while (i < len && (isDigit(input[i]) || input[i] === '.')) num += input[i++];
                if (i < len && (input[i] === 'e' || input[i] === 'E')) {
                    num += input[i++];
                    if (i < len && (input[i] === '+' || input[i] === '-')) num += input[i++];
                    while (i < len && isDigit(input[i])) num += input[i++];
                }
                tokens.push({ type: TokenType.NUMBER, value: parseFloat(num) });
                continue;
            }

            // Identifier
            if (isAlpha(input[i]) || input[i] === '_') {
                let ident = '';
                while (i < len && (isAlpha(input[i]) || isDigit(input[i]) || input[i] === '_')) ident += input[i++];
                tokens.push({ type: TokenType.IDENT, value: ident });
                continue;
            }

            // Two-char operators
            if (i + 1 < len) {
                const two = input[i] + input[i + 1];
                if (two === '<=') { tokens.push({ type: TokenType.LTE }); i += 2; continue; }
                if (two === '>=') { tokens.push({ type: TokenType.GTE }); i += 2; continue; }
                if (two === '==') { tokens.push({ type: TokenType.EQEQ }); i += 2; continue; }
                if (two === '!=') { tokens.push({ type: TokenType.NEQ }); i += 2; continue; }
                if (two === '&&') { tokens.push({ type: TokenType.AND }); i += 2; continue; }
                if (two === '||') { tokens.push({ type: TokenType.OR }); i += 2; continue; }
            }

            // Single-char operators
            const opMap = {
                '+': TokenType.PLUS, '-': TokenType.MINUS, '*': TokenType.STAR,
                '/': TokenType.SLASH, '%': TokenType.PERCENT, '^': TokenType.CARET,
                '(': TokenType.LPAREN, ')': TokenType.RPAREN, ',': TokenType.COMMA,
                '!': TokenType.BANG, '|': TokenType.PIPE, '=': TokenType.EQUALS,
                '<': TokenType.LT, '>': TokenType.GT, ';': TokenType.SEMICOLON,
                '[': TokenType.LBRACKET, ']': TokenType.RBRACKET,
            };
            if (opMap[input[i]]) { tokens.push({ type: opMap[input[i]] }); i++; continue; }
            i++; // skip unknown
        }

        tokens.push({ type: TokenType.EOF });
        return tokens;
    }

    function isDigit(c) { return c >= '0' && c <= '9'; }
    function isAlpha(c) { return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z'); }

    // ========================================================
    // PARSER — Recursive Descent
    // ========================================================
    function parse(tokens) {
        let pos = 0;
        function peek() { return tokens[pos]; }
        function advance() { return tokens[pos++]; }
        function expect(type) {
            const tok = advance();
            if (tok.type !== type) throw new Error(`Expected ${type}, got ${tok.type}`);
            return tok;
        }
        function match(type) { if (peek().type === type) { advance(); return true; } return false; }
        function canStartTerm(tok) {
            return tok.type === TokenType.NUMBER || tok.type === TokenType.IDENT ||
                tok.type === TokenType.LPAREN || tok.type === TokenType.PIPE;
        }

        function expression() { return comparison(); }

        function comparison() {
            let left = additive();
            const compOps = [TokenType.LT, TokenType.GT, TokenType.LTE, TokenType.GTE, TokenType.EQEQ, TokenType.NEQ];
            while (compOps.includes(peek().type)) {
                const opTok = advance();
                const opMap = { '<': '<', '>': '>', '<=': '<=', '>=': '>=', '==': '==', '!=': '!=' };
                const op = opMap[opTok.type] || opTok.type;
                const right = additive();
                left = { type: 'comparison', op, left, right };
            }
            return left;
        }

        function additive() {
            let left = multiplicative();
            while (peek().type === TokenType.PLUS || peek().type === TokenType.MINUS) {
                const op = advance().type === TokenType.PLUS ? '+' : '-';
                left = { type: 'binary', op, left, right: multiplicative() };
            }
            return left;
        }

        function multiplicative() {
            let left = implicitMult();
            while (peek().type === TokenType.STAR || peek().type === TokenType.SLASH || peek().type === TokenType.PERCENT) {
                const tok = advance();
                const op = tok.type === TokenType.STAR ? '*' : (tok.type === TokenType.SLASH ? '/' : '%');
                left = { type: 'binary', op, left, right: implicitMult() };
            }
            return left;
        }

        function implicitMult() {
            let left = power();
            while (canStartTerm(peek())) {
                left = { type: 'binary', op: '*', left, right: power() };
            }
            return left;
        }

        function power() {
            let base = unary();
            if (match(TokenType.CARET)) {
                base = { type: 'binary', op: '^', left: base, right: power() };
            }
            return base;
        }

        function unary() {
            if (peek().type === TokenType.MINUS) { advance(); return { type: 'unary', op: '-', operand: postfix() }; }
            if (peek().type === TokenType.PLUS) { advance(); return postfix(); }
            return postfix();
        }

        function postfix() {
            let node = call();
            while (peek().type === TokenType.BANG) {
                advance();
                node = { type: 'call', name: 'factorial', args: [node] };
            }
            return node;
        }

        function call() {
            if (peek().type === TokenType.IDENT) {
                const name = peek().value;
                const lowerName = name.toLowerCase();

                // Special: sum(expr, var, start, end) and prod(expr, var, start, end)
                if ((lowerName === 'sum' || lowerName === 'prod') &&
                    pos + 1 < tokens.length && tokens[pos + 1].type === TokenType.LPAREN) {
                    advance(); advance(); // consume name and (
                    const exprAst = expression();
                    expect(TokenType.COMMA);
                    const varTok = advance();
                    if (varTok.type !== TokenType.IDENT) throw new Error('Expected variable name');
                    expect(TokenType.COMMA);
                    const startAst = expression();
                    expect(TokenType.COMMA);
                    const endAst = expression();
                    expect(TokenType.RPAREN);
                    return { type: 'special', name: lowerName, body: exprAst, varName: varTok.value.toLowerCase(), start: startAst, end: endAst };
                }

                // Special: deriv(expr) or diff(expr) — numerical derivative w.r.t. x
                if ((lowerName === 'deriv' || lowerName === 'diff' || lowerName === 'derivative') &&
                    pos + 1 < tokens.length && tokens[pos + 1].type === TokenType.LPAREN) {
                    advance(); advance();
                    const exprAst = expression();
                    let varName = 'x';
                    if (peek().type === TokenType.COMMA) {
                        advance();
                        const varTok = advance();
                        if (varTok.type === TokenType.IDENT) varName = varTok.value.toLowerCase();
                    }
                    expect(TokenType.RPAREN);
                    return { type: 'derivative', body: exprAst, varName };
                }

                // Special: nderiv2(expr) or diff2(expr) — second derivative
                if ((lowerName === 'deriv2' || lowerName === 'diff2' || lowerName === 'dderiv') &&
                    pos + 1 < tokens.length && tokens[pos + 1].type === TokenType.LPAREN) {
                    advance(); advance();
                    const exprAst = expression();
                    let varName = 'x';
                    if (peek().type === TokenType.COMMA) {
                        advance();
                        const varTok = advance();
                        if (varTok.type === TokenType.IDENT) varName = varTok.value.toLowerCase();
                    }
                    expect(TokenType.RPAREN);
                    return { type: 'derivative2', body: exprAst, varName };
                }

                // Special: integral(expr, a, b) or int(expr, a, b)
                if ((lowerName === 'integral' || lowerName === 'int' || lowerName === 'integrate') &&
                    pos + 1 < tokens.length && tokens[pos + 1].type === TokenType.LPAREN) {
                    advance(); advance();
                    const exprAst = expression();
                    expect(TokenType.COMMA);
                    const lower = expression();
                    expect(TokenType.COMMA);
                    const upper = expression();
                    let varName = 'x';
                    if (peek().type === TokenType.COMMA) {
                        advance();
                        const varTok = advance();
                        if (varTok.type === TokenType.IDENT) varName = varTok.value.toLowerCase();
                    }
                    expect(TokenType.RPAREN);
                    return { type: 'integral', body: exprAst, lower, upper, varName };
                }

                // Regular function call
                if (pos + 1 < tokens.length && tokens[pos + 1].type === TokenType.LPAREN) {
                    advance(); advance();
                    const args = [];
                    if (peek().type !== TokenType.RPAREN) {
                        args.push(expression());
                        while (match(TokenType.COMMA)) args.push(expression());
                    }
                    expect(TokenType.RPAREN);
                    return { type: 'call', name: lowerName, args };
                }
            }
            return atom();
        }

        function atom() {
            if (peek().type === TokenType.NUMBER) return { type: 'number', value: advance().value };
            if (peek().type === TokenType.IDENT) return { type: 'variable', name: advance().value.toLowerCase() };
            if (match(TokenType.LPAREN)) {
                const expr = expression();
                expect(TokenType.RPAREN);
                return expr;
            }
            if (match(TokenType.PIPE)) {
                const expr = expression();
                expect(TokenType.PIPE);
                return { type: 'call', name: 'abs', args: [expr] };
            }
            throw new Error(`Unexpected token: ${peek().type}`);
        }

        const ast = expression();
        return ast;
    }

    // ========================================================
    // EVALUATOR
    // ========================================================
    function evaluate(ast, vars = {}) {
        if (!ast) return NaN;

        switch (ast.type) {
            case 'number': return ast.value;
            case 'variable': {
                const name = ast.name;
                if (name in vars) return vars[name];
                if (name in CONSTANTS) return CONSTANTS[name];
                return NaN;
            }
            case 'binary': {
                const l = evaluate(ast.left, vars);
                const r = evaluate(ast.right, vars);
                switch (ast.op) {
                    case '+': return l + r;
                    case '-': return l - r;
                    case '*': return l * r;
                    case '/': return r === 0 ? NaN : l / r;
                    case '%': return r === 0 ? NaN : ((l % r) + r) % r;
                    case '^': return Math.pow(l, r);
                    default: return NaN;
                }
            }
            case 'unary': {
                const val = evaluate(ast.operand, vars);
                return ast.op === '-' ? -val : val;
            }
            case 'comparison': {
                const l = evaluate(ast.left, vars);
                const r = evaluate(ast.right, vars);
                switch (ast.op) {
                    case '<': return l < r ? 1 : 0;
                    case '>': return l > r ? 1 : 0;
                    case '<=': return l <= r ? 1 : 0;
                    case '>=': return l >= r ? 1 : 0;
                    case '==': return Math.abs(l - r) < 1e-10 ? 1 : 0;
                    case '!=': return Math.abs(l - r) >= 1e-10 ? 1 : 0;
                    default: return NaN;
                }
            }
            case 'call': {
                const name = ast.name;
                const argVals = ast.args.map(a => evaluate(a, vars));
                if (name in FUNCTIONS) return FUNCTIONS[name].fn(...argVals);
                return NaN;
            }
            case 'special': {
                // sum or prod
                const a = Math.floor(evaluate(ast.start, vars));
                const b = Math.floor(evaluate(ast.end, vars));
                if (ast.name === 'sum') return computeSum(ast.body, ast.varName, a, b, vars);
                if (ast.name === 'prod') return computeProduct(ast.body, ast.varName, a, b, vars);
                return NaN;
            }
            case 'derivative': {
                const varName = ast.varName;
                const point = vars[varName] || 0;
                return numericalDerivative(ast.body, point, vars, varName);
            }
            case 'derivative2': {
                const varName = ast.varName;
                const point = vars[varName] || 0;
                return numericalSecondDerivative(ast.body, point, vars, varName);
            }
            case 'integral': {
                const a = evaluate(ast.lower, vars);
                const b = evaluate(ast.upper, vars);
                return numericalIntegral(ast.body, a, b, vars, ast.varName);
            }
            default: return NaN;
        }
    }

    // ========================================================
    // COMPLEX EVALUATOR
    // ========================================================
    function evaluateComplex(ast, z) {
        // z is a Complex number mapped to the variable
        if (!ast) return new Complex(NaN);

        switch (ast.type) {
            case 'number': return new Complex(ast.value);
            case 'variable': {
                const name = ast.name;
                if (name === 'z' || name === 'x') return z;
                if (name === 'i') return new Complex(0, 1);
                if (name in CONSTANTS) return new Complex(CONSTANTS[name]);
                return new Complex(NaN);
            }
            case 'binary': {
                const l = evaluateComplex(ast.left, z);
                const r = evaluateComplex(ast.right, z);
                switch (ast.op) {
                    case '+': return l.add(r);
                    case '-': return l.sub(r);
                    case '*': return l.mul(r);
                    case '/': return l.div(r);
                    case '^': return l.pow(r);
                    default: return new Complex(NaN);
                }
            }
            case 'unary': {
                const val = evaluateComplex(ast.operand, z);
                return ast.op === '-' ? val.neg() : val;
            }
            case 'call': {
                const name = ast.name;
                const args = ast.args.map(a => evaluateComplex(a, z));
                // Complex-aware function dispatch
                switch (name) {
                    case 'sin': return Complex.sin(args[0]);
                    case 'cos': return Complex.cos(args[0]);
                    case 'tan': return Complex.tan(args[0]);
                    case 'exp': return Complex.exp(args[0]);
                    case 'ln':
                    case 'log': return Complex.log(args[0]);
                    case 'sqrt': return Complex.sqrt(args[0]);
                    case 'sinh': return Complex.sinh(args[0]);
                    case 'cosh': return Complex.cosh(args[0]);
                    case 'abs':
                    case 'cabs': return new Complex(args[0].abs());
                    case 're': return new Complex(args[0].re);
                    case 'im': return new Complex(args[0].im);
                    case 'conj': return args[0].conj();
                    case 'carg':
                    case 'arg': return new Complex(args[0].arg());
                    case 'gamma': return Complex.gamma(args[0]);
                    case 'zeta': return Complex.zeta(args[0]);
                    default: {
                        // Fallback: evaluate as real if all args are real-ish
                        if (name in FUNCTIONS) {
                            const realArgs = args.map(a => a.toReal());
                            const result = FUNCTIONS[name].fn(...realArgs);
                            return new Complex(result);
                        }
                        return new Complex(NaN);
                    }
                }
            }
            default: return new Complex(NaN);
        }
    }

    function normalizeLatexInput(input) {
        if (typeof input !== 'string') return '';
        let out = input;

        out = out
            .replace(/◻|□/g, '1')
            .replace(/\\left|\\right/g, '')
            .replace(/\\cdot|\\times/g, '*')
            .replace(/\\pi/g, 'pi')
            .replace(/\\theta/g, 'theta')
            .replace(/\\phi/g, 'phi')
            .replace(/\\tau/g, 'tau')
            .replace(/π/g, 'pi')
            .replace(/θ/g, 'theta')
            .replace(/ϕ|φ/g, 'phi')
            .replace(/τ/g, 'tau')
            .replace(/·/g, '*');

        let prev = '';
        while (out !== prev) {
            prev = out;
            out = out.replace(/\\frac\s*\{([^{}]+)\}\s*\{([^{}]+)\}/g, '(($1)/($2))');
            out = out.replace(/\\sqrt\s*\{([^{}]+)\}/g, 'sqrt($1)');
            out = out.replace(/\^\s*\{([^{}]+)\}/g, '^($1)');
        }

        out = out.replace(/\\([a-zA-Z]+)\s*\{/g, '$1(');
        out = out.replace(/[{}]/g, m => (m === '{' ? '(' : ')'));
        out = out.replace(/\\([a-zA-Z]+)/g, '$1');

        return out;
    }

    // ========================================================
    // EXPRESSION ANALYSIS
    // ========================================================
    function parseExpression(input) {
        const rawInput = String(input ?? '');
        input = normalizeLatexInput(rawInput);
        input = input.trim();
        if (!input) return null;

        // Detect polar: r = f(theta) or r = f(t)
        const polarMatch = input.match(/^r\s*=\s*(.+)$/i);
        if (polarMatch) {
            try {
                const tokens = tokenize(polarMatch[1]);
                const ast = parse(tokens);
                return {
                    type: 'polar',
                    ast,
                    raw: rawInput,
                    isDiscrete: false,
                    variables: findVariables(ast),
                };
            } catch (e) {
                return { type: 'error', error: e.message, raw: rawInput };
            }
        }

        // Detect parametric: (expr, expr) — a tuple for x(t), y(t)
        const paramMatch = input.match(/^\(\s*(.+?)\s*,\s*(.+?)\s*\)$/);
        if (paramMatch) {
            try {
                const tokensX = tokenize(paramMatch[1]);
                const tokensY = tokenize(paramMatch[2]);
                const astX = parse(tokensX);
                const astY = parse(tokensY);
                return {
                    type: 'parametric',
                    astX, astY,
                    raw: rawInput,
                    isDiscrete: false,
                    variables: new Set([...findVariables(astX), ...findVariables(astY)]),
                };
            } catch (e) {
                return { type: 'error', error: e.message, raw: rawInput };
            }
        }

        // Detect implicit: f(x,y) = g(x,y) (contains both x and y, with =)
        const implicitMatch = input.match(/^(.+?)\s*=\s*(.+)$/);
        if (implicitMatch) {
            const leftStr = implicitMatch[1];
            const rightStr = implicitMatch[2];

            try {
                const tokL = tokenize(leftStr);
                const astL = parse(tokL);
                const tokR = tokenize(rightStr);
                const astR = parse(tokR);

                const varsL = findVariables(astL);
                const varsR = findVariables(astR);
                const allVars = new Set([...varsL, ...varsR]);

                // If it has both x and y, it's implicit
                if (allVars.has('x') && allVars.has('y')) {
                    // f(x,y) - g(x,y) = 0
                    const ast = { type: 'binary', op: '-', left: astL, right: astR };
                    return {
                        type: 'implicit',
                        ast,
                        raw: rawInput,
                        isDiscrete: false,
                        variables: allVars,
                    };
                }

                // Single variable assignment: might be slider
                const varName = leftStr.trim().toLowerCase();
                if (!allVars.has('x') && !allVars.has('y') && !allVars.has('n') &&
                    !(varName in FUNCTIONS) && !(varName in CONSTANTS) &&
                    /^[a-z_]\w*$/.test(varName) && varName !== 'y' && varName !== 'r') {
                    // Check if right side uses variables
                    if (varsR.size === 0) {
                        const value = evaluate(astR, {});
                        return { type: 'slider', name: varName, value, ast: astR, raw: rawInput };
                    }
                    // Has variables — treat as function
                    return { type: 'function', ast: astR, raw: rawInput, isDiscrete: detectDiscrete(astR), variables: varsR };
                }

                // y = f(x) or standard assignment
                return { type: 'function', ast: astR, raw: rawInput, isDiscrete: detectDiscrete(astR), variables: varsR };
            } catch (e) {
                return { type: 'error', error: e.message, raw: rawInput };
            }
        }

        // Regular expression
        try {
            const tokens = tokenize(input);
            const ast = parse(tokens);
            const usedVars = findVariables(ast);

            // Check if it only uses z (complex mode hint)
            const isComplexHint = usedVars.has('z') || (usedVars.has('i') && !usedVars.has('x'));

            return {
                type: 'function',
                ast,
                raw: rawInput,
                isDiscrete: detectDiscrete(ast),
                variables: usedVars,
                isComplexHint,
            };
        } catch (e) {
            return { type: 'error', error: e.message, raw: rawInput };
        }
    }

    function findVariables(ast) {
        const vars = new Set();
        function walk(node) {
            if (!node) return;
            if (node.type === 'variable' && !(node.name in CONSTANTS) && node.name !== 'i') vars.add(node.name);
            if (node.left) walk(node.left);
            if (node.right) walk(node.right);
            if (node.operand) walk(node.operand);
            if (node.args) node.args.forEach(walk);
            if (node.body) walk(node.body);
            if (node.start) walk(node.start);
            if (node.end) walk(node.end);
            if (node.lower) walk(node.lower);
            if (node.upper) walk(node.upper);
            if (node.astX) walk(node.astX);
            if (node.astY) walk(node.astY);
        }
        walk(ast);
        return vars;
    }

    function findFunctions(ast) {
        const funcs = new Set();
        function walk(node) {
            if (!node) return;
            if (node.type === 'call') funcs.add(node.name);
            if (node.left) walk(node.left);
            if (node.right) walk(node.right);
            if (node.operand) walk(node.operand);
            if (node.args) node.args.forEach(walk);
            if (node.body) walk(node.body);
        }
        walk(ast);
        return funcs;
    }

    function detectDiscrete(ast) {
        const funcs = findFunctions(ast);
        for (const f of funcs) {
            if (f in FUNCTIONS && FUNCTIONS[f].discrete) return true;
        }
        return false;
    }

    // ========================================================
    // AUTOCOMPLETE
    // ========================================================
    function getAutocompleteSuggestions(prefix) {
        prefix = prefix.toLowerCase();
        const results = [];
        for (const [name, entry] of Object.entries(FUNCTIONS)) {
            if (name.startsWith(prefix)) {
                results.push({
                    name,
                    desc: entry.desc,
                    args: typeof entry.args === 'number' ? entry.args : entry.args[0],
                });
            }
        }
        // Special forms
        const specials = [
            { name: 'sum', desc: 'Summation Σ(expr,var,start,end)', args: 4 },
            { name: 'prod', desc: 'Product Π(expr,var,start,end)', args: 4 },
            { name: 'deriv', desc: 'Derivative d/dx f(x)', args: 1 },
            { name: 'diff', desc: 'Derivative d/dx f(x)', args: 1 },
            { name: 'integral', desc: 'Integral ∫ f(x) dx', args: 3 },
            { name: 'diff2', desc: 'Second derivative', args: 1 },
        ];
        for (const s of specials) {
            if (s.name.startsWith(prefix)) results.push(s);
        }
        for (const [name] of Object.entries(CONSTANTS)) {
            if (name.startsWith(prefix)) {
                results.push({ name, desc: `${CONSTANTS[name].toFixed(6)}`, args: 0 });
            }
        }
        return results.slice(0, 15);
    }

    // ========================================================
    // PUBLIC API
    // ========================================================
    return {
        parseExpression,
        evaluate,
        evaluateComplex,
        Complex,
        findVariables,
        findFunctions,
        detectDiscrete,
        normalizeLatexInput,
        getAutocompleteSuggestions,
        FUNCTIONS,
        CONSTANTS,
        tokenize,
        parse,
        numericalDerivative,
        numericalIntegral,
        computeSum,
    };
})();
