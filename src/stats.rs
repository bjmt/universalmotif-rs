//! Statistical helpers: binomial and hypergeometric probabilities, Fisher's
//! exact test and P-value adjustment.
//!
//! The binomial density uses Loader's saddle-point method, and hypergeometric
//! tails are summed as a series relative to the density at the cut-off, so
//! tails are computed directly rather than as one minus the other tail, and
//! stay accurate down to the smallest P-values.

const LN_2PI: f64 = 1.837_877_066_409_345_5;
const LN_SQRT_2PI: f64 = 0.918_938_533_204_672_8;

/// The error of Stirling's approximation, `ln(n!) - ln(sqrt(2 pi n) (n/e)^n)`.
fn stirlerr(n: f64) -> f64 {
    const S0: f64 = 1.0 / 12.0;
    const S1: f64 = 1.0 / 360.0;
    const S2: f64 = 1.0 / 1260.0;
    const S3: f64 = 1.0 / 1680.0;
    const S4: f64 = 1.0 / 1188.0;
    if n <= 15.0 {
        return libm::lgamma(n + 1.0) - (n + 0.5) * n.ln() + n - LN_SQRT_2PI;
    }
    let nn = n * n;
    if n > 500.0 {
        return (S0 - S1 / nn) / n;
    }
    if n > 80.0 {
        return (S0 - (S1 - S2 / nn) / nn) / n;
    }
    if n > 35.0 {
        return (S0 - (S1 - (S2 - S3 / nn) / nn) / nn) / n;
    }
    (S0 - (S1 - (S2 - (S3 - S4 / nn) / nn) / nn) / nn) / n
}

/// The deviance term `x ln(x / np) + np - x`, computed without cancellation
/// when `x` is close to `np`.
fn bd0(x: f64, np: f64) -> f64 {
    if (x - np).abs() < 0.1 * (x + np) {
        let mut v = (x - np) / (x + np);
        let mut s = (x - np) * v;
        if s.abs() < f64::MIN_POSITIVE {
            return s;
        }
        let mut ej = 2.0 * x * v;
        v *= v;
        for j in 1..1000 {
            ej *= v;
            let s1 = s + ej / (2 * j + 1) as f64;
            if s1 == s {
                return s1;
            }
            s = s1;
        }
    }
    x * (x / np).ln() + np - x
}

/// Binomial probability of `x` successes in `n` trials by Loader's saddle-point
/// method, with `q = 1 - p` given separately for accuracy.
fn dbinom_raw(x: f64, n: f64, p: f64, q: f64) -> f64 {
    if p == 0.0 {
        return if x == 0.0 { 1.0 } else { 0.0 };
    }
    if q == 0.0 {
        return if x == n { 1.0 } else { 0.0 };
    }
    if x == 0.0 {
        if n == 0.0 {
            return 1.0;
        }
        return if p > q {
            q.powf(n)
        } else {
            (n * (-p).ln_1p()).exp()
        };
    }
    if x == n {
        return if p > q {
            (n * (-q).ln_1p()).exp()
        } else {
            p.powf(n)
        };
    }
    if x < 0.0 || x > n {
        return 0.0;
    }
    let lc = stirlerr(n) - stirlerr(x) - stirlerr(n - x) - bd0(x, n * p) - bd0(n - x, n * q);
    let lf = LN_2PI + x.ln() + (-x / n).ln_1p();
    (lc - 0.5 * lf).exp()
}

/// Binomial probability of exactly `x` successes in `n` trials.
pub fn dbinom(x: u64, n: u64, p: f64) -> f64 {
    dbinom_raw(x as f64, n as f64, p, 1.0 - p)
}

/// Probability of at least `k` successes in `n` trials with success
/// probability `p`.
pub fn pbinom_upper(k: u64, n: u64, p: f64) -> f64 {
    if k == 0 {
        return 1.0;
    }
    if k > n || p <= 0.0 {
        return 0.0;
    }
    if p >= 1.0 {
        return 1.0;
    }
    let q = 1.0 - p;
    let (nf, ratio) = (n as f64, p / q);
    // Terms are summed until they no longer change the total.
    if k as f64 > nf * p {
        // The upper tail, where the terms only fall.
        let mut term = dbinom_raw(k as f64, nf, p, q);
        let mut sum = term;
        let mut i = k as f64;
        while i < nf && term > f64::EPSILON * 1e-3 * sum {
            term *= (nf - i) / (i + 1.0) * ratio;
            sum += term;
            i += 1.0;
        }
        sum.min(1.0)
    } else {
        // The upper tail is at least about a half, so it is one minus the
        // lower tail, summed downwards from k - 1.
        let mut term = dbinom_raw((k - 1) as f64, nf, p, q);
        let mut sum = term;
        let mut i = (k - 1) as f64;
        while i > 0.0 && term > f64::EPSILON * 1e-3 * sum {
            term *= i / (nf - i + 1.0) / ratio;
            sum += term;
            i -= 1.0;
        }
        (1.0 - sum).max(0.0)
    }
}

/// Probability of drawing `x` red balls in `n` draws without replacement
/// from `r` red and `b` black ones.
pub fn dhyper(x: u64, r: u64, b: u64, n: u64) -> f64 {
    if n < x || r < x || n - x > b {
        return 0.0;
    }
    if n == 0 {
        return if x == 0 { 1.0 } else { 0.0 };
    }
    let (x, r, b, n) = (x as f64, r as f64, b as f64, n as f64);
    let p = n / (r + b);
    let q = (r + b - n) / (r + b);
    let p1 = dbinom_raw(x, r, p, q);
    let p2 = dbinom_raw(n - x, b, p, q);
    let p3 = dbinom_raw(n, r + b, p, q);
    p1 * p2 / p3
}

/// The ratio of the lower tail to the density at `x`.
fn pdhyper(mut x: f64, nr: f64, nb: f64, n: f64) -> f64 {
    let mut sum = 0.0;
    let mut term = 1.0;
    while x > 0.0 && term >= f64::EPSILON * sum {
        term *= x * (nb - n + x) / (n + 1.0 - x) / (nr + 1.0 - x);
        sum += term;
        x -= 1.0;
    }
    1.0 + sum
}

/// Hypergeometric distribution function: the probability of at most `x` red
/// balls in `n` draws from `nr` red and `nb` black ones, or of more than `x` if
/// `lower` is false. Negative `x` is allowed.
pub fn phyper(x: i64, nr: u64, nb: u64, n: u64, lower: bool) -> f64 {
    let (mut x, mut nr, mut nb, n) = (x as f64, nr as f64, nb as f64, n as f64);
    let mut lower = lower;
    if n > nr + nb {
        return f64::NAN;
    }
    if x * (nr + nb) > n * nr {
        std::mem::swap(&mut nr, &mut nb);
        x = n - x - 1.0;
        lower = !lower;
    }
    let (zero, one) = if lower { (0.0, 1.0) } else { (1.0, 0.0) };
    if x < 0.0 || x < n - nb {
        return zero;
    }
    if x >= nr || x >= n {
        return one;
    }
    let d = dhyper(x as u64, nr as u64, nb as u64, n as u64);
    if d == 0.0 {
        return zero;
    }
    let pd = pdhyper(x, nr, nb, n);
    if lower { d * pd } else { 0.5 - d * pd + 0.5 }
}

/// One-sided Fisher's exact test of the 2x2 table `[[a, b], [c, d]]` for the
/// "greater" alternative: the probability of `a` or more in the top-left cell
/// given the margins.
pub fn fisher_greater(a: u64, b: u64, c: u64, d: u64) -> f64 {
    phyper(a as i64 - 1, a + c, b + d, a + b, false)
}

/// The natural log of [`fisher_greater`], which stays finite for P-values too
/// small for a double.
pub fn ln_fisher_greater(a: u64, b: u64, c: u64, d: u64) -> f64 {
    let p = fisher_greater(a, b, c, d);
    if p > 1e-300 {
        return p.ln();
    }
    // The upper tail from `a` to the largest possible count, summed in log
    // space from the hypergeometric densities.
    let (r, n, k) = (a + c, a + b + c + d, a + b);
    let lchoose = |n: u64, k: u64| {
        libm::lgamma(n as f64 + 1.0)
            - libm::lgamma(k as f64 + 1.0)
            - libm::lgamma((n - k) as f64 + 1.0)
    };
    let denom = lchoose(n, k);
    let mut log_sum = f64::NEG_INFINITY;
    for x in a..=r.min(k) {
        let term = lchoose(r, x) + lchoose(n - r, k - x) - denom;
        if term < log_sum - 40.0 {
            break;
        }
        log_sum = if log_sum > term {
            log_sum + (term - log_sum).exp().ln_1p()
        } else {
            term + (log_sum - term).exp().ln_1p()
        };
    }
    log_sum.min(0.0)
}

/// Benjamini-Hochberg q-values, in natural log, for P-values given as natural
/// logs, counting `m` tests in all; tests not given are taken to have larger
/// P-values. With `m` equal to the number given, these are the logs of
/// [`p_adjust_bh`].
pub fn ln_p_adjust_bh(ln_p: &[f64], m: f64) -> Vec<f64> {
    let n = ln_p.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| ln_p[a].total_cmp(&ln_p[b]));
    let mut q = vec![0.0; n];
    let mut prev = 0.0f64;
    for i in (0..n).rev() {
        let qv = (ln_p[order[i]] + (m / (i + 1) as f64).ln()).min(prev);
        prev = qv;
        q[order[i]] = qv;
    }
    q
}

/// Benjamini-Hochberg adjusted P-values.
pub fn p_adjust_bh(p: &[f64]) -> Vec<f64> {
    let n = p.len();
    let mut order: Vec<usize> = (0..n).collect();
    // Decreasing P-values; the sort is stable, so ties keep input order.
    order.sort_by(|&a, &b| p[b].total_cmp(&p[a]));
    let mut out = vec![0.0; n];
    let mut running = f64::INFINITY;
    for (rank, &i) in order.iter().enumerate() {
        let k = (n - rank) as f64;
        running = running.min(n as f64 / k * p[i]);
        out[i] = running.min(1.0);
    }
    out
}

/// Bonferroni adjusted P-values, capped at 1.
pub fn p_adjust_bonferroni(p: &[f64]) -> Vec<f64> {
    let n = p.len() as f64;
    p.iter().map(|&x| (x * n).min(1.0)).collect()
}

/// The median, or `None` for no values. For an even count it is the mean of the
/// middle two values.
pub fn median(values: &[f64]) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    let mut v = values.to_vec();
    v.sort_by(f64::total_cmp);
    let m = v.len() / 2;
    Some(if v.len() % 2 == 1 {
        v[m]
    } else {
        (v[m - 1] + v[m]) / 2.0
    })
}

#[cfg(test)]
mod tests {
    #[test]
    fn bh_matches_r() {
        // p.adjust(c(0.01, 0.04, 0.03, 0.005, 0.5), "BH")
        let got = super::p_adjust_bh(&[0.01, 0.04, 0.03, 0.005, 0.5]);
        let want = [0.025, 0.05, 0.05, 0.025, 0.5];
        for (g, w) in got.iter().zip(want) {
            assert!((g - w).abs() < 1e-15, "{got:?}");
        }
    }
}
