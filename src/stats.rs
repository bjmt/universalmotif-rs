//! Statistical helpers.

/// Benjamini-Hochberg adjusted P-values, as R's `p.adjust(method = "BH")`.
pub fn p_adjust_bh(p: &[f64]) -> Vec<f64> {
    let n = p.len();
    let mut order: Vec<usize> = (0..n).collect();
    // Decreasing P-values; R's order() is stable, so ties keep input order.
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
