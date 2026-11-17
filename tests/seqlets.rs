//! Calling seqlets in contribution scores.

use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};
use universalmotif::seqlets::{SeqletOptions, SeqletSign, call_seqlets};

/// Noise around zero, with a block of `value` over `start..end` in every
/// `every`th sequence.
fn tracks(n: usize, len: usize, blocks: &[(usize, usize, f64)], every: usize) -> Vec<Vec<f64>> {
    let mut rng = StdRng::seed_from_u64(4);
    (0..n)
        .map(|i| {
            let mut x: Vec<f64> = (0..len).map(|_| rng.random_range(-0.02..0.02)).collect();
            if i % every == 0 {
                for &(a, b, v) in blocks {
                    x[a..b].iter_mut().for_each(|p| *p += v);
                }
            }
            x
        })
        .collect()
}

#[test]
fn blocks_of_high_scores_are_called() {
    let x = tracks(200, 300, &[(100, 110, 0.3)], 2);
    let opts = SeqletOptions {
        flank: 0,
        ..Default::default()
    };
    let s = call_seqlets(&x, None, &opts).unwrap();
    let pos: Vec<_> = s.iter().filter(|s| s.positive).collect();
    assert!(pos.len() >= 95, "{} seqlets", pos.len());
    // Seqlets may take in a few positions of the noise around a block.
    let on_block = pos
        .iter()
        .filter(|s| s.sequence % 2 == 0 && s.core_start >= 95 && s.core_end <= 115)
        .count();
    assert!(
        on_block >= 95 && on_block as f64 >= 0.9 * pos.len() as f64,
        "{on_block} of {}",
        pos.len()
    );
    assert!(pos.iter().all(|s| s.score > 0.0 && s.pvalue < 0.01));
    // Seqlets come in order of sequence and position.
    assert!(
        s.windows(2)
            .all(|w| (w[0].sequence, w[0].core_start) <= (w[1].sequence, w[1].core_start))
    );
}

#[test]
fn noise_alone_gives_almost_nothing_with_the_shuffle_filter() {
    let x = tracks(200, 300, &[], 1);
    let filtered = call_seqlets(&x, None, &SeqletOptions::default()).unwrap();
    let raw = call_seqlets(
        &x,
        None,
        &SeqletOptions {
            fdr: None,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(
        raw.len() > 10 * filtered.len().max(1),
        "{} raw, {} kept",
        raw.len(),
        filtered.len()
    );
}

#[test]
fn a_quiet_stretch_splits_two_blocks() {
    // Two blocks with a run of zeros between them: no stretch spanning both
    // is a seqlet, since the run inside it is unremarkable.
    let x = tracks(100, 200, &[(50, 58, 0.3), (68, 76, 0.3)], 1);
    let opts = SeqletOptions {
        flank: 0,
        fdr: None,
        ..Default::default()
    };
    let s = call_seqlets(&x, None, &opts).unwrap();
    let spanning = s
        .iter()
        .filter(|s| s.positive && s.core_start < 58 && s.core_end > 68)
        .count();
    assert_eq!(spanning, 0);
    let first = s
        .iter()
        .filter(|s| s.positive && s.core_start >= 48 && s.core_end <= 60)
        .count();
    assert!(first >= 95);
}

#[test]
fn signs_and_flanks() {
    let x = tracks(100, 200, &[(20, 30, -0.3), (150, 160, 0.3)], 1);
    let opts = SeqletOptions {
        flank: 3,
        ..Default::default()
    };
    let s = call_seqlets(&x, None, &opts).unwrap();
    let neg: Vec<_> = s.iter().filter(|s| !s.positive).collect();
    assert!(
        neg.len() >= 95
            && neg
                .iter()
                .all(|s| s.score < 0.0 && s.core_start >= 18 && s.core_end <= 32)
    );
    assert!(
        s.iter()
            .all(|s| s.start + 3 == s.core_start || s.start == 0)
    );
    assert!(s.iter().all(|s| s.end == s.core_end + 3 || s.end == 200));
    let pos_only = call_seqlets(
        &x,
        None,
        &SeqletOptions {
            sign: SeqletSign::Positive,
            ..opts.clone()
        },
    )
    .unwrap();
    assert!(pos_only.iter().all(|s| s.positive));
    assert_eq!(pos_only.len(), s.iter().filter(|s| s.positive).count());
}

#[test]
fn bad_options_and_empty_tracks() {
    let x = tracks(5, 50, &[], 1);
    for opts in [
        SeqletOptions {
            threshold: 0.0,
            ..Default::default()
        },
        SeqletOptions {
            min_len: 10,
            max_len: 5,
            ..Default::default()
        },
        SeqletOptions {
            fdr: Some(1.5),
            ..Default::default()
        },
    ] {
        assert!(call_seqlets(&x, None, &opts).is_err());
    }
    assert!(
        call_seqlets(&[vec![0.0; 40]], None, &SeqletOptions::default())
            .unwrap()
            .is_empty()
    );
    assert!(
        call_seqlets(&[vec![1.0, f64::NAN, 2.0]], None, &SeqletOptions::default())
            .unwrap()
            .is_empty()
    );
}

#[test]
fn null_scores_set_the_null() {
    let x = tracks(200, 300, &[(100, 110, 0.3)], 2);
    let opts = SeqletOptions {
        flank: 0,
        sign: SeqletSign::Positive,
        ..Default::default()
    };
    // Null scores with blocks as strong as the real ones make every seqlet
    // look false; quiet null scores, twice as many, keep the real ones.
    let loud = tracks(200, 300, &[(50, 60, 0.3)], 1);
    let quiet = tracks(400, 300, &[], 1);
    let none = call_seqlets(&x, Some(&loud), &opts).unwrap();
    let kept = call_seqlets(&x, Some(&quiet), &opts).unwrap();
    assert!(none.len() < 10, "{}", none.len());
    assert!(kept.len() >= 95, "{}", kept.len());
}

#[test]
fn single_precision_scores_give_the_same_seqlets() {
    // Scores that f32 holds exactly, so that both precisions see the same values.
    let x: Vec<Vec<f64>> = tracks(100, 200, &[(80, 90, 0.25)], 2)
        .into_iter()
        .map(|s| s.into_iter().map(|v| v as f32 as f64).collect())
        .collect();
    let x32: Vec<Vec<f32>> = x
        .iter()
        .map(|s| s.iter().map(|&v| v as f32).collect())
        .collect();
    let rows: Vec<&[f32]> = x32.iter().map(Vec::as_slice).collect();
    let opts = SeqletOptions::default();
    let a = call_seqlets(&x, None, &opts).unwrap();
    let b = call_seqlets(&rows, None, &opts).unwrap();
    assert!(!a.is_empty());
    assert_eq!(
        a.iter()
            .map(|s| (s.sequence, s.start, s.end))
            .collect::<Vec<_>>(),
        b.iter()
            .map(|s| (s.sequence, s.start, s.end))
            .collect::<Vec<_>>()
    );
}

#[test]
fn one_null_makes_batches_agree() {
    let x = tracks(200, 300, &[(100, 110, 0.3)], 2);
    let opts = SeqletOptions {
        fdr: None,
        ..Default::default()
    };
    // The scores' own shuffles have the null of the scores themselves.
    assert_eq!(
        call_seqlets(&x, Some(&x), &opts).unwrap(),
        call_seqlets(&x, None, &opts).unwrap()
    );
    // With the same null scores for every batch, calling in batches gives the
    // seqlets of one call.
    let null = tracks(60, 300, &[], 1);
    let whole = call_seqlets(&x, Some(&null), &opts).unwrap();
    let mut batched = Vec::new();
    for (k, batch) in x.chunks(50).enumerate() {
        for mut s in call_seqlets(batch, Some(&null), &opts).unwrap() {
            s.sequence += 50 * k;
            batched.push(s);
        }
    }
    assert_eq!(batched, whole);
    // Without them, each batch has a null of its own and the P-values differ.
    let own = call_seqlets(&x[50..100], None, &opts).unwrap();
    let given = call_seqlets(&x[50..100], Some(&null), &opts).unwrap();
    assert!(own.iter().zip(&given).any(|(a, b)| a.pvalue != b.pvalue));
}
