//! Choosing motifs by the targets they add, and refining them from their hits.

use universalmotif::Motif;
use universalmotif::create::{CreateOptions, from_strings};
use universalmotif::select::{RefineOptions, SelectOptions, Site, refine_motifs, select_motifs};
use universalmotif::shuffle::create_sequences;

fn motif(site: &str, name: &str) -> Motif {
    let opts = CreateOptions {
        name: name.into(),
        pseudocount: 1.0,
        nsites: Some(20.0),
        ..Default::default()
    };
    from_strings(&[site], None, &opts).unwrap()
}

/// Random sequences with `site` planted at position `at` of every `every`th
/// one, starting from `offset`.
fn plant(seqs: &mut [String], site: &str, every: usize, offset: usize, at: usize) -> Vec<usize> {
    let mut planted = Vec::new();
    for (i, s) in seqs.iter_mut().enumerate().skip(offset).step_by(every) {
        s.replace_range(at..at + site.len(), site);
        planted.push(i);
    }
    planted
}

fn sequences() -> Vec<String> {
    let mut seqs = create_sequences(&['A', 'C', 'G', 'T'], 300, 250, None, 3, 1).unwrap();
    plant(&mut seqs, "TTGACTCATG", 2, 0, 100);
    plant(&mut seqs, "CAGATAAGGC", 5, 1, 180);
    seqs
}

#[test]
fn motifs_are_ranked_by_the_sequences_they_add() {
    let seqs = sequences();
    let motifs = vec![
        motif("AGATAAG", "gata"),
        motif("TGACTCA", "ap1"),
        motif("TGACTCAT", "ap1_long"),
        motif("CCCGGG", "rare"),
    ];
    let sel = select_motifs(&motifs, &seqs, None, None, &SelectOptions::default()).unwrap();
    assert_eq!(sel.n_targets, 300);
    // AP-1 is in half the sequences, GATA in a fifth.
    assert_eq!(sel.rank[1], Some(1));
    assert_eq!(sel.rank[0], Some(2));
    assert!(sel.targets[1] >= 150 && sel.targets[0] >= 60);
    // The longer AP-1 reaches the same sequences, so it adds little.
    assert!(sel.gain[2] < 0.01 && !sel.kept[2]);
    assert!(sel.kept[0] && sel.kept[1]);
    let total: f64 = sel.gain.iter().sum();
    let last = sel.cumulative.iter().cloned().fold(0.0, f64::max);
    assert!((total - last).abs() < 1e-12);
    // Gains never grow along the order.
    let mut by_rank: Vec<(usize, f64)> = sel
        .rank
        .iter()
        .zip(&sel.gain)
        .filter_map(|(r, g)| r.map(|r| (r, *g)))
        .collect();
    by_rank.sort_by_key(|x| x.0);
    assert!(by_rank.windows(2).all(|w| w[0].1 >= w[1].1));
}

#[test]
fn sites_are_targets_when_given() {
    let seqs = sequences();
    let sites: Vec<Site> = (0..seqs.len())
        .step_by(5)
        .map(|i| Site {
            sequence: i + 1,
            start: 178,
            end: 192,
        })
        .filter(|s| s.sequence < seqs.len())
        .collect();
    let motifs = vec![motif("TGACTCA", "ap1"), motif("AGATAAG", "gata")];
    let opts = SelectOptions {
        min_overlap: 5,
        ..Default::default()
    };
    let sel = select_motifs(&motifs, &seqs, Some(&sites), None, &opts).unwrap();
    assert_eq!(sel.n_targets, sites.len());
    // The sites hold the GATA sites, which AP-1 hits only by chance.
    assert_eq!(sel.rank[1], Some(1));
    assert!(sel.targets[1] as f64 >= 0.95 * sites.len() as f64);
    assert!(sel.gain[0] < 0.1);
    let bad = [Site {
        sequence: 999,
        start: 0,
        end: 5,
    }];
    assert!(select_motifs(&motifs, &seqs, Some(&bad), None, &opts).is_err());
}

#[test]
fn similar_motifs_are_left_out() {
    let seqs = sequences();
    let motifs = vec![
        motif("TGACTCA", "ap1"),
        motif("TGAGTCA", "ap1_other"),
        motif("AGATAAG", "gata"),
    ];
    // Say the two AP-1 spellings are near copies, as their CWMs might be.
    let sim = vec![
        vec![1.0, 0.95, 0.1],
        vec![0.95, 1.0, 0.1],
        vec![0.1, 0.1, 1.0],
    ];
    let opts = SelectOptions {
        max_similarity: Some(0.9),
        min_gain: 0.0,
        ..Default::default()
    };
    let sel = select_motifs(&motifs, &seqs, None, Some(&sim), &opts).unwrap();
    assert_eq!(sel.similar_to[1], Some(0));
    assert_eq!(sel.rank[1], None);
    assert!(!sel.kept[1] && sel.kept[0] && sel.kept[2]);
    assert!(select_motifs(&motifs, &seqs, None, Some(&sim[..2]), &opts).is_err());
}

#[test]
fn refinement_finds_the_whole_site() {
    let mut seqs = create_sequences(&['A', 'C', 'G', 'T'], 400, 250, None, 5, 1).unwrap();
    plant(&mut seqs, "GCTTGACTCATGAC", 2, 0, 100);
    // Start from the middle of the site.
    let start = motif("TGACTCA", "ap1");
    let out = refine_motifs(
        std::slice::from_ref(&start),
        &seqs,
        None,
        &RefineOptions::default(),
    )
    .unwrap();
    let r = &out[0];
    assert!(r.sites >= 190, "{} sites", r.sites);
    assert_eq!(r.motif.name, "ap1");
    assert_eq!(r.motif.motif_type, start.motif_type);
    let c = r.motif.consensus.clone().unwrap();
    // Read on either strand, the refined motif spans the planted flanks.
    let rc: String = c
        .chars()
        .rev()
        .map(|x| match x {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            o => o,
        })
        .collect();
    assert!(
        c.contains("GCTTGACTCATGAC") || rc.contains("GCTTGACTCATGAC"),
        "{c}"
    );
    assert!(r.motif.width() <= 14 + 2);
}

#[test]
fn refinement_reads_letters_outside_the_alphabet_as_ambiguity_codes() {
    // Sites with an N or an R, on either strand, once gave columns of
    // different sums, which no motif accepts.
    let mut seqs = create_sequences(&['A', 'C', 'G', 'T'], 400, 250, None, 5, 1).unwrap();
    plant(&mut seqs, "GCTTGACTCATGAC", 2, 0, 100);
    plant(&mut seqs, "GCTTGACTCANGAC", 8, 0, 100);
    plant(&mut seqs, "GTCRTGAGTCAAGC", 8, 4, 100);
    let start = motif("TGACTCA", "ap1");
    let out = refine_motifs(&[start], &seqs, None, &RefineOptions::default()).unwrap();
    let sums = out[0].motif.matrix.col_sums();
    assert!(sums.iter().all(|s| (s - sums[0]).abs() < 1e-9), "{sums:?}");
}

#[test]
fn refinement_respects_its_limits() {
    let mut seqs = create_sequences(&['A', 'C', 'G', 'T'], 200, 250, None, 7, 1).unwrap();
    plant(&mut seqs, "GCTTGACTCATGAC", 2, 0, 100);
    let start = motif("TGACTCA", "ap1");
    let rare = motif("CGCGCGCG", "rare");
    let opts = RefineOptions {
        max_width: 9,
        ..Default::default()
    };
    let out = refine_motifs(&[start, rare.clone()], &seqs, None, &opts).unwrap();
    assert!(out[0].motif.width() <= 9);
    // Too few hits to rebuild from.
    assert_eq!(out[1].sites, 0);
    assert_eq!(out[1].motif, rare);
    // Only hits in the sites count.
    let sites = [Site {
        sequence: 0,
        start: 95,
        end: 120,
    }];
    let few = refine_motifs(
        &[motif("TGACTCA", "ap1")],
        &seqs,
        Some(&sites),
        &RefineOptions::default(),
    )
    .unwrap();
    assert_eq!(few[0].sites, 0);
    assert!(
        refine_motifs(
            &[rare],
            &seqs,
            None,
            &RefineOptions {
                max_width: 0,
                ..Default::default()
            }
        )
        .is_err()
    );
}

#[test]
fn competitive_refinement_keeps_a_weak_motif_apart() {
    let mut seqs = create_sequences(&['A', 'C', 'G', 'T'], 400, 250, None, 9, 1).unwrap();
    plant(&mut seqs, "TTGACTCATG", 2, 0, 100);
    plant(&mut seqs, "AGATAAGA", 4, 1, 180);
    // A loose motif that matches some AP-1 sites as well as its own.
    let weak = motif("AGATAAG", "gata");
    let ap1 = motif("TGACTCA", "ap1");
    let alone = refine_motifs(
        &[ap1.clone(), weak.clone()],
        &seqs,
        None,
        &RefineOptions {
            competitive: false,
            ..Default::default()
        },
    )
    .unwrap();
    let shared = refine_motifs(&[ap1, weak], &seqs, None, &RefineOptions::default()).unwrap();
    let c = shared[1].motif.consensus.clone().unwrap();
    assert!(c.contains("GATAAG") || c.contains("CTTATC"), "{c}");
    assert!(alone[1].sites > 0);
}

#[test]
fn gapped_motifs_are_refined_part_by_part() {
    let mut seqs = create_sequences(&['A', 'C', 'G', 'T'], 300, 250, None, 11, 1).unwrap();
    // Two parts with flanks, a spacer of 3 to 6 random letters between them.
    let spacers = ["ACG", "TTCA", "GACTT", "CATGCA"];
    for (i, s) in seqs.iter_mut().enumerate().step_by(2) {
        let site = format!("GCTGACTCA{}AGATAAGCG", spacers[(i / 2) % 4]);
        s.replace_range(100..100 + site.len(), &site);
    }
    let start =
        universalmotif::ops::add_gap(&motif("TGACTCAAGATAAG", "pair"), Some(&[7]), &[3], &[6])
            .unwrap();
    let out = refine_motifs(
        std::slice::from_ref(&start),
        &seqs,
        None,
        &RefineOptions::default(),
    )
    .unwrap();
    let r = &out[0].motif;
    assert!(out[0].sites >= 140, "{} sites", out[0].sites);
    assert!(r.gap.isgapped);
    assert_eq!(
        (r.gap.mingap.clone(), r.gap.maxgap.clone()),
        (vec![3], vec![6])
    );
    // Each part grew to take in its planted flank.
    assert_eq!(r.width(), 18, "{:?}", r.consensus);
    assert_eq!(r.gap.gaploc, vec![9]);
    let c = r.consensus.clone().unwrap();
    assert!(
        c.starts_with("GCTGACTCA") && c.ends_with("AGATAAGCG"),
        "{c}"
    );
}

#[test]
fn stability_counts_reruns_that_find_a_motif_again() {
    use universalmotif::finder::{FinderOptions, motif_finder};
    use universalmotif::select::{StabilityOptions, motif_stability};

    let mut seqs = create_sequences(&['A', 'C', 'G', 'T'], 300, 100, None, 21, 1).unwrap();
    plant(&mut seqs, "TTGACTCATG", 2, 0, 40);
    let finder = FinderOptions {
        nmotifs: 3,
        seed: 1,
        ..Default::default()
    };
    let found = motif_finder(&seqs, None, &finder).unwrap();
    assert!(!found.is_empty());
    let mut motifs: Vec<Motif> = found.into_iter().map(|f| f.motif).collect();
    motifs.push(motif("CGCGTACGCG", "absent"));
    let opts = StabilityOptions::default();
    let st = motif_stability(&motifs, &seqs, None, &finder, &opts).unwrap();
    assert_eq!(st.len(), motifs.len());
    assert_eq!(st[0], 1.0, "{st:?}");
    assert_eq!(*st.last().unwrap(), 0.0);
    assert_eq!(
        st,
        motif_stability(&motifs, &seqs, None, &finder, &opts).unwrap()
    );
    assert!(
        motif_stability(
            &motifs,
            &seqs,
            None,
            &finder,
            &StabilityOptions {
                rounds: 0,
                ..opts.clone()
            }
        )
        .is_err()
    );
    assert!(
        motif_stability(
            &motifs,
            &seqs,
            None,
            &finder,
            &StabilityOptions {
                fraction: 1.0,
                ..opts
            }
        )
        .is_err()
    );
}
