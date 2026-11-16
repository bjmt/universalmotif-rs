//! Progress reporting and cancellation.

use universalmotif::Error;
use universalmotif::create::{CreateOptions, from_strings};
use universalmotif::enrich::{EnrichOptions, enrich_motifs};
use universalmotif::finder::{FinderOptions, motif_finder};
use universalmotif::gapped::{GappedFinderOptions, gapped_motif_finder};
use universalmotif::progress::{Progress, with_progress};
use universalmotif::scan::{ScanOptions, scan_sequences};
use universalmotif::shuffle::{ShuffleOptions, create_sequences, shuffle_sequences};

fn seqs(n: usize, len: usize) -> Vec<String> {
    create_sequences(&['A', 'C', 'G', 'T'], len, n, None, 1, 1).unwrap()
}

#[test]
fn shuffling_counts_every_byte() {
    let s = seqs(20, 300);
    let p = Progress::new();
    with_progress(&p, || shuffle_sequences(&s, &ShuffleOptions::default())).unwrap();
    assert_eq!(p.stage(), "shuffling");
    assert_eq!(p.total(), 6000);
    assert_eq!(p.done(), 6000);
}

#[test]
fn scanning_counts_bases_times_motifs() {
    let s = seqs(10, 1000);
    let motifs = vec![
        from_strings(&["TGACTCA"], None, &CreateOptions::default()).unwrap(),
        from_strings(&["GATAAG"], None, &CreateOptions::default()).unwrap(),
    ];
    let p = Progress::new();
    let opts = ScanOptions {
        nthreads: 2,
        ..Default::default()
    };
    let hits = with_progress(&p, || scan_sequences(&motifs, &s, &opts)).unwrap();
    assert_eq!(p.stage(), "scanning");
    assert_eq!(p.total(), 20_000);
    assert_eq!(p.done(), p.total());
    // Reporting changes nothing about the result.
    assert_eq!(hits, scan_sequences(&motifs, &s, &opts).unwrap());
}

#[test]
fn the_finders_count_their_searches() {
    let s = seqs(60, 120);
    let p = Progress::new();
    let opts = FinderOptions {
        min_width: 6,
        max_width: 8,
        nmotifs: 2,
        ..Default::default()
    };
    with_progress(&p, || motif_finder(&s, None, &opts)).unwrap();
    assert_eq!(p.stage(), "searching");
    assert_eq!((p.done(), p.total()), (6, 6));

    let p = Progress::new();
    let opts = GappedFinderOptions {
        k: vec![4, 5],
        nmotifs: 2,
        ..Default::default()
    };
    with_progress(&p, || gapped_motif_finder(&s, None, &opts)).unwrap();
    // The scans the search makes of each candidate do not take over.
    assert_eq!(p.stage(), "searching");
    assert_eq!((p.done(), p.total()), (4, 4));
}

#[test]
fn enrichment_reports_its_shuffle_and_its_scan() {
    let s = seqs(30, 200);
    let m = from_strings(&["TGACTCA"], None, &CreateOptions::default()).unwrap();
    let p = Progress::new();
    with_progress(&p, || {
        enrich_motifs(&[m], &s, None, &EnrichOptions::default())
    })
    .unwrap();
    assert_eq!(p.stage(), "scanning");
    assert_eq!(p.done(), p.total());
}

#[test]
fn a_cancelled_handle_stops_the_work() {
    let s = seqs(20, 500);
    let m = from_strings(&["TGACTCA"], None, &CreateOptions::default()).unwrap();
    let p = Progress::new();
    p.cancel();
    let err = |r: Result<(), Error>| r.unwrap_err();
    assert_eq!(
        err(with_progress(&p, || shuffle_sequences(&s, &ShuffleOptions::default())).map(|_| ())),
        Error::Cancelled
    );
    assert_eq!(
        err(with_progress(&p, || scan_sequences(
            std::slice::from_ref(&m),
            &s,
            &ScanOptions::default()
        ))
        .map(|_| ())),
        Error::Cancelled
    );
    assert_eq!(
        err(with_progress(&p, || motif_finder(&s, None, &FinderOptions::default())).map(|_| ())),
        Error::Cancelled
    );
    assert_eq!(
        err(with_progress(&p, || gapped_motif_finder(
            &s,
            None,
            &GappedFinderOptions::default()
        ))
        .map(|_| ())),
        Error::Cancelled
    );
    assert_eq!(Error::Cancelled.to_string(), "cancelled");
}

#[test]
fn handles_nest_and_are_put_back() {
    let s = seqs(5, 100);
    let (outer, inner) = (Progress::new(), Progress::new());
    with_progress(&outer, || {
        with_progress(&inner, || shuffle_sequences(&s, &ShuffleOptions::default())).unwrap();
        assert_eq!(outer.total(), 0);
        shuffle_sequences(&s, &ShuffleOptions::default()).unwrap();
    });
    assert_eq!((inner.done(), outer.done()), (500, 500));
    // Outside any handle, nothing is counted.
    shuffle_sequences(&s, &ShuffleOptions::default()).unwrap();
    assert_eq!(outer.done(), 500);
}
