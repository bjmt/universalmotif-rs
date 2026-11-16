//! Progress of long-running functions, and a way to stop them.
//!
//! Functions that can take a while, such as scanning a genome or searching
//! for motifs, report their progress to a [`Progress`] handle when they are
//! called inside [`with_progress`]. They name each stage of the work, say
//! how much of it there is, and count it as it is done, so that another
//! thread can read the counts and show them. The crate itself never prints
//! anything. Calling [`Progress::cancel`] makes the function stop soon
//! after, returning [`Error::Cancelled`].
//!
//! The functions that report are
//! [`scan_table`](crate::scan::scan_table) and the others built on the same
//! scanner (which count bases times motifs), [`motif_finder`] and
//! [`gapped_motif_finder`] (motifs searched for), and
//! [`shuffle_sequences`](crate::shuffle::shuffle_sequences) (bytes).
//! Functions that call these, such as
//! [`enrich_motifs`](crate::enrich::enrich_motifs), report each call as a
//! stage of its own.
//!
//! [`motif_finder`]: crate::finder::motif_finder
//! [`gapped_motif_finder`]: crate::gapped::gapped_motif_finder
//!
//! ```
//! use universalmotif::progress::{Progress, with_progress};
//! use universalmotif::shuffle::{ShuffleOptions, shuffle_sequences};
//!
//! let progress = Progress::new();
//! let seqs = vec!["ACGTACGTAAACCCGGGTTT".to_string(); 10];
//! with_progress(&progress, || shuffle_sequences(&seqs, &ShuffleOptions::default())).unwrap();
//! assert_eq!(progress.stage(), "shuffling");
//! assert_eq!(progress.done(), progress.total());
//! ```

use std::cell::RefCell;
use std::fmt;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex};

use crate::error::{Error, Result};

#[derive(Default)]
struct Inner {
    done: AtomicU64,
    total: AtomicU64,
    cancelled: AtomicBool,
    stage: Mutex<String>,
}

/// A handle on the progress of a running function, shared by cloning.
///
/// Clones refer to the same counts, so the function can report to one
/// while another thread reads them from the other.
#[derive(Clone, Default)]
pub struct Progress(Arc<Inner>);

impl Progress {
    /// A new handle, with nothing done yet.
    pub fn new() -> Self {
        Self::default()
    }

    /// The units of work done in the current stage.
    pub fn done(&self) -> u64 {
        self.0.done.load(Ordering::Relaxed)
    }

    /// The units of work in the current stage, or 0 before one starts.
    pub fn total(&self) -> u64 {
        self.0.total.load(Ordering::Relaxed)
    }

    /// The name of the current stage, such as `"scanning"`. When the total
    /// is not known ahead, the name ends with what is counted, such as
    /// `"scanning bases"`.
    pub fn stage(&self) -> String {
        self.0.stage.lock().map(|s| s.clone()).unwrap_or_default()
    }

    /// Asks the function to stop. It returns [`Error::Cancelled`] soon
    /// after, at the next point where it counts its progress.
    pub fn cancel(&self) {
        self.0.cancelled.store(true, Ordering::Relaxed);
    }

    /// Whether [`cancel`](Self::cancel) has been called.
    pub fn is_cancelled(&self) -> bool {
        self.0.cancelled.load(Ordering::Relaxed)
    }

    fn begin(&self, stage: &str, total: u64) {
        if let Ok(mut s) = self.0.stage.lock() {
            *s = stage.to_string();
        }
        self.0.done.store(0, Ordering::Relaxed);
        self.0.total.store(total, Ordering::Relaxed);
    }
}

impl fmt::Debug for Progress {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Progress")
            .field("stage", &self.stage())
            .field("done", &self.done())
            .field("total", &self.total())
            .field("cancelled", &self.is_cancelled())
            .finish()
    }
}

thread_local! {
    static CURRENT: RefCell<Option<Progress>> = const { RefCell::new(None) };
}

/// Runs `f` with `progress` as the handle that the functions it calls on
/// this thread report to, restoring the previous one afterwards.
pub fn with_progress<T>(progress: &Progress, f: impl FnOnce() -> T) -> T {
    struct Restore(Option<Progress>);
    impl Drop for Restore {
        fn drop(&mut self) {
            CURRENT.with(|c| *c.borrow_mut() = self.0.take());
        }
    }
    let _restore = Restore(CURRENT.with(|c| c.borrow_mut().replace(progress.clone())));
    f()
}

/// Runs `f` with no handle, so that functions used internally as parts of
/// a larger one do not report stages of their own.
pub(crate) fn without_progress<T>(f: impl FnOnce() -> T) -> T {
    let previous = CURRENT.with(|c| c.borrow_mut().take());
    let out = f();
    CURRENT.with(|c| *c.borrow_mut() = previous);
    out
}

/// The handle set by [`with_progress`] on this thread, if any. Functions
/// read it as they start, before handing work to other threads.
pub(crate) fn current() -> Option<Progress> {
    CURRENT.with(|c| c.borrow().clone())
}

/// Starts a stage of `total` units, or of an unknown number with a total of
/// 0, named as in [`Progress::stage`].
pub(crate) fn begin(p: &Option<Progress>, stage: &str, total: u64) {
    if let Some(p) = p {
        p.begin(stage, total);
    }
}

/// Counts `n` units as done.
pub(crate) fn add(p: &Option<Progress>, n: u64) {
    if let Some(p) = p {
        p.0.done.fetch_add(n, Ordering::Relaxed);
    }
}

/// Whether the caller has asked the function to stop.
pub(crate) fn cancelled(p: &Option<Progress>) -> bool {
    p.as_ref().is_some_and(Progress::is_cancelled)
}

/// [`Error::Cancelled`] if the caller has asked the function to stop.
pub(crate) fn check(p: &Option<Progress>) -> Result<()> {
    if cancelled(p) {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
