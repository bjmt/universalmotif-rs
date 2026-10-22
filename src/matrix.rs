//! Dense motif matrices.

use crate::error::{Result, arg_err};

/// A dense `nrow x ncol` matrix of `f64`, where rows are alphabet letters and
/// columns are motif positions.
///
/// Values are stored position-major, with each column contiguous, which is
/// column-major order for an array of shape `(nrow, ncol)`.
#[derive(Debug, Clone, PartialEq)]
pub struct Matrix {
    nrow: usize,
    ncol: usize,
    data: Vec<f64>,
}

impl Matrix {
    /// Wrap position-major data, given column after column.
    pub fn from_vec(nrow: usize, ncol: usize, data: Vec<f64>) -> Result<Self> {
        if nrow == 0 {
            return arg_err("a motif matrix needs at least one row");
        }
        if data.len() != nrow * ncol {
            return arg_err(format!(
                "matrix data has {} values, expected {} ({} rows x {} columns)",
                data.len(),
                nrow * ncol,
                nrow,
                ncol
            ));
        }
        Ok(Matrix { nrow, ncol, data })
    }

    /// Build a matrix from a list of columns, which must all have the same
    /// length.
    pub fn from_columns<C: AsRef<[f64]>>(columns: &[C]) -> Result<Self> {
        let Some(first) = columns.first() else {
            return arg_err("a motif matrix needs at least one column");
        };
        let nrow = first.as_ref().len();
        let mut data = Vec::with_capacity(nrow * columns.len());
        for c in columns {
            let c = c.as_ref();
            if c.len() != nrow {
                return arg_err("all matrix columns must have the same length");
            }
            data.extend_from_slice(c);
        }
        Matrix::from_vec(nrow, columns.len(), data)
    }

    pub fn nrow(&self) -> usize {
        self.nrow
    }

    pub fn ncol(&self) -> usize {
        self.ncol
    }

    /// The raw position-major values.
    pub fn as_slice(&self) -> &[f64] {
        &self.data
    }

    pub fn into_vec(self) -> Vec<f64> {
        self.data
    }

    pub fn get(&self, row: usize, col: usize) -> f64 {
        self.data[col * self.nrow + row]
    }

    pub fn col(&self, j: usize) -> &[f64] {
        &self.data[j * self.nrow..(j + 1) * self.nrow]
    }

    pub fn col_mut(&mut self, j: usize) -> &mut [f64] {
        &mut self.data[j * self.nrow..(j + 1) * self.nrow]
    }

    pub fn columns(&self) -> impl ExactSizeIterator<Item = &[f64]> {
        self.data.chunks_exact(self.nrow)
    }

    pub fn col_sums(&self) -> Vec<f64> {
        self.columns().map(sum).collect()
    }

    /// Apply a column-wise transformation, producing a new matrix of the same
    /// shape.
    pub fn map_columns(&self, mut f: impl FnMut(&[f64]) -> Vec<f64>) -> Matrix {
        let mut data = Vec::with_capacity(self.data.len());
        for c in self.columns() {
            let out = f(c);
            debug_assert_eq!(out.len(), self.nrow);
            data.extend(out);
        }
        Matrix {
            nrow: self.nrow,
            ncol: self.ncol,
            data,
        }
    }

    /// A copy of the columns in `range`.
    pub fn slice_columns(&self, range: std::ops::Range<usize>) -> Matrix {
        Matrix {
            nrow: self.nrow,
            ncol: range.len(),
            data: self.data[range.start * self.nrow..range.end * self.nrow].to_vec(),
        }
    }

    /// Reverse the order of both rows and columns. For DNA and RNA this gives
    /// the reverse complement.
    pub fn reversed(&self) -> Matrix {
        let mut data = self.data.clone();
        data.reverse();
        Matrix {
            nrow: self.nrow,
            ncol: self.ncol,
            data,
        }
    }
}

/// Plain left-to-right sum, with no reordering or compensation.
pub(crate) fn sum(x: &[f64]) -> f64 {
    x.iter().fold(0.0, |a, &b| a + b)
}

/// Index of the first maximum. A value replaces the best so far only if it
/// is strictly greater, so NaN never replaces a value.
pub(crate) fn which_max(x: &[f64]) -> usize {
    let mut best = 0;
    for (i, &v) in x.iter().enumerate().skip(1) {
        if x[best] < v {
            best = i;
        }
    }
    best
}
