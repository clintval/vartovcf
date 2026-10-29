//! Folds the partial rows VarDict writes for a multi-nucleotide variant (MNV) back into it.
//!
//! VarDict grows an MNV through a read only while the read goes on with well-called bases, so a read
//! that ends inside the MNV, has a no-call (N) there or a base below VarDict's `-q` plus 5 adds to a
//! shorter allele instead: an SNV or a shorter MNV carrying the MNV's ALT bases over its own span.
//! VarDict folds some of these back itself, the partials ending at the MNV's last base and the ones
//! at its first base whose reads end, on average, within the partial, and writes the rest as rows
//! of their own.
//!
//! A partial row is folded when its reads are plausible leakage from the MNV. Each read of an MNV
//! molecule lands in any one partial row with a small probability, the leak rate, which depends on
//! the library: read length, fragment size, base-quality and no-call rates. A partial with `c` reads
//! under an MNV with `n` is then a draw from a binomial over the `n + c` reads, and it is folded
//! when `c < n` and `P(X >= c) >= LEAK_TAIL` at [`LeakRate`]. A row with at least as many reads as
//! the MNV is a separate allele, or the MNV's reads counted again, at any rate.
//!
//! With the reads VarDict counted, the partial's own reads decide instead ([`Rule::Reads`]): each
//! read's template is followed to the MNV's other changed bases, as [`crate::reads`] does, and the
//! reads that are the MNV's molecules, or consistent with leakage, move to the MNV while the reads
//! of another allele stay on the partial, which is then written with only them.

use std::collections::HashMap;
use std::path::PathBuf;

use crate::reads::ReadOptions;
use crate::record::TumorOnlyVariant;
use crate::vcf::Numeric;

/// The least `P(X >= c)` a partial's `c` reads may have under the leak binomial and still be folded.
pub const LEAK_TAIL: f64 = 0.05;

/// The probability that a read of an MNV molecule is counted toward any one partial row.
///
/// The rate varies with the library and the MNV's sequence context; the README describes how to
/// estimate it from a run without `--fix-leaky-mnv`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LeakRate(f64);

impl Default for LeakRate {
    fn default() -> Self {
        LeakRate(0.05)
    }
}

impl LeakRate {
    /// A leak rate, which must be at least 0.001 and below 0.5.
    pub fn new(rate: f64) -> Result<Self, String> {
        if (0.001..0.5).contains(&rate) {
            Ok(LeakRate(rate))
        } else {
            Err(format!(
                "The leak rate must be at least 0.001 and below 0.5, not {rate}!"
            ))
        }
    }

    /// The rate as a number.
    pub fn value(self) -> f64 {
        self.0
    }

    /// Whether an MNV with `mnv_reads` ALT reads explains a partial with `partial_reads` as leakage:
    /// `partial_reads < mnv_reads` and `P(X >= partial_reads) >= LEAK_TAIL` for
    /// `X ~ Binomial(mnv_reads + partial_reads, rate)`.
    pub fn explains(self, mnv_reads: i32, partial_reads: i32) -> bool {
        partial_reads < mnv_reads
            && binomial_upper_tail(
                i64::from(mnv_reads) + i64::from(partial_reads),
                i64::from(partial_reads),
                self.0,
            ) >= LEAK_TAIL
    }

    /// The fewest reads an MNV needs to explain a partial with `partial_reads`, found by bisection,
    /// as more MNV reads never stop explaining a partial.
    pub fn fewest_mnv_reads(self, partial_reads: i32) -> i32 {
        let mut low = partial_reads;
        let mut high = partial_reads + 1;
        while !self.explains(high, partial_reads) {
            low = high;
            high = high.saturating_mul(2);
        }
        while high - low > 1 {
            let middle = low + (high - low) / 2;
            if self.explains(middle, partial_reads) {
                high = middle;
            } else {
                low = middle;
            }
        }
        high
    }
}

/// `P(X >= successes)` for `X ~ Binomial(trials, p)`: 1 minus the lower tail, each term computed in
/// log space.
fn binomial_upper_tail(trials: i64, successes: i64, p: f64) -> f64 {
    if successes <= 0 {
        return 1.0;
    }
    if successes > trials {
        return 0.0;
    }
    let log_odds = p.ln() - (-p).ln_1p();
    let mut log_pmf = trials as f64 * (-p).ln_1p();
    let mut lower = log_pmf.exp();
    for k in 0..successes - 1 {
        log_pmf += ((trials - k) as f64 / (k + 1) as f64).ln() + log_odds;
        lower += log_pmf.exp();
    }
    (1.0 - lower).max(0.0)
}

/// How the partial rows of a region are folded: by their counts, or by their reads.
#[derive(Clone, Debug, PartialEq)]
pub enum LeakyMnvFix {
    /// Fold a partial whose reads are plausible leakage at this rate.
    Counts(LeakRate),
    /// Fold the reads of a partial that the BAMs VarDict read show to be the MNV's molecules or
    /// consistent with leakage, one BAM per sample.
    Reads {
        /// The indexed BAMs, the tumor (or only) sample's first.
        bams: Vec<PathBuf>,
        /// How VarDict read them.
        options: ReadOptions,
    },
}

/// What the reads say about a partial row under an MNV, per sample, the tumor (or only) sample
/// first: a call holding the partial's reads that go to the MNV, and one holding the reads that
/// stay.
#[derive(Clone, Debug, PartialEq)]
pub struct ReadFate<'a> {
    /// The reads that are the MNV's molecules or consistent with leakage, as a call of the partial.
    pub moved: Vec<TumorOnlyVariant<'a>>,
    /// The reads of another allele, as a call of the partial.
    pub kept: Vec<TumorOnlyVariant<'a>>,
}

/// The read fates of a region's rows, keyed by (MNV row, partial row).
pub type ReadFates<'a> = HashMap<(usize, usize), ReadFate<'a>>;

/// How a partial row is told from a separate allele.
#[derive(Clone, Copy, Debug)]
pub enum Rule<'r, 'a> {
    /// By its reads' count against the MNV's, at this leak rate.
    Counts(LeakRate),
    /// By what its reads' templates hold at the MNV's other changed bases.
    Reads(&'r ReadFates<'a>),
}

/// How a row is written once the partial rows of its region are folded.
#[derive(Clone, Debug, PartialEq)]
pub enum Fold<'a> {
    /// The row is written as VarDict wrote it.
    Keep,
    /// The row is a partial of the row at this index and is not written.
    Into(usize),
    /// The row is a partial of the row at this index that kept some reads: it is written with
    /// these calls, one per sample, holding the reads that stayed, and `moved` reads per sample
    /// went to the MNV.
    Split {
        /// The MNV row the reads went to.
        into: usize,
        /// The row's calls holding the reads that stayed, one per sample.
        calls: Vec<TumorOnlyVariant<'a>>,
        /// The reads that went to the MNV, one count per sample.
        moved: Vec<i32>,
    },
    /// The row is an MNV that absorbed partial rows, written with these calls, one per sample.
    Merged(Merged<'a>),
}

/// An MNV's calls after its partial rows are folded into it.
#[derive(Clone, Debug, PartialEq)]
pub struct Merged<'a> {
    /// The merged call of each sample, the tumor (or only) sample first.
    pub calls: Vec<TumorOnlyVariant<'a>>,
    /// The number of partial rows folded in.
    pub partials: usize,
    /// The ALT reads the partial rows added to each sample's AD[1].
    pub partial_reads: Vec<i32>,
}

/// Returns how to write each row of one VarDict region, given each row's calls, one per sample.
///
/// A row takes part when it is a substitution (an SNV or an MNV) with ALT reads in some sample. A
/// child is a sub-event of a parent when it lies within the parent's span and its ALT bases equal
/// the parent's over that span, so it changes no base the parent does not and changes the bases it
/// covers as the parent does. By counts, a parent can claim a child when, in every sample with
/// reads of the child, its own reads explain the child's as leakage at the leak rate; the test uses
/// only the parent's own reads, not those of the children it took, and a claimed child is folded
/// whole. By reads, a parent can claim a child when some sample moves at least one of the child's
/// reads to it; the child is folded whole when every sample moves every read, and is otherwise
/// written with the reads that stayed.
///
/// Rows are visited from the most to the least supported: by ALT reads summed over the samples,
/// then the longer row, the earlier start and the ALT allele, so the fold does not depend on the
/// order VarDict wrote the rows in. An MNV row that no root takes becomes a root, including one a
/// root could claim but has no DP room for at its first base; an SNV no root takes is written as it
/// is. A row that several roots can claim goes to the one that needs the fewest breaks to leave it,
/// one for each end of the child inside the root's span, since every break is a rare event in a
/// read; then to the root visited first, the one with the most ALT reads.
///
/// A child starting at a root's first base is taken only while the root's ALT reads and those of
/// the children it took there fit in that base's DP, as VarDict counts each read's base there once;
/// otherwise it goes to the next root that can claim it, so what happens to it depends on the
/// siblings visited before it. A row marked `written`, already written for an earlier region, can
/// be a root but is never claimed; a child such a root takes is folded into a record that is not
/// written again, so it is dropped.
pub fn fold_partials<'a>(
    rows: &[&[TumorOnlyVariant<'a>]],
    written: &[bool],
    rule: Rule<'_, 'a>,
) -> Vec<Fold<'a>> {
    let mut folds = vec![Fold::Keep; rows.len()];
    let mut nodes: Vec<usize> = (0..rows.len())
        .filter(|&row| {
            is_substitution(&rows[row][0]) && rows[row].iter().any(|call| call.alt_depth > 0)
        })
        .collect();
    if !nodes.iter().any(|&row| rows[row][0].ref_allele.len() > 1) {
        return folds;
    }
    nodes.sort_by_key(|&row| {
        let reads: i64 = rows[row].iter().map(|call| call.alt_depth as i64).sum();
        (
            -reads,
            -(rows[row][0].ref_allele.len() as i64),
            rows[row][0].start,
            rows[row][0].alt_allele,
            row,
        )
    });
    let mut roots: Vec<Root> = Vec::new();
    for &node in &nodes {
        let mut candidates: Vec<(usize, usize, Vec<TumorOnlyVariant<'a>>)> = roots
            .iter()
            .enumerate()
            .filter(|_| !written[node])
            .filter_map(|(rank, root)| {
                claimed(rule, rows, root.row, node)
                    .map(|claimed| (breaks(&rows[root.row][0], &rows[node][0]), rank, claimed))
            })
            .collect();
        if candidates.is_empty() {
            if rows[node][0].ref_allele.len() > 1 {
                roots.push(Root::new(node, rows[node]));
            }
            continue;
        }
        candidates.sort_by_key(|(breaks, rank, _)| (*breaks, *rank));
        let taker = candidates
            .into_iter()
            .find(|(_, rank, claimed)| roots[*rank].admits(rows, node, claimed));
        if let Some((_, rank, claimed)) = taker {
            roots[rank].take(rows, node, claimed);
        } else if rows[node][0].ref_allele.len() > 1 {
            roots.push(Root::new(node, rows[node]));
        }
    }
    for root in roots.iter_mut().filter(|root| !root.children.is_empty()) {
        root.children.sort_by_key(|(child, _)| *child);
        for (child, claimed) in &root.children {
            folds[*child] = match rule {
                Rule::Counts(_) => Fold::Into(root.row),
                Rule::Reads(fates) => {
                    let kept = &fates[&(root.row, *child)].kept;
                    if kept.iter().all(|call| call.alt_depth == 0) {
                        Fold::Into(root.row)
                    } else {
                        Fold::Split {
                            into: root.row,
                            calls: kept.clone(),
                            moved: claimed.iter().map(|call| call.alt_depth).collect(),
                        }
                    }
                }
            };
        }
        let children: Vec<&[TumorOnlyVariant<'a>]> = root
            .children
            .iter()
            .map(|(_, claimed)| claimed.as_slice())
            .collect();
        folds[root.row] = Fold::Merged(merge(rows[root.row], &children));
    }
    folds
}

/// The calls of a child's reads a parent row claims, one per sample, or `None` when it claims
/// none: by counts, the child is the parent's sub-event and, in every sample with reads of the
/// child, the parent's reads explain them as leakage, and every read is claimed; by reads, some
/// sample moves a read of the child to the parent, and the moved reads are claimed.
fn claimed<'a>(
    rule: Rule<'_, 'a>,
    rows: &[&[TumorOnlyVariant<'a>]],
    parent: usize,
    child: usize,
) -> Option<Vec<TumorOnlyVariant<'a>>> {
    let (parent_calls, child_calls) = (rows[parent], rows[child]);
    match rule {
        Rule::Counts(leak_rate) => {
            let explained = is_sub_event(&parent_calls[0], &child_calls[0])
                && parent_calls
                    .iter()
                    .zip(child_calls)
                    .all(|(p, c)| c.alt_depth <= 0 || leak_rate.explains(p.alt_depth, c.alt_depth));
            explained.then(|| child_calls.to_vec())
        }
        Rule::Reads(fates) => fates
            .get(&(parent, child))
            .filter(|fate| fate.moved.iter().any(|call| call.alt_depth > 0))
            .map(|fate| fate.moved.clone()),
    }
}

/// An MNV row that partial rows fold into, with the ALT reads it holds at its first base and the
/// calls it claimed of each child.
struct Root<'a> {
    row: usize,
    children: Vec<(usize, Vec<TumorOnlyVariant<'a>>)>,
    first_base_reads: Vec<i64>,
}

impl<'a> Root<'a> {
    /// A root without children yet.
    fn new(row: usize, calls: &[TumorOnlyVariant]) -> Self {
        Root {
            row,
            children: Vec::new(),
            first_base_reads: calls.iter().map(|call| call.alt_depth as i64).collect(),
        }
    }

    /// Whether a child's claimed reads fit beside this root's: always for a child starting after the
    /// root's first base, and otherwise only while every sample's ALT reads at that base stay within
    /// its DP.
    fn admits(
        &self,
        rows: &[&[TumorOnlyVariant]],
        child: usize,
        claimed: &[TumorOnlyVariant],
    ) -> bool {
        rows[child][0].start != rows[self.row][0].start
            || rows[self.row]
                .iter()
                .zip(claimed)
                .zip(&self.first_base_reads)
                .all(|((root, call), &reads)| {
                    call.alt_depth == 0 || reads + call.alt_depth as i64 <= root.depth as i64
                })
    }

    /// Takes a child's claimed reads into this root.
    fn take(
        &mut self,
        rows: &[&[TumorOnlyVariant]],
        child: usize,
        claimed: Vec<TumorOnlyVariant<'a>>,
    ) {
        if rows[child][0].start == rows[self.row][0].start {
            for (reads, call) in self.first_base_reads.iter_mut().zip(&claimed) {
                *reads += call.alt_depth as i64;
            }
        }
        self.children.push((child, claimed));
    }
}

/// The breaks a read of `parent` needs to leave only `child`: one for each end of the child that
/// falls inside the parent's span.
fn breaks(parent: &TumorOnlyVariant, child: &TumorOnlyVariant) -> usize {
    usize::from(child.start != parent.start) + usize::from(child.end != parent.end)
}

/// Whether a call is a substitution of one or more bases with a simple ALT allele.
pub(crate) fn is_substitution(call: &TumorOnlyVariant) -> bool {
    call.ref_allele.len() == call.alt_allele.len()
        && call.ref_allele != call.alt_allele
        && call
            .alt_allele
            .bytes()
            .all(|base| base.is_ascii_alphabetic())
}

/// Whether `child` lies within `parent`'s span and carries `parent`'s ALT bases over its own span.
pub(crate) fn is_sub_event(parent: &TumorOnlyVariant, child: &TumorOnlyVariant) -> bool {
    if child.ref_allele.len() >= parent.ref_allele.len() || child.start < parent.start {
        return false;
    }
    let offset = (child.start - parent.start) as usize;
    let span = offset..offset + child.ref_allele.len();
    parent.alt_allele.get(span.clone()) == Some(child.alt_allele)
        && parent.ref_allele.get(span) == Some(child.ref_allele)
}

/// Cuts a region's rows, sorted by start, into the clusters they fold within: a row starting past
/// the REF of every earlier row of its cluster starts the next, so an MNV's cluster holds every row
/// inside its REF.
#[derive(Default)]
pub(crate) struct Clusters {
    end: Option<u64>,
}

impl Clusters {
    /// Takes the next row by its start and REF length, returning whether it starts a cluster.
    pub(crate) fn starts_cluster(&mut self, start: u64, ref_length: usize) -> bool {
        let last = start
            .saturating_add(ref_length as u64)
            .saturating_sub(1)
            .max(start);
        match self.end {
            Some(end) if start <= end => {
                self.end = Some(end.max(last));
                false
            }
            _ => {
                self.end = Some(last);
                true
            }
        }
    }
}

/// Folds the partial rows' calls into the MNV's, one sample at a time.
fn merge<'a>(parent: &[TumorOnlyVariant<'a>], children: &[&[TumorOnlyVariant<'a>]]) -> Merged<'a> {
    let calls: Vec<TumorOnlyVariant<'a>> = parent
        .iter()
        .enumerate()
        .map(|(sample, call)| {
            let sample_children: Vec<&TumorOnlyVariant<'a>> =
                children.iter().map(|row| &row[sample]).collect();
            merge_sample(call, &sample_children)
        })
        .collect();
    let partial_reads = calls
        .iter()
        .zip(parent)
        .map(|(merged, call)| merged.alt_depth - call.alt_depth)
        .collect();
    Merged {
        calls,
        partials: children.len(),
        partial_reads,
    }
}

/// Folds one sample's partial calls into the MNV's call.
///
/// The reads of a partial starting at the MNV's first base are in that base's DP already, as VarDict
/// counts each read's base there into one allele and into the depth, so only AD grows; a partial
/// starting later adds its reads to AD and DP alike, as VarDict's own MNV merging does.
fn merge_sample<'a>(
    parent: &TumorOnlyVariant<'a>,
    children: &[&TumorOnlyVariant<'a>],
) -> TumorOnlyVariant<'a> {
    let rows: Vec<&TumorOnlyVariant> = std::iter::once(parent)
        .chain(children.iter().copied())
        .collect();
    let later_reads: i32 = children
        .iter()
        .filter(|child| child.start != parent.start)
        .map(|child| child.alt_depth)
        .sum();
    let added_reads: i32 = children.iter().map(|child| child.alt_depth).sum();
    let alt_depth = parent.alt_depth + added_reads;
    let depth = parent.depth + later_reads;
    let alt_forward: i32 = rows.iter().map(|row| row.alt_forward).sum();
    let alt_reverse: i32 = rows.iter().map(|row| row.alt_reverse).sum();
    let with_reads: Vec<&TumorOnlyVariant> = rows
        .iter()
        .copied()
        .filter(|row| row.alt_depth > 0)
        .collect();
    let mean = |value: &dyn Fn(&TumorOnlyVariant) -> f32| {
        read_weighted_mean(&with_reads, value).unwrap_or(value(parent))
    };
    let positions_vary = with_reads
        .iter()
        .any(|row| row.stdev_position_in_read != 0.0)
        || with_reads
            .windows(2)
            .any(|pair| pair[0].mean_position_in_read != pair[1].mean_position_in_read);
    let high_quality_variant_reads = if rows
        .iter()
        .any(|row| row.high_quality_variant_reads.is_missing())
    {
        i32::missing()
    } else {
        rows.iter().map(|row| row.high_quality_variant_reads).sum()
    };
    let realigned_reads = parent.af_adjusted as f64 * parent.depth as f64 + added_reads as f64;
    TumorOnlyVariant {
        depth,
        alt_depth,
        alt_forward,
        alt_reverse,
        mean_position_in_read: mean(&|row| row.mean_position_in_read),
        stdev_position_in_read: if positions_vary { 1.0 } else { 0.0 },
        base_quality_mean: mean(&|row| row.base_quality_mean),
        mean_mapping_quality: mean(&|row| row.mean_mapping_quality),
        mean_mismatches_in_reads: mean(&|row| row.mean_mismatches_in_reads),
        high_quality_variant_reads,
        af_adjusted: if depth > 0 {
            (realigned_reads / depth as f64) as f32
        } else {
            0.0
        },
        ..parent.clone()
    }
}

/// The mean of a per-read average over the rows' ALT reads, each row weighted by its ALT reads as
/// VarDict sums these averages when it merges reads into an allele; `None` without ALT reads.
fn read_weighted_mean(
    rows: &[&TumorOnlyVariant],
    value: &dyn Fn(&TumorOnlyVariant) -> f32,
) -> Option<f32> {
    let reads: f64 = rows.iter().map(|row| row.alt_depth as f64).sum();
    let total: f64 = rows
        .iter()
        .map(|row| row.alt_depth as f64 * value(row) as f64)
        .sum();
    (reads > 0.0).then(|| (total / reads) as f32)
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::rstest;

    use super::*;
    use crate::fisher::fisher_exact_two_sided;

    /// A call of `alt_depth` ALT reads, split evenly by strand, with fixed per-read summaries.
    fn call(
        start: u64,
        ref_allele: &'static str,
        alt_allele: &'static str,
        alt_depth: i32,
        depth: i32,
    ) -> TumorOnlyVariant<'static> {
        TumorOnlyVariant {
            start,
            end: start + ref_allele.len() as u64 - 1,
            ref_allele,
            alt_allele,
            depth,
            alt_depth,
            alt_forward: alt_depth / 2,
            alt_reverse: alt_depth - alt_depth / 2,
            ref_forward: 50,
            ref_reverse: 50,
            mean_position_in_read: 20.0,
            stdev_position_in_read: 1.0,
            base_quality_mean: 80.0,
            mean_mapping_quality: 60.0,
            mean_mismatches_in_reads: 2.0,
            high_quality_variant_reads: alt_depth,
            variant_type: if ref_allele.len() == 1 {
                "SNV"
            } else {
                "Complex"
            },
            ..Default::default()
        }
    }

    /// Folds tumor-only rows.
    fn fold(rows: &[TumorOnlyVariant<'static>]) -> Vec<Fold<'static>> {
        let rows: Vec<&[TumorOnlyVariant]> = rows.iter().map(std::slice::from_ref).collect();
        fold_partials(
            &rows,
            &vec![false; rows.len()],
            Rule::Counts(LeakRate::default()),
        )
    }

    /// The rows each row folds into, `None` for a row that is written.
    fn targets(folds: &[Fold]) -> Vec<Option<usize>> {
        folds
            .iter()
            .map(|fold| match fold {
                Fold::Into(root) => Some(*root),
                _ => None,
            })
            .collect()
    }

    /// The merged call of the tumor (or only) sample of a row.
    fn merged<'a>(folds: &'a [Fold], row: usize) -> &'a Merged<'a> {
        match &folds[row] {
            Fold::Merged(merged) => merged,
            other => panic!("row {row} is not merged: {other:?}"),
        }
    }

    #[test]
    fn test_a_partial_at_the_first_base_adds_to_ad_but_not_dp() {
        let folds = fold(&[call(10, "GG", "CA", 100, 200), call(10, "G", "C", 5, 200)]);
        assert_eq!(targets(&folds), vec![None, Some(0)]);
        let merged = merged(&folds, 0);
        assert_eq!(merged.partials, 1);
        assert_eq!(merged.partial_reads, vec![5]);
        assert_eq!(
            (merged.calls[0].alt_depth, merged.calls[0].depth),
            (105, 200)
        );
    }

    #[test]
    fn test_a_shorter_mnv_left_by_one_no_call_folds_into_the_mnv() {
        // Reference CGTGGTT; three reads CGAAATT and one CGAANTT, whose N stops the MNV a base short.
        let folds = fold(&[call(3, "TGG", "AAA", 3, 4), call(3, "TG", "AA", 1, 4)]);
        assert_eq!(targets(&folds), vec![None, Some(0)]);
        let merged = &merged(&folds, 0).calls[0];
        assert_eq!((merged.alt_depth, merged.depth), (4, 4));
        assert_eq!(merged.af_value(), 1.0);
    }

    #[test]
    fn test_an_anchored_shorter_mnv_folds_into_the_anchored_mnv() {
        let folds = fold(&[call(2, "GTGG", "GAAA", 3, 4), call(2, "GTG", "GAA", 1, 4)]);
        assert_eq!(targets(&folds), vec![None, Some(0)]);
        let merged = &merged(&folds, 0).calls[0];
        assert_eq!((merged.alt_depth, merged.depth), (4, 4));
    }

    #[test]
    fn test_partials_of_a_gapped_mnv_fold_at_its_changed_bases() {
        let folds = fold(&[
            call(20, "TTGG", "CTGA", 80, 200),
            call(20, "T", "C", 4, 200),
            call(23, "G", "A", 3, 203),
            call(21, "T", "A", 2, 201),
            call(22, "G", "C", 2, 202),
        ]);
        assert_eq!(targets(&folds), vec![None, Some(0), Some(0), None, None]);
        let merged = &merged(&folds, 0).calls[0];
        assert_eq!((merged.alt_depth, merged.depth), (87, 203));
    }

    #[test]
    fn test_a_gapped_mnv_folds_into_a_longer_mnv_that_carries_it() {
        let folds = fold(&[
            call(20, "TTGGA", "CTGAC", 80, 200),
            call(20, "TTGG", "CTGA", 4, 200),
            call(23, "GA", "AT", 3, 203),
        ]);
        assert_eq!(targets(&folds), vec![None, Some(0), None]);
    }

    #[test]
    fn test_a_partial_at_a_later_base_adds_to_ad_and_dp() {
        let folds = fold(&[call(10, "GG", "CA", 100, 200), call(11, "G", "A", 5, 210)]);
        assert_eq!(targets(&folds), vec![None, Some(0)]);
        let merged = merged(&folds, 0);
        assert_eq!(
            (merged.calls[0].alt_depth, merged.calls[0].depth),
            (105, 205)
        );
    }

    #[test]
    fn test_a_partial_folds_whether_it_comes_before_or_after_its_mnv() {
        let folds = fold(&[call(11, "G", "A", 5, 210), call(10, "GG", "CA", 100, 200)]);
        assert_eq!(targets(&folds), vec![Some(1), None]);
    }

    #[test]
    fn test_a_strong_snv_under_a_weak_mnv_stays_apart() {
        let folds = fold(&[call(10, "GG", "CA", 5, 200), call(10, "G", "C", 49, 200)]);
        assert_eq!(folds, vec![Fold::Keep, Fold::Keep]);
    }

    #[rstest]
    #[case::other_alt_base(call(10, "G", "T", 5, 200))]
    #[case::outside_the_span(call(12, "G", "A", 5, 200))]
    #[case::straddling_the_span(call(11, "GG", "AT", 5, 200))]
    #[case::reference_inside_the_span(call(10, "GGG", "CGA", 5, 200))]
    fn test_a_row_that_is_not_a_sub_event_stays_apart(#[case] other: TumorOnlyVariant<'static>) {
        let mnv = if other.ref_allele.len() == 3 {
            call(10, "GGG", "CAA", 100, 200)
        } else {
            call(10, "GG", "CA", 100, 200)
        };
        assert_eq!(fold(&[mnv, other]), vec![Fold::Keep, Fold::Keep]);
    }

    #[rstest]
    #[case(2, 1, Some(0))]
    #[case(1, 1, None)]
    #[case(3, 2, None)]
    #[case(9, 5, None)]
    #[case(5, 2, None)]
    #[case(6, 2, Some(0))]
    #[case(100, 9, Some(0))]
    #[case(100, 10, None)]
    fn test_a_partial_folds_only_when_its_reads_are_plausible_leakage(
        #[case] mnv_reads: i32,
        #[case] snv_reads: i32,
        #[case] expected: Option<usize>,
    ) {
        let folds = fold(&[
            call(10, "GG", "CA", mnv_reads, 200),
            call(10, "G", "C", snv_reads, 200),
        ]);
        assert_eq!(targets(&folds)[1], expected);
    }

    #[rstest]
    #[case(2, 1, 0.03, 0.0591)]
    #[case(4, 1, 0.03, 0.11471)]
    #[case(5, 2, 0.03, 0.00847)]
    #[case(10, 0, 0.03, 1.0)]
    #[case(10, 11, 0.03, 0.0)]
    #[case(1041, 41, 0.03, 0.05057)]
    #[case(1042, 42, 0.03, 0.03595)]
    #[case(3, 3, 0.5, 0.125)]
    fn test_the_binomial_upper_tail(
        #[case] trials: i64,
        #[case] successes: i64,
        #[case] p: f64,
        #[case] expected: f64,
    ) {
        let tail = binomial_upper_tail(trials, successes, p);
        assert!((tail - expected).abs() < 5e-5, "{tail} is not {expected}");
    }

    /// The most reads a partial may have and still be explained by an MNV with `mnv_reads`.
    fn bound(rate: LeakRate, mnv_reads: i32) -> i32 {
        (0..mnv_reads)
            .rev()
            .find(|&reads| rate.explains(mnv_reads, reads))
            .unwrap_or(0)
    }

    #[test]
    fn test_the_default_leak_rate_bounds_a_partials_reads_by_the_mnvs() {
        let rate = LeakRate::default();
        let bounds: Vec<(i32, i32)> = [1, 2, 5, 6, 10, 15, 20, 25, 50, 100, 200, 1000]
            .into_iter()
            .map(|reads| (reads, bound(rate, reads)))
            .collect();
        assert_eq!(
            bounds,
            vec![
                (1, 0),
                (2, 1),
                (5, 1),
                (6, 2),
                (10, 2),
                (15, 3),
                (20, 3),
                (25, 4),
                (50, 6),
                (100, 9),
                (200, 16),
                (1000, 65)
            ]
        );
    }

    #[rstest]
    #[case(1, 1)]
    #[case(2, 2)]
    #[case(5, 7)]
    fn test_a_partial_with_at_least_the_mnvs_reads_is_never_leakage(
        #[case] mnv_reads: i32,
        #[case] partial_reads: i32,
    ) {
        let rate = LeakRate::new(0.49).unwrap();
        assert!(!rate.explains(mnv_reads, partial_reads));
        assert!(rate.explains(mnv_reads + partial_reads, partial_reads));
    }

    #[rstest]
    #[case(0.0)]
    #[case(-0.1)]
    #[case(0.0009)]
    #[case(0.5)]
    #[case(1.0)]
    #[case(f64::NAN)]
    fn test_a_leak_rate_below_the_floor_or_at_the_ceiling_is_refused(#[case] rate: f64) {
        assert!(LeakRate::new(rate).is_err());
    }

    #[rstest]
    #[case(0.001)]
    #[case(0.05)]
    #[case(0.49)]
    fn test_the_fewest_mnv_reads_is_where_a_count_up_first_explains_the_partial(#[case] rate: f64) {
        let rate = LeakRate::new(rate).unwrap();
        for partial_reads in [1, 2, 5, 10] {
            let fewest = rate.fewest_mnv_reads(partial_reads);
            assert!(rate.explains(fewest, partial_reads));
            assert!(!rate.explains(fewest - 1, partial_reads));
            if rate.value() > 0.01 {
                let counted = (partial_reads + 1..)
                    .find(|&reads| rate.explains(reads, partial_reads))
                    .unwrap();
                assert_eq!(fewest, counted);
            }
        }
    }

    #[test]
    fn test_a_higher_leak_rate_folds_more_reads() {
        let rows = [call(10, "GG", "CA", 100, 200), call(10, "G", "C", 9, 200)];
        let calls: Vec<&[TumorOnlyVariant]> = rows.iter().map(std::slice::from_ref).collect();
        let strict = fold_partials(
            &calls,
            &[false, false],
            Rule::Counts(LeakRate::new(0.03).unwrap()),
        );
        let loose = fold_partials(
            &calls,
            &[false, false],
            Rule::Counts(LeakRate::new(0.05).unwrap()),
        );
        assert_eq!(targets(&strict), vec![None, None]);
        assert_eq!(targets(&loose), vec![None, Some(0)]);
    }

    #[test]
    fn test_partials_do_not_share_a_leak_budget() {
        let folds = fold(&[
            call(10, "GGG", "CAT", 100, 200),
            call(10, "GG", "CA", 6, 200),
            call(10, "G", "C", 6, 200),
            call(11, "G", "A", 6, 206),
        ]);
        assert_eq!(targets(&folds), vec![None, Some(0), Some(0), Some(0)]);
    }

    #[test]
    fn test_a_partial_of_a_claimed_sub_mnv_folds_into_the_root_that_explains_it() {
        // The 3-read sub-MNV cannot explain a 2-read SNV, but the root it leaked from can.
        let folds = fold(&[
            call(10, "GGGG", "CATA", 200, 400),
            call(10, "GGG", "CAT", 3, 400),
            call(10, "G", "C", 2, 400),
        ]);
        assert_eq!(targets(&folds), vec![None, Some(0), Some(0)]);
    }

    #[test]
    fn test_a_refused_sub_mnv_is_a_root_and_competes_for_a_partial() {
        // The sub-MNV has too many reads to be leakage, so it is a root; the SNV at its first base
        // needs one break to leave either MNV and goes to the one with more reads.
        let folds = fold(&[
            call(10, "GGGG", "CATA", 100, 400),
            call(10, "GGG", "CAT", 40, 400),
            call(10, "G", "C", 3, 400),
        ]);
        assert_eq!(targets(&folds), vec![None, None, Some(0)]);
    }

    #[test]
    fn test_chained_partials_fold_into_the_longest_mnv() {
        let folds = fold(&[
            call(10, "GGG", "CAT", 100, 200),
            call(10, "GG", "CA", 5, 200),
            call(10, "G", "C", 3, 200),
            call(11, "G", "A", 2, 205),
            call(12, "G", "T", 4, 206),
        ]);
        assert_eq!(
            targets(&folds),
            vec![None, Some(0), Some(0), Some(0), Some(0)]
        );
        let merged = merged(&folds, 0);
        assert_eq!(merged.partials, 4);
        assert_eq!(
            (merged.calls[0].alt_depth, merged.calls[0].depth),
            (114, 206)
        );
    }

    #[test]
    fn test_a_partial_two_mnvs_share_goes_to_the_better_supported() {
        let folds = fold(&[
            call(10, "GGG", "CTT", 40, 200),
            call(10, "GGG", "CAA", 100, 200),
            call(10, "G", "C", 3, 200),
        ]);
        assert_eq!(targets(&folds), vec![None, None, Some(1)]);
    }

    #[test]
    fn test_a_partial_goes_to_the_mnv_that_needs_the_fewest_breaks() {
        let folds = fold(&[
            call(9, "GGG", "TCA", 100, 200),
            call(10, "GG", "CT", 40, 200),
            call(10, "G", "C", 3, 200),
        ]);
        assert_eq!(targets(&folds), vec![None, None, Some(1)]);
    }

    #[test]
    fn test_a_partial_at_the_first_base_folds_only_while_it_fits_in_dp() {
        let folds = fold(&[
            call(10, "GGG", "CAT", 100, 105),
            call(10, "GG", "CA", 4, 105),
            call(10, "G", "C", 3, 105),
        ]);
        assert_eq!(targets(&folds), vec![None, Some(0), None]);
        let merged = merged(&folds, 0);
        assert!(merged.calls[0].alt_depth <= merged.calls[0].depth);
    }

    #[test]
    fn test_an_mnv_no_root_has_room_for_takes_its_own_partials() {
        let folds = fold(&[
            call(10, "GGG", "CAT", 100, 103),
            call(10, "GG", "CA", 4, 103),
            call(11, "G", "A", 1, 103),
        ]);
        assert_eq!(targets(&folds), vec![None, None, Some(1)]);
    }

    #[test]
    fn test_the_order_of_the_rows_does_not_change_the_fold() {
        let rows = [
            call(10, "GGG", "CAT", 100, 200),
            call(12, "G", "T", 4, 206),
            call(10, "G", "C", 3, 200),
            call(11, "G", "A", 2, 205),
        ];
        let reversed: Vec<TumorOnlyVariant> = rows.iter().rev().cloned().collect();
        let forward = fold(&rows);
        let backward = fold(&reversed);
        assert_eq!(merged(&forward, 0), merged(&backward, 3));
    }

    #[test]
    fn test_two_mnvs_tied_on_reads_take_a_partial_in_alt_allele_order_whatever_the_input_order() {
        let rows = [
            call(10, "GG", "CT", 50, 200),
            call(10, "GG", "CA", 50, 200),
            call(10, "G", "C", 2, 200),
        ];
        let swapped = [rows[1].clone(), rows[0].clone(), rows[2].clone()];
        assert_eq!(targets(&fold(&rows)), vec![None, None, Some(1)]);
        assert_eq!(targets(&fold(&swapped)), vec![None, None, Some(0)]);
    }

    #[test]
    fn test_a_row_written_before_is_not_claimed_but_can_be_a_root() {
        let rows = [
            call(10, "GG", "CA", 100, 200),
            call(10, "G", "C", 5, 200),
            call(11, "G", "A", 4, 205),
        ];
        let calls: Vec<&[TumorOnlyVariant]> = rows.iter().map(std::slice::from_ref).collect();
        let folds = fold_partials(
            &calls,
            &[false, true, false],
            Rule::Counts(LeakRate::default()),
        );
        assert_eq!(targets(&folds), vec![None, None, Some(0)]);
        let folds = fold_partials(
            &calls,
            &[true, false, false],
            Rule::Counts(LeakRate::default()),
        );
        assert_eq!(targets(&folds), vec![None, Some(0), Some(0)]);
    }

    #[test]
    fn test_without_an_mnv_every_row_is_kept() {
        let folds = fold(&[call(10, "G", "C", 100, 200), call(11, "G", "A", 5, 200)]);
        assert_eq!(folds, vec![Fold::Keep, Fold::Keep]);
    }

    #[test]
    fn test_the_merged_call_combines_the_rows_summaries() {
        let mnv = TumorOnlyVariant {
            alt_forward: 65,
            alt_reverse: 30,
            mean_position_in_read: 20.0,
            stdev_position_in_read: 0.0,
            base_quality_mean: 90.0,
            mean_mapping_quality: 60.0,
            mean_mismatches_in_reads: 1.0,
            af_adjusted: 0.05,
            ..call(10, "GG", "CA", 95, 200)
        };
        let partial = TumorOnlyVariant {
            alt_forward: 5,
            alt_reverse: 0,
            mean_position_in_read: 2.0,
            stdev_position_in_read: 0.0,
            base_quality_mean: 70.0,
            mean_mapping_quality: 40.0,
            mean_mismatches_in_reads: 3.0,
            ..call(10, "G", "C", 5, 200)
        };
        let folds = fold(&[mnv, partial]);
        let merged = &merged(&folds, 0).calls[0];
        assert_eq!((merged.alt_forward, merged.alt_reverse), (70, 30));
        assert_eq!(merged.mean_position_in_read, 19.1);
        assert_eq!(merged.stdev_position_in_read, 1.0);
        assert_eq!(merged.base_quality_mean, 89.0);
        assert_eq!(merged.mean_mapping_quality, 59.0);
        assert_eq!(merged.mean_mismatches_in_reads, 1.1);
        assert_eq!(merged.high_quality_variant_reads, 100);
        assert_eq!(merged.af_adjusted, 0.075);
        let expected = fisher_exact_two_sided(50, 50, 70, 30) as f32;
        assert_eq!(merged.strand_bias_fisher_p_value(), expected);
    }

    #[test]
    fn test_read_positions_vary_only_when_the_rows_differ() {
        let same = |start, ref_allele, alt_allele, reads| TumorOnlyVariant {
            stdev_position_in_read: 0.0,
            ..call(start, ref_allele, alt_allele, reads, 200)
        };
        let folds = fold(&[same(10, "GG", "CA", 90), same(10, "G", "C", 3)]);
        assert_eq!(merged(&folds, 0).calls[0].stdev_position_in_read, 0.0);
    }

    #[test]
    fn test_a_missing_hicnt_makes_the_merged_hicnt_missing() {
        let partial = TumorOnlyVariant {
            high_quality_variant_reads: i32::missing(),
            ..call(10, "G", "C", 3, 200)
        };
        let folds = fold(&[call(10, "GG", "CA", 90, 200), partial]);
        assert!(
            merged(&folds, 0).calls[0]
                .high_quality_variant_reads
                .is_missing()
        );
    }

    /// Folds tumor-normal rows given as (tumor, normal) pairs.
    fn fold_pairs(rows: &[[TumorOnlyVariant<'static>; 2]]) -> Vec<Fold<'static>> {
        let rows: Vec<&[TumorOnlyVariant]> = rows.iter().map(|row| &row[..]).collect();
        fold_partials(
            &rows,
            &vec![false; rows.len()],
            Rule::Counts(LeakRate::default()),
        )
    }

    #[test]
    fn test_a_tumor_normal_partial_folds_in_both_samples() {
        let folds = fold_pairs(&[
            [
                call(10, "GG", "CA", 100, 200),
                call(10, "GG", "CA", 50, 100),
            ],
            [call(11, "G", "A", 4, 210), call(11, "G", "A", 2, 105)],
        ]);
        assert_eq!(targets(&folds), vec![None, Some(0)]);
        let merged = merged(&folds, 0);
        assert_eq!(merged.partial_reads, vec![4, 2]);
        let counts: Vec<(i32, i32)> = merged
            .calls
            .iter()
            .map(|call| (call.alt_depth, call.depth))
            .collect();
        assert_eq!(counts, vec![(104, 204), (52, 102)]);
    }

    #[test]
    fn test_a_tumor_normal_partial_must_be_leakage_in_every_sample_with_its_reads() {
        let folds = fold_pairs(&[
            [call(10, "GG", "CA", 100, 200), call(10, "GG", "CA", 0, 100)],
            [call(10, "G", "C", 4, 200), call(10, "G", "C", 3, 100)],
        ]);
        assert_eq!(folds, vec![Fold::Keep, Fold::Keep]);
    }

    #[test]
    fn test_a_sample_without_reads_of_the_partial_does_not_object() {
        let folds = fold_pairs(&[
            [call(10, "GG", "CA", 100, 200), call(10, "GG", "CA", 0, 100)],
            [call(10, "G", "C", 4, 200), call(10, "G", "C", 0, 100)],
        ]);
        assert_eq!(targets(&folds), vec![None, Some(0)]);
    }

    /// A fate of one sample: `moved` and `kept` reads of a row's call.
    fn fate(row: &TumorOnlyVariant<'static>, moved: i32, kept: i32) -> ReadFate<'static> {
        let with = |reads: i32| TumorOnlyVariant {
            alt_depth: reads,
            alt_forward: reads,
            alt_reverse: 0,
            high_quality_variant_reads: reads,
            ..row.clone()
        };
        ReadFate {
            moved: vec![with(moved)],
            kept: vec![with(kept)],
        }
    }

    /// Folds tumor-only rows by read fates keyed by (MNV row, partial row).
    fn fold_by_reads(
        rows: &[TumorOnlyVariant<'static>],
        fates: &[((usize, usize), ReadFate<'static>)],
    ) -> Vec<Fold<'static>> {
        let rows: Vec<&[TumorOnlyVariant]> = rows.iter().map(std::slice::from_ref).collect();
        let fates: ReadFates = fates.iter().cloned().collect();
        fold_partials(&rows, &vec![false; rows.len()], Rule::Reads(&fates))
    }

    #[test]
    fn test_by_reads_a_partial_whose_every_read_moves_is_folded_whole() {
        let rows = [call(10, "GG", "CA", 20, 200), call(10, "G", "C", 15, 200)];
        let folds = fold_by_reads(&rows, &[((0, 1), fate(&rows[1], 15, 0))]);
        assert_eq!(targets(&folds), vec![None, Some(0)]);
        let merged = merged(&folds, 0);
        assert_eq!(
            (merged.calls[0].alt_depth, merged.calls[0].depth),
            (35, 200)
        );
        assert_eq!(merged.partial_reads, vec![15]);
    }

    #[test]
    fn test_by_reads_a_partial_that_keeps_reads_is_written_with_them() {
        let rows = [call(10, "GG", "CA", 100, 200), call(11, "G", "A", 5, 210)];
        let folds = fold_by_reads(&rows, &[((0, 1), fate(&rows[1], 2, 3))]);
        match &folds[1] {
            Fold::Split { into, calls, moved } => {
                assert_eq!((*into, moved), (0, &vec![2]));
                assert_eq!((calls[0].alt_depth, calls[0].depth), (3, 210));
            }
            other => panic!("not split: {other:?}"),
        }
        let merged = merged(&folds, 0);
        assert_eq!(
            (merged.calls[0].alt_depth, merged.calls[0].depth),
            (102, 202)
        );
        assert_eq!(
            (merged.partials, merged.partial_reads.clone()),
            (1, vec![2])
        );
    }

    #[test]
    fn test_by_reads_a_partial_with_no_read_moving_is_kept_and_claims_nothing() {
        let rows = [call(10, "GG", "CA", 100, 200), call(10, "G", "C", 2, 200)];
        let folds = fold_by_reads(&rows, &[((0, 1), fate(&rows[1], 0, 2))]);
        assert_eq!(folds, vec![Fold::Keep, Fold::Keep]);
    }

    #[test]
    fn test_by_reads_the_count_rule_does_not_apply() {
        // Counts would refuse 40 reads under 100, and would fold 1 read under 100.
        let rows = [
            call(10, "GG", "CA", 100, 200),
            call(10, "G", "C", 40, 200),
            call(11, "G", "A", 1, 201),
        ];
        let folds = fold_by_reads(
            &rows,
            &[
                ((0, 1), fate(&rows[1], 40, 0)),
                ((0, 2), fate(&rows[2], 0, 1)),
            ],
        );
        assert_eq!(targets(&folds), vec![None, Some(0), None]);
        assert_eq!(folds[2], Fold::Keep);
    }

    #[test]
    fn test_by_reads_a_partial_at_the_first_base_fits_in_dp_with_its_moved_reads() {
        let rows = [call(10, "GG", "CA", 100, 103), call(10, "G", "C", 5, 103)];
        let fits = fold_by_reads(&rows, &[((0, 1), fate(&rows[1], 3, 2))]);
        assert!(matches!(fits[1], Fold::Split { .. }));
        let overflows = fold_by_reads(&rows, &[((0, 1), fate(&rows[1], 4, 1))]);
        assert_eq!(overflows, vec![Fold::Keep, Fold::Keep]);
    }

    #[test]
    fn test_by_reads_a_tumor_normal_partial_is_folded_whole_only_when_both_samples_move_every_read()
    {
        let rows: Vec<[TumorOnlyVariant; 2]> = vec![
            [
                call(10, "GG", "CA", 100, 200),
                call(10, "GG", "CA", 50, 100),
            ],
            [call(11, "G", "A", 4, 210), call(11, "G", "A", 2, 105)],
        ];
        let calls: Vec<&[TumorOnlyVariant]> = rows.iter().map(|row| &row[..]).collect();
        let fates = |normal_moved: i32, normal_kept: i32| {
            let tumor = fate(&rows[1][0], 4, 0);
            let normal = fate(&rows[1][1], normal_moved, normal_kept);
            ReadFates::from([(
                (0usize, 1usize),
                ReadFate {
                    moved: vec![tumor.moved[0].clone(), normal.moved[0].clone()],
                    kept: vec![tumor.kept[0].clone(), normal.kept[0].clone()],
                },
            )])
        };
        let whole = fold_partials(&calls, &[false, false], Rule::Reads(&fates(2, 0)));
        assert_eq!(targets(&whole), vec![None, Some(0)]);
        assert_eq!(merged(&whole, 0).partial_reads, vec![4, 2]);
        let split = fold_partials(&calls, &[false, false], Rule::Reads(&fates(0, 2)));
        match &split[1] {
            Fold::Split { calls, moved, .. } => {
                assert_eq!(moved, &vec![4, 0]);
                assert_eq!(
                    calls.iter().map(|call| call.alt_depth).collect::<Vec<_>>(),
                    vec![0, 2]
                );
            }
            other => panic!("not split: {other:?}"),
        }
        assert_eq!(merged(&split, 0).partial_reads, vec![4, 0]);
    }

    #[test]
    fn test_by_reads_a_row_written_before_is_not_claimed() {
        let rows = [call(10, "GG", "CA", 100, 200), call(10, "G", "C", 5, 200)];
        let calls: Vec<&[TumorOnlyVariant]> = rows.iter().map(std::slice::from_ref).collect();
        let fates = ReadFates::from([((0usize, 1usize), fate(&rows[1], 5, 0))]);
        let folds = fold_partials(&calls, &[false, true], Rule::Reads(&fates));
        assert_eq!(folds, vec![Fold::Keep, Fold::Keep]);
    }

    #[test]
    fn test_a_cluster_holds_the_rows_starting_inside_an_earlier_rows_ref() {
        let mut clusters = Clusters::default();
        let rows = [
            (10, 3),
            (12, 1),
            (13, 2),
            (15, 1),
            (16, 1),
            (u64::MAX, 3),
            (u64::MAX, 1),
        ];
        let starts: Vec<bool> = rows
            .into_iter()
            .map(|(start, length)| clusters.starts_cluster(start, length))
            .collect();
        assert_eq!(starts, vec![true, false, true, true, true, true, false]);
    }

    /// A stream of pseudo-random numbers, the same for the same seed.
    struct Lcg(u64);

    impl Lcg {
        fn below(&mut self, n: u64) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (self.0 >> 33) % n
        }

        /// A base other than `base`.
        fn other_than(&mut self, base: u8) -> u8 {
            let index = b"ACGT".iter().position(|&b| b == base).unwrap();
            b"ACGT"[(index + 1 + self.below(3) as usize) % 4]
        }
    }

    /// A region's rows over a short reference: MNVs, partials of them with a few of their reads, and
    /// unrelated SNVs, sorted by start as they come from the sort.
    fn random_rows(rng: &mut Lcg) -> Vec<TumorOnlyVariant<'static>> {
        const REFERENCE: &[u8] = b"ACGTTGCAAGCTGACCTGATCGGATCCATGCAGTCAGTACGATCGATGCATGCTAGC";
        let allele = |bases: &[u8]| -> &'static str {
            Box::leak(String::from_utf8(bases.to_vec()).unwrap().into_boxed_str())
        };
        let mut rows = Vec::new();
        for _ in 0..1 + rng.below(6) {
            let length = 2 + rng.below(4) as usize;
            let offset = rng.below((REFERENCE.len() - length) as u64) as usize;
            let reference = &REFERENCE[offset..offset + length];
            let alt: Vec<u8> = reference.iter().map(|&base| rng.other_than(base)).collect();
            let reads = 20 + rng.below(200) as i32;
            let start = 100 + offset as u64;
            rows.push(call(start, allele(reference), allele(&alt), reads, 1000));
            for _ in 0..rng.below(4) {
                let from = rng.below(length as u64) as usize;
                let to = from + 1 + rng.below((length - from) as u64) as usize;
                let partial_reads = 1 + rng.below(reads as u64 / 3) as i32;
                let (reference, alt) = (allele(&reference[from..to]), allele(&alt[from..to]));
                rows.push(call(
                    start + from as u64,
                    reference,
                    alt,
                    partial_reads,
                    1000,
                ));
            }
        }
        for _ in 0..rng.below(5) {
            let offset = rng.below(REFERENCE.len() as u64) as usize;
            let alt = rng.other_than(REFERENCE[offset]);
            let (reference, alt) = (allele(&REFERENCE[offset..=offset]), allele(&[alt]));
            rows.push(call(
                100 + offset as u64,
                reference,
                alt,
                1 + rng.below(50) as i32,
                1000,
            ));
        }
        rows.sort_by_key(|row| row.start);
        rows
    }

    #[test]
    fn test_folding_a_regions_clusters_one_at_a_time_folds_it_as_a_whole() {
        let mut rng = Lcg(7);
        let mut folded = 0;
        for _ in 0..500 {
            let rows = random_rows(&mut rng);
            let written: Vec<bool> = rows.iter().map(|_| rng.below(10) == 0).collect();
            let calls: Vec<&[TumorOnlyVariant]> = rows.iter().map(std::slice::from_ref).collect();
            let rule = Rule::Counts(LeakRate::default());
            let whole = fold_partials(&calls, &written, rule);
            let mut clusters = Clusters::default();
            let mut bounds: Vec<usize> = (0..rows.len())
                .filter(|&row| clusters.starts_cluster(rows[row].start, rows[row].ref_allele.len()))
                .collect();
            bounds.push(rows.len());
            let mut by_cluster = Vec::new();
            for cluster in bounds.windows(2) {
                let (first, last) = (cluster[0], cluster[1]);
                let folds = fold_partials(&calls[first..last], &written[first..last], rule);
                by_cluster.extend(folds.into_iter().map(|fold| match fold {
                    Fold::Into(root) => Fold::Into(root + first),
                    Fold::Split { into, calls, moved } => Fold::Split {
                        into: into + first,
                        calls,
                        moved,
                    },
                    other => other,
                }));
            }
            folded += whole.iter().filter(|fold| **fold != Fold::Keep).count();
            assert_eq!(by_cluster, whole, "{rows:?}");
        }
        assert!(folded > 200, "only {folded} rows folded");
    }
}
