//! Folds the partial rows VarDict writes for a multi-nucleotide variant (MNV) back into it.
//!
//! VarDict grows an MNV through a read only while the read goes on with well-called bases, so a read
//! that ends inside the MNV, has a no-call (N) there or a base below VarDict's `-q` plus 5 adds to a
//! shorter allele instead: an SNV or a shorter MNV carrying the MNV's ALT bases over its own span.
//! VarDict folds some of these back itself, the partials ending at the MNV's last base and the ones
//! at its first base whose reads end, on average, within the partial, and writes the rest as rows
//! of their own.
//!
//! Partial rows are grouped with their MNV by the directional adjacency rule of UMI-tools (Smith,
//! Heger and Sudbery 2017): a parent with `n` ALT reads can claim a sub-event with at most
//! `(n + 1) / 2` of them.

use rust_htslib::bcf::record::Numeric;

use crate::fisher::fisher_exact_two_sided;
use crate::record::TumorOnlyVariant;

/// How a row is written once the partial rows of its region are folded.
#[derive(Clone, Debug, PartialEq)]
pub enum Fold<'a> {
    /// The row is written as VarDict wrote it.
    Keep,
    /// The row is a partial of the row at this index and is not written.
    Into(usize),
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
/// covers as the parent does; it can claim the child when, in every sample, it has at least
/// `2 * child - 1` ALT reads. Rows are visited from the most to the least supported, and a row no
/// earlier root can claim is a root itself. A row that several roots can claim goes to the one that
/// needs the fewest breaks to leave it, one for each end of the child inside the root's span, since
/// every break is a rare event in a read; then to the root with the most ALT reads. A child starting
/// at the root's first base is folded only while the root's ALT reads and those of the children it
/// took there fit in that base's DP, as VarDict counts each read's base there once. A row marked
/// `written`, already written for an earlier region, can be a root but is never claimed.
pub fn fold_partials<'a>(rows: &[&[TumorOnlyVariant<'a>]], written: &[bool]) -> Vec<Fold<'a>> {
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
            row,
        )
    });
    let mut roots: Vec<Root> = Vec::new();
    for &node in &nodes {
        let mut candidates: Vec<(usize, usize)> = roots
            .iter()
            .enumerate()
            .filter(|(_, root)| !written[node] && can_claim(rows[root.row], rows[node]))
            .map(|(rank, root)| (breaks(&rows[root.row][0], &rows[node][0]), rank))
            .collect();
        if candidates.is_empty() {
            if rows[node][0].ref_allele.len() > 1 {
                roots.push(Root::new(node, rows[node]));
            }
            continue;
        }
        candidates.sort_unstable();
        if let Some(&(_, rank)) = candidates
            .iter()
            .find(|(_, rank)| roots[*rank].admits(rows, node))
        {
            roots[rank].take(rows, node);
        } else if rows[node][0].ref_allele.len() > 1 {
            roots.push(Root::new(node, rows[node]));
        }
    }
    for root in roots.iter_mut().filter(|root| !root.children.is_empty()) {
        root.children.sort_unstable();
        for &child in &root.children {
            folds[child] = Fold::Into(root.row);
        }
        let children: Vec<&[TumorOnlyVariant<'a>]> =
            root.children.iter().map(|&child| rows[child]).collect();
        folds[root.row] = Fold::Merged(merge(rows[root.row], &children));
    }
    folds
}

/// An MNV row that partial rows fold into, with the ALT reads it holds at its first base.
struct Root {
    row: usize,
    children: Vec<usize>,
    first_base_reads: Vec<i64>,
}

impl Root {
    /// A root without children yet.
    fn new(row: usize, calls: &[TumorOnlyVariant]) -> Self {
        Root {
            row,
            children: Vec::new(),
            first_base_reads: calls.iter().map(|call| call.alt_depth as i64).collect(),
        }
    }

    /// Whether a child's reads fit beside this root's: always for a child starting after the root's
    /// first base, and otherwise only while every sample's ALT reads at that base stay within its DP.
    fn admits(&self, rows: &[&[TumorOnlyVariant]], child: usize) -> bool {
        rows[child][0].start != rows[self.row][0].start
            || rows[self.row]
                .iter()
                .zip(rows[child])
                .zip(&self.first_base_reads)
                .all(|((root, call), &reads)| {
                    call.alt_depth == 0 || reads + call.alt_depth as i64 <= root.depth as i64
                })
    }

    /// Takes a child into this root.
    fn take(&mut self, rows: &[&[TumorOnlyVariant]], child: usize) {
        if rows[child][0].start == rows[self.row][0].start {
            for (reads, call) in self.first_base_reads.iter_mut().zip(rows[child]) {
                *reads += call.alt_depth as i64;
            }
        }
        self.children.push(child);
    }
}

/// Whether a parent row can claim a child row: the child is its sub-event and, in every sample, the
/// parent has at least `2 * child - 1` ALT reads.
fn can_claim(parent: &[TumorOnlyVariant], child: &[TumorOnlyVariant]) -> bool {
    is_sub_event(&parent[0], &child[0])
        && parent
            .iter()
            .zip(child)
            .all(|(p, c)| p.alt_depth >= 2 * c.alt_depth - 1)
}

/// The breaks a read of `parent` needs to leave only `child`: one for each end of the child that
/// falls inside the parent's span.
fn breaks(parent: &TumorOnlyVariant, child: &TumorOnlyVariant) -> usize {
    usize::from(child.start != parent.start) + usize::from(child.end != parent.end)
}

/// Whether a call is a substitution of one or more bases with a simple ALT allele.
fn is_substitution(call: &TumorOnlyVariant) -> bool {
    call.ref_allele.len() == call.alt_allele.len()
        && call.ref_allele != call.alt_allele
        && call
            .alt_allele
            .bytes()
            .all(|base| base.is_ascii_alphabetic())
}

/// Whether `child` lies within `parent`'s span and carries `parent`'s ALT bases over its own span.
fn is_sub_event(parent: &TumorOnlyVariant, child: &TumorOnlyVariant) -> bool {
    if child.ref_allele.len() >= parent.ref_allele.len() || child.start < parent.start {
        return false;
    }
    let offset = (child.start - parent.start) as usize;
    let span = offset..offset + child.ref_allele.len();
    parent.alt_allele.get(span.clone()) == Some(child.alt_allele)
        && parent.ref_allele.get(span) == Some(child.ref_allele)
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
        strand_bias_p_value: fisher_exact_two_sided(
            parent.ref_forward.max(0) as u64,
            parent.ref_reverse.max(0) as u64,
            alt_forward.max(0) as u64,
            alt_reverse.max(0) as u64,
        ) as f32,
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
        fold_partials(&rows, &vec![false; rows.len()])
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
            call(20, "TTGG", "CTGA", 6, 200),
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
    #[case(9, 5, Some(0))]
    #[case(9, 6, None)]
    #[case(1, 1, Some(0))]
    #[case(2, 2, None)]
    fn test_the_gradient_is_two_n_minus_one(
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

    #[test]
    fn test_chained_partials_fold_into_the_longest_mnv() {
        let folds = fold(&[
            call(10, "GGG", "CAT", 100, 200),
            call(10, "GG", "CA", 10, 200),
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
            (119, 206)
        );
    }

    #[test]
    fn test_a_partial_two_mnvs_share_goes_to_the_better_supported() {
        let folds = fold(&[
            call(10, "GGG", "CTT", 40, 200),
            call(10, "GGG", "CAA", 100, 200),
            call(10, "G", "C", 5, 200),
        ]);
        assert_eq!(targets(&folds), vec![None, None, Some(1)]);
    }

    #[test]
    fn test_a_partial_goes_to_the_mnv_that_needs_the_fewest_breaks() {
        let folds = fold(&[
            call(9, "GGG", "TCA", 100, 200),
            call(10, "GG", "CT", 40, 200),
            call(10, "G", "C", 5, 200),
        ]);
        assert_eq!(targets(&folds), vec![None, None, Some(1)]);
    }

    #[test]
    fn test_a_partial_at_the_first_base_folds_only_while_it_fits_in_dp() {
        let folds = fold(&[
            call(10, "GGG", "CAT", 100, 110),
            call(10, "GG", "CA", 8, 110),
            call(10, "G", "C", 5, 110),
        ]);
        assert_eq!(targets(&folds), vec![None, Some(0), None]);
        let merged = merged(&folds, 0);
        assert!(merged.calls[0].alt_depth <= merged.calls[0].depth);
    }

    #[test]
    fn test_an_mnv_no_root_has_room_for_takes_its_own_partials() {
        let folds = fold(&[
            call(10, "GGG", "CAT", 100, 104),
            call(10, "GG", "CA", 8, 104),
            call(11, "G", "A", 3, 104),
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
    fn test_a_row_written_before_is_not_claimed_but_can_be_a_root() {
        let rows = [
            call(10, "GG", "CA", 100, 200),
            call(10, "G", "C", 5, 200),
            call(11, "G", "A", 4, 205),
        ];
        let calls: Vec<&[TumorOnlyVariant]> = rows.iter().map(std::slice::from_ref).collect();
        let folds = fold_partials(&calls, &[false, true, false]);
        assert_eq!(targets(&folds), vec![None, None, Some(0)]);
        let folds = fold_partials(&calls, &[true, false, false]);
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
            alt_forward: 60,
            alt_reverse: 30,
            mean_position_in_read: 20.0,
            stdev_position_in_read: 0.0,
            base_quality_mean: 90.0,
            mean_mapping_quality: 60.0,
            mean_mismatches_in_reads: 1.0,
            af_adjusted: 0.05,
            ..call(10, "GG", "CA", 90, 200)
        };
        let partial = TumorOnlyVariant {
            alt_forward: 10,
            alt_reverse: 0,
            mean_position_in_read: 2.0,
            stdev_position_in_read: 0.0,
            base_quality_mean: 70.0,
            mean_mapping_quality: 40.0,
            mean_mismatches_in_reads: 3.0,
            ..call(10, "G", "C", 10, 200)
        };
        let folds = fold(&[mnv, partial]);
        let merged = &merged(&folds, 0).calls[0];
        assert_eq!((merged.alt_forward, merged.alt_reverse), (70, 30));
        assert_eq!(merged.mean_position_in_read, 18.2);
        assert_eq!(merged.stdev_position_in_read, 1.0);
        assert_eq!(merged.base_quality_mean, 88.0);
        assert_eq!(merged.mean_mapping_quality, 58.0);
        assert_eq!(merged.mean_mismatches_in_reads, 1.2);
        assert_eq!(merged.high_quality_variant_reads, 100);
        assert_eq!(merged.af_adjusted, 0.1);
        let expected = fisher_exact_two_sided(50, 50, 70, 30) as f32;
        assert_eq!(merged.strand_bias_p_value, expected);
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
        fold_partials(&rows, &vec![false; rows.len()])
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
    fn test_a_tumor_normal_partial_must_follow_the_gradient_in_every_sample() {
        let folds = fold_pairs(&[
            [call(10, "GG", "CA", 100, 200), call(10, "GG", "CA", 0, 100)],
            [call(10, "G", "C", 4, 200), call(10, "G", "C", 3, 100)],
        ]);
        assert_eq!(folds, vec![Fold::Keep, Fold::Keep]);
    }
}
