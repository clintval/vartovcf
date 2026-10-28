//! What VarDict counted at each position, in VarDict's own coordinates, so a record can carry the
//! depth at its POS and a position without a reference-only record can be filled. Each row's
//! counts go into the sort under the position VarDict counted it at, and the writer walks them in
//! position order.
//!
//! VarDict counts each row's reads at one position: the first deleted base of a deletion it writes
//! with the base before, the anchor of an insertion followed by matched bases, which it writes from
//! the base after, and the first base of every other row. A deletion or MNV also covers each base
//! after that one its REF spans, where its reads have no base VarDict counts.

use std::cmp::Reverse;
use std::collections::BinaryHeap;

use crate::fill::Fill;
use crate::record::{TumorOnlyVariant, variant_class};

/// The most samples a row holds.
pub const SAMPLES: usize = 2;

/// One value per sample, the tumor (or only) sample first.
pub type PerSample = [i32; SAMPLES];

/// A row's insertion or deletion as VarDict's genotype column spells it after the slash: `-` and a
/// deletion's length, `+` and the length of an insertion's description, or `None` for any other
/// allele or a column not in VarDict's notation.
fn vardict_indel<'a>(call: &TumorOnlyVariant<'a>) -> Option<&'a str> {
    let (_, allele) = call.gt.split_once('/')?;
    (allele.starts_with('-') || allele.starts_with('+')).then_some(allele)
}

/// The position VarDict counted a row's reads at. VarDict writes a plain deletion, spelled `-` and a
/// length, from the base before its first deleted base, where it counted it, and an insertion
/// followed by matched bases, spelled `+` and its description's length, which is ALT's plus one,
/// from the base after its anchor; every other row starts where it was counted. A row whose genotype
/// column is not in VarDict's notation is read by its alleles' shape.
pub fn vardict_position(call: &TumorOnlyVariant) -> u64 {
    let digits = |length: &str| length.parse::<usize>().ok();
    match vardict_indel(call) {
        Some(deletion) if deletion.starts_with('-') && digits(&deletion[1..]).is_some() => {
            call.start + 1
        }
        Some(insertion)
            if insertion.starts_with('+')
                && !call.alt_allele.starts_with('<')
                && digits(&insertion[1..]) == Some(call.alt_allele.len() + 1) =>
        {
            call.start - 1
        }
        Some(_) => call.start,
        None => {
            let (ref_allele, alt_allele) = (call.ref_allele.as_bytes(), call.alt_allele.as_bytes());
            let anchored = ref_allele.len() > alt_allele.len()
                && alt_allele.len() == 1
                && ref_allele[0].eq_ignore_ascii_case(&alt_allele[0]);
            if anchored { call.start + 1 } else { call.start }
        }
    }
}

/// Whether VarDict took a row's reads from the REF reads of the base it counted the row at, as it
/// does for an insertion's.
fn inserts(call: &TumorOnlyVariant) -> bool {
    match vardict_indel(call) {
        Some(allele) => allele.starts_with('+') && !call.alt_allele.starts_with('<'),
        None => variant_class(call.ref_allele, call.alt_allele) == Some("Insertion"),
    }
}

/// What VarDict counted in one row, at the position it counted the row at, kept in the sort until
/// the writer reaches that position.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Counts {
    /// Whether VarDict labels the row Complex, as it overcounts such a row's REF reads.
    pub complex: bool,
    /// The last base past the position that a deletion's or MNV's REF covers, or 0 for any other
    /// row.
    pub span_to: u64,
    /// The forward REF reads.
    pub ref_forward: PerSample,
    /// The reverse REF reads.
    pub ref_reverse: PerSample,
    /// The depth.
    pub depth: PerSample,
    /// The ALT reads.
    pub alt_reads: PerSample,
    /// The forward ALT reads of an insertion, which VarDict takes from REF's, or 0 for any other row.
    pub inserted_forward: PerSample,
    /// The reverse ALT reads of an insertion, or 0 for any other row.
    pub inserted_reverse: PerSample,
}

impl Counts {
    /// The position VarDict counted a row at, and what it counted there, from the row's calls, one
    /// per sample with the tumor (or only) sample's first.
    pub fn of(row: &[TumorOnlyVariant]) -> (u64, Counts) {
        let call = &row[0];
        let position = vardict_position(call);
        let each = |value: fn(&TumorOnlyVariant) -> i32| -> PerSample {
            let mut values = [0; SAMPLES];
            for (slot, call) in values.iter_mut().zip(row.iter()) {
                *slot = value(call);
            }
            values
        };
        let span_to = match call.alt_allele {
            "<DEL>" => call.end,
            alt if !alt.starts_with('<') => call.start + call.ref_allele.len() as u64 - 1,
            _ => 0,
        };
        let inserted = inserts(call);
        let zero = [0; SAMPLES];
        let counts = Counts {
            complex: call.variant_type == "Complex",
            span_to: if span_to > position { span_to } else { 0 },
            ref_forward: each(|call| call.ref_forward),
            ref_reverse: each(|call| call.ref_reverse),
            depth: each(|call| call.depth),
            alt_reads: each(|call| call.alt_depth),
            inserted_forward: if inserted {
                each(|call| call.alt_forward)
            } else {
                zero
            },
            inserted_reverse: if inserted {
                each(|call| call.alt_reverse)
            } else {
                zero
            },
        };
        (position, counts)
    }

    /// Appends the counts as bytes, which [`Counts::decode`] reads back.
    pub fn encode(&self, out: &mut Vec<u8>) {
        out.push(u8::from(self.complex));
        out.extend_from_slice(&self.span_to.to_le_bytes());
        for values in self.values() {
            for value in values {
                out.extend_from_slice(&value.to_le_bytes());
            }
        }
    }

    /// Reads back counts that [`Counts::encode`] wrote.
    pub fn decode(bytes: &[u8]) -> Counts {
        let mut counts = Counts {
            complex: bytes[0] == 1,
            span_to: u64::from_le_bytes(bytes[1..9].try_into().unwrap()),
            ..Counts::default()
        };
        let (words, _) = bytes[9..].as_chunks::<4>();
        let mut words = words.iter().map(|word| i32::from_le_bytes(*word));
        for values in counts.values_mut() {
            for value in values.iter_mut() {
                *value = words.next().expect("Counts are encoded whole!");
            }
        }
        counts
    }

    fn values(&self) -> [PerSample; 6] {
        [
            self.ref_forward,
            self.ref_reverse,
            self.depth,
            self.alt_reads,
            self.inserted_forward,
            self.inserted_reverse,
        ]
    }

    fn values_mut(&mut self) -> [&mut PerSample; 6] {
        [
            &mut self.ref_forward,
            &mut self.ref_reverse,
            &mut self.depth,
            &mut self.alt_reads,
            &mut self.inserted_forward,
            &mut self.inserted_reverse,
        ]
    }
}

/// The rows VarDict counted at one position, combined.
#[derive(Clone, Copy)]
struct Counted {
    /// Whether the row whose REF reads stand for the position's is labeled Complex.
    complex: bool,
    /// That row's forward REF reads.
    ref_forward: PerSample,
    /// That row's reverse REF reads.
    ref_reverse: PerSample,
    /// That row's depth.
    depth: PerSample,
    /// The most depth of any row here.
    max_depth: PerSample,
    /// The ALT reads of every row here.
    alt_reads: PerSample,
    /// The forward ALT reads of the insertions here.
    inserted_forward: PerSample,
    /// The reverse ALT reads of the insertions here.
    inserted_reverse: PerSample,
}

impl Counted {
    fn new(counts: &Counts) -> Self {
        Counted {
            complex: counts.complex,
            ref_forward: counts.ref_forward,
            ref_reverse: counts.ref_reverse,
            depth: counts.depth,
            max_depth: counts.depth,
            alt_reads: counts.alt_reads,
            inserted_forward: counts.inserted_forward,
            inserted_reverse: counts.inserted_reverse,
        }
    }

    /// Adds another row counted at the same position: its reads join the position's, and its REF
    /// reads stand for the position's when it is not Complex and the chosen row is, or when it has
    /// more depth in the tumor (or only) sample with neither or both Complex.
    fn add(&mut self, other: &Counts) {
        if (!other.complex, other.depth[0]) > (!self.complex, self.depth[0]) {
            self.complex = other.complex;
            self.ref_forward = other.ref_forward;
            self.ref_reverse = other.ref_reverse;
            self.depth = other.depth;
        }
        for sample in 0..SAMPLES {
            self.max_depth[sample] = self.max_depth[sample].max(other.depth[sample]);
            self.alt_reads[sample] += other.alt_reads[sample];
            self.inserted_forward[sample] += other.inserted_forward[sample];
            self.inserted_reverse[sample] += other.inserted_reverse[sample];
        }
    }
}

/// What is known at one position.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Site {
    /// How a reference-only record there is filled.
    pub fill: Fill,
    /// The forward REF reads of the row VarDict counted there whose REF reads stand for the
    /// position's, with those of each insertion there, or 0 without a row.
    pub ref_forward: PerSample,
    /// The reverse REF reads, counted as the forward ones are.
    pub ref_reverse: PerSample,
    /// The depth of that row, or of the deletion or MNV spanning the position with the most in the
    /// tumor (or only) sample, or 0 without either.
    pub depth: PerSample,
    /// The reads VarDict counted at the position as REF or as any allele, an insertion's once, with
    /// the reads of each deletion or MNV spanning it, no more than the most depth of a row counted
    /// there.
    pub site_dp: PerSample,
}

impl Site {
    /// What is known at a position VarDict counted no row at and no deletion or MNV spans.
    pub const NO_READS: Site = Site {
        fill: Fill::NoReads,
        ref_forward: [0; SAMPLES],
        ref_reverse: [0; SAMPLES],
        depth: [0; SAMPLES],
        site_dp: [0; SAMPLES],
    };
}

/// A walk over one contig's positions in increasing order, fed what VarDict counted at each, which
/// holds only the rows of the position it is at and the deletions and MNVs spanning it.
#[derive(Default)]
pub struct Walk {
    position: u64,
    counted: Option<Counted>,
    pending: Vec<(u64, PerSample, PerSample)>,
    by_end: BinaryHeap<Reverse<(u64, PerSample)>>,
    by_depth: BinaryHeap<(i32, Reverse<u64>, PerSample)>,
    reads: PerSample,
}

impl Walk {
    /// Moves to a position, never below the one it is at; the deletions and MNVs counted at the
    /// position it leaves then span the bases after it.
    pub fn at(&mut self, position: u64) {
        if position <= self.position {
            return;
        }
        for (to, depth, reads) in self.pending.drain(..) {
            if to >= position {
                self.by_end.push(Reverse((to, reads)));
                self.by_depth.push((depth[0], Reverse(to), depth));
                for (total, reads) in self.reads.iter_mut().zip(reads) {
                    *total += reads;
                }
            }
        }
        self.counted = None;
        self.position = position;
        while let Some(&Reverse((to, reads))) = self.by_end.peek() {
            if to >= position {
                break;
            }
            self.by_end.pop();
            for (total, reads) in self.reads.iter_mut().zip(reads) {
                *total -= reads;
            }
        }
        while self
            .by_depth
            .peek()
            .is_some_and(|&(_, Reverse(to), _)| to < position)
        {
            self.by_depth.pop();
        }
    }

    /// Adds a row VarDict counted at the position the walk is at.
    pub fn add(&mut self, counts: &Counts) {
        match &mut self.counted {
            Some(counted) => counted.add(counts),
            None => self.counted = Some(Counted::new(counts)),
        }
        if counts.span_to > self.position {
            self.pending
                .push((counts.span_to, counts.depth, counts.alt_reads));
        }
    }

    /// What is known at the position the walk is at.
    pub fn site(&self) -> Site {
        match (&self.counted, self.by_depth.peek()) {
            (Some(row), _) => {
                let mut site = Site {
                    fill: Fill::Rows,
                    ref_forward: row.ref_forward,
                    ref_reverse: row.ref_reverse,
                    depth: row.depth,
                    site_dp: [0; SAMPLES],
                };
                for sample in 0..SAMPLES {
                    let reads = row.ref_forward[sample]
                        + row.ref_reverse[sample]
                        + row.alt_reads[sample]
                        + self.reads[sample];
                    site.site_dp[sample] = reads.min(row.max_depth[sample]);
                    site.ref_forward[sample] += row.inserted_forward[sample];
                    site.ref_reverse[sample] += row.inserted_reverse[sample];
                }
                site
            }
            (None, Some(&(_, _, depth))) => Site {
                fill: Fill::Spanned,
                ref_forward: [0; SAMPLES],
                ref_reverse: [0; SAMPLES],
                depth,
                site_dp: self.reads,
            },
            (None, None) => Site::NO_READS,
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::rstest;

    use super::*;

    /// A tumor-only row's call with these alleles, depths and REF and ALT reads, each split evenly
    /// across strands.
    fn call(
        start: u64,
        ref_allele: &'static str,
        alt_allele: &'static str,
        gt: &'static str,
        depth: i32,
        ref_reads: i32,
        alt_reads: i32,
    ) -> TumorOnlyVariant<'static> {
        TumorOnlyVariant {
            contig: "chr1",
            start,
            end: start + ref_allele.len() as u64 - 1,
            ref_allele,
            alt_allele,
            gt,
            depth,
            ref_forward: ref_reads / 2,
            ref_reverse: ref_reads - ref_reads / 2,
            alt_depth: alt_reads,
            alt_forward: alt_reads / 2,
            alt_reverse: alt_reads - alt_reads / 2,
            variant_type: match variant_class(ref_allele, alt_allele) {
                Some("SNV") => "SNV",
                Some("Insertion") => "Insertion",
                Some("Deletion") => "Deletion",
                Some(_) if alt_allele.starts_with('<') => "DEL",
                Some(_) => "Complex",
                None => "",
            },
            ..TumorOnlyVariant::default()
        }
    }

    /// What is known at each of these positions, walking from the first row's counted position,
    /// or the first of these positions when it is lower.
    fn walk(rows: &[&[TumorOnlyVariant]], positions: std::ops::RangeInclusive<u64>) -> Vec<Site> {
        let mut counts: Vec<(u64, Counts)> = rows.iter().map(|row| Counts::of(row)).collect();
        counts.sort_by_key(|&(position, _)| position);
        let first = counts
            .first()
            .map_or(*positions.start(), |&(position, _)| position);
        let mut walk = Walk::default();
        (first.min(*positions.start())..=*positions.end())
            .map(|position| {
                walk.at(position);
                for (_, row) in counts.iter().filter(|(at, _)| *at == position) {
                    walk.add(row);
                }
                (position, walk.site())
            })
            .filter(|(position, _)| positions.contains(position))
            .map(|(_, site)| site)
            .collect()
    }

    /// The SITE_DP of the tumor (or only) sample of tumor-only rows at each position.
    fn site_dps(rows: &[TumorOnlyVariant], positions: std::ops::RangeInclusive<u64>) -> Vec<i32> {
        let rows: Vec<&[TumorOnlyVariant]> = rows.iter().map(std::slice::from_ref).collect();
        walk(&rows, positions)
            .into_iter()
            .map(|site| site.site_dp[0])
            .collect()
    }

    #[test]
    fn test_counts_read_back_as_they_were_written() {
        let tumor = TumorOnlyVariant {
            alt_forward: 7,
            alt_reverse: 3,
            ..call(5, "C", "CA", "C/+1", 100, 60, 10)
        };
        let normal = call(5, "C", "CA", "C/+1", -1, i32::MAX, 0);
        let (position, counts) = Counts::of(&[tumor, normal]);
        assert_eq!(position, 5);
        assert_eq!(counts.inserted_forward, [7, 0]);
        let mut bytes = Vec::new();
        counts.encode(&mut bytes);
        assert_eq!(Counts::decode(&bytes), counts);
    }

    #[test]
    fn test_vardict_position_is_the_first_deleted_base_of_an_anchored_deletion() {
        assert_eq!(vardict_position(&call(10, "CA", "C", "", 9, 0, 0)), 11);
        assert_eq!(vardict_position(&call(10, "C", "CA", "", 9, 0, 0)), 10);
        assert_eq!(vardict_position(&call(10, "C", "T", "", 9, 0, 0)), 10);
        assert_eq!(vardict_position(&call(10, "CA", "TG", "", 9, 0, 0)), 10);
        assert_eq!(vardict_position(&call(10, "CAG", "T", "", 9, 0, 0)), 10);
        assert_eq!(vardict_position(&call(10, "C", "C", "", 9, 9, 0)), 10);
        assert_eq!(vardict_position(&call(10, "C", "<DEL>", "", 9, 0, 0)), 10);
    }

    #[rstest]
    #[case("CA", "C", "A/-1", 11)]
    #[case("C", "<DEL>", "C/-1000", 11)]
    #[case("TCG", "T", "T/-2T", 10)]
    #[case("TCG", "T", "T/-2mT", 10)]
    #[case("TCG", "T", "-2T/-2T", 10)]
    #[case("C", "CA", "C/+1", 10)]
    #[case("G", "AT", "C/+3", 9)]
    #[case("CC", "AATG", "C/+5", 9)]
    #[case("C", "<DUP>", "C/+160", 10)]
    fn test_vardict_position_reads_an_indel_from_vardicts_allele_notation(
        #[case] ref_allele: &'static str,
        #[case] alt_allele: &'static str,
        #[case] gt: &'static str,
        #[case] expected: u64,
    ) {
        let row = call(10, ref_allele, alt_allele, gt, 9, 0, 0);
        assert_eq!(vardict_position(&row), expected);
    }

    #[test]
    fn test_site_dp_counts_ref_and_every_allele_at_a_position() {
        let rows = [
            call(5, "C", "A", "C/A", 100, 70, 20),
            call(5, "C", "T", "C/T", 100, 70, 10),
        ];
        assert_eq!(site_dps(&rows, 4..=6), vec![0, 100, 0]);
    }

    #[test]
    fn test_site_dp_counts_an_insertions_reads_once() {
        let rows = [
            call(5, "C", "CA", "C/+1", 100, 60, 40),
            call(5, "C", "T", "C/T", 100, 60, 0),
        ];
        assert_eq!(site_dps(&rows, 5..=5), vec![100]);
    }

    #[test]
    fn test_site_dp_counts_the_reads_of_an_mnv_or_deletion_at_each_base_it_spans() {
        let rows = [
            call(5, "G", "G", "G/G", 300, 300, 0),
            call(5, "GCAT", "G", "G/-3", 300, 200, 100),
            call(9, "AC", "GT", "AC/GT", 300, 150, 150),
            call(10, "C", "C", "C/C", 300, 150, 0),
        ];
        assert_eq!(
            site_dps(&rows, 5..=11),
            vec![300, 300, 100, 100, 300, 300, 0]
        );
    }

    #[test]
    fn test_site_dp_is_no_more_than_the_most_depth_counted_at_the_position() {
        let rows = [call(5, "CG", "AT", "CG/AT", 90, 80, 30)];
        assert_eq!(site_dps(&rows, 5..=6), vec![90, 30]);
    }

    #[test]
    fn test_site_dp_spans_a_symbolic_deletion_to_its_end() {
        let rows = [TumorOnlyVariant {
            end: 15,
            ..call(10, "C", "<DEL>", "C/-5", 200, 0, 200)
        }];
        assert_eq!(
            site_dps(&rows, 10..=16),
            vec![0, 200, 200, 200, 200, 200, 0]
        );
    }

    #[test]
    fn test_site_dp_is_kept_for_each_sample() {
        let tumor = call(5, "C", "A", "C/A", 100, 70, 30);
        let normal = call(5, "C", "A", "C/A", 80, 79, 1);
        let rows = [[tumor, normal]];
        let rows: Vec<&[TumorOnlyVariant]> = rows.iter().map(|row| &row[..]).collect();
        assert_eq!(walk(&rows, 5..=5)[0].site_dp, [100, 80]);
    }

    /// What is known at each position, swept in order, as the tumor (or only) sample's fill, REF
    /// reads, depth and SITE_DP.
    fn sites_at(
        rows: &[TumorOnlyVariant],
        positions: std::ops::RangeInclusive<u64>,
    ) -> Vec<(Fill, i32, i32, i32)> {
        let rows: Vec<&[TumorOnlyVariant]> = rows.iter().map(std::slice::from_ref).collect();
        walk(&rows, positions)
            .into_iter()
            .map(|site| {
                let ref_reads = site.ref_forward[0] + site.ref_reverse[0];
                (site.fill, ref_reads, site.depth[0], site.site_dp[0])
            })
            .collect()
    }

    #[test]
    fn test_a_sites_ref_reads_are_a_row_not_labeled_complex_with_the_most_depth() {
        let rows = [
            call(5, "CG", "AT", "CG/AT", 90, 80, 10),
            call(5, "C", "A", "C/A", 60, 58, 1),
            call(5, "C", "T", "C/T", 61, 59, 1),
        ];
        assert_eq!(sites_at(&rows, 5..=5), vec![(Fill::Rows, 59, 61, 71)]);
        let complex = [call(5, "CG", "AT", "CG/AT", 90, 80, 10)];
        assert_eq!(
            sites_at(&complex, 5..=6),
            vec![(Fill::Rows, 80, 90, 90), (Fill::Spanned, 0, 90, 10)]
        );
    }

    #[test]
    fn test_a_sites_ref_reads_count_the_reads_of_each_insertion_there() {
        let rows = [
            call(5, "G", "GG", "G/+1", 105, 80, 20),
            call(5, "G", "GGG", "G/+2", 105, 80, 5),
        ];
        assert_eq!(sites_at(&rows, 5..=5), vec![(Fill::Rows, 105, 105, 105)]);
    }

    #[test]
    fn test_an_insertion_before_matched_bases_is_counted_at_its_anchor() {
        let rows = [call(6, "G", "AT", "C/+3", 50, 40, 10)];
        assert_eq!(
            sites_at(&rows, 5..=7),
            vec![
                (Fill::Rows, 50, 50, 50),
                (Fill::Spanned, 0, 50, 10),
                (Fill::NoReads, 0, 0, 0),
            ]
        );
    }

    #[test]
    fn test_a_spanned_site_takes_the_depth_of_the_span_with_the_most() {
        let rows = [
            call(10, "CGGAT", "C", "C/-4", 50, 0, 50),
            call(11, "GG", "G", "G/-1", 80, 0, 80),
            call(13, "AT", "A", "A/-1", 60, 0, 60),
        ];
        assert_eq!(
            sites_at(&rows, 13..=15),
            vec![
                (Fill::Spanned, 0, 50, 50),
                (Fill::Rows, 0, 60, 60),
                (Fill::NoReads, 0, 0, 0),
            ]
        );
    }

    #[test]
    fn test_a_site_is_filled_by_what_vardict_counted_there() {
        let rows = [
            call(5, "GCAT", "G", "G/-3", 300, 200, 100),
            call(6, "C", "T", "C/T", 301, 280, 20),
        ];
        let fills: Vec<(Fill, i32, i32)> = sites_at(&rows, 5..=9)
            .into_iter()
            .map(|(fill, ref_reads, depth, _)| (fill, ref_reads, depth))
            .collect();
        assert_eq!(
            fills,
            vec![
                (Fill::NoReads, 0, 0),
                (Fill::Rows, 280, 301),
                (Fill::Spanned, 0, 300),
                (Fill::Spanned, 0, 300),
                (Fill::NoReads, 0, 0),
            ]
        );
    }
}
