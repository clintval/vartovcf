//! FILTER labels for VarDict calls.

use crate::reads::ReadEvidence;
use crate::record::TumorOnlyVariant;
use crate::vcf::Numeric;

/// The FILTER label for calls whose ALT reads carry the variant near their read ends.
pub const NEAR_READ_END: &str = "NEAR_READ_END";

/// The FILTER label for calls whose ALT reads have a low mean mapping quality.
pub const LOW_MEAN_MAPQ: &str = "LOW_MEAN_MAPQ";

/// The FILTER label for one-base insertions or deletions in long homopolymers.
pub const HOMOPOLYMER_INDEL: &str = "HOMOPOLYMER_INDEL";

/// The FILTER label for one-unit insertions or deletions in long 2-6 bp tandem repeats.
pub const TANDEM_REPEAT_INDEL: &str = "TANDEM_REPEAT_INDEL";

/// The FILTER label for calls whose ALT reads carry many other mismatches.
pub const HIGH_MEAN_MISMATCHES: &str = "HIGH_MEAN_MISMATCHES";

/// The FILTER label for calls whose ALT reads all place the variant at the same read position.
pub const SAME_READ_POSITION: &str = "SAME_READ_POSITION";

/// The FILTER label for calls whose ALT reads split across strands differently from the REF reads.
pub const STRAND_BIAS: &str = "STRAND_BIAS";

/// The FILTER label for calls with a low allele frequency.
pub const LOW_AF: &str = "LOW_AF";

/// The FILTER label for calls whose ALT reads have a low mean variant base quality.
pub const LOW_QMEAN: &str = "LOW_QMEAN";

/// The FILTER label for calls with a low read depth.
pub const LOW_DP: &str = "LOW_DP";

/// The FILTER label for calls with few high-quality ALT reads.
pub const LOW_HICNT: &str = "LOW_HICNT";

/// The FILTER label for SNVs that VarDict counted at or above the read floor, by no more than the
/// reads it reassigned could add, but whose reads fall short of it when recounted.
pub const FOREIGN_CLIP_SUPPORT: &str = "FOREIGN_CLIP_SUPPORT";

/// The FILTER label for long insertions with few high-quality ALT reads that no read carries
/// whole.
pub const LOW_SUPPORT_LONG_INSERTION: &str = "LOW_SUPPORT_LONG_INSERTION";

/// The rule of the LOW_SUPPORT_LONG_INSERTION label: which insertions it tests, by the bases they
/// add, and the high-quality ALT reads one needs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LowSupportLongInsertion {
    /// The fewest added bases of an insertion the label tests.
    pub min_length: usize,
    /// The number of high-quality ALT reads below which a tested insertion is labelled.
    pub min_hicnt: i32,
}

impl LowSupportLongInsertion {
    /// Whether a call is one the label tests, and so is looked for in the reads when they are
    /// given: an insertion, or a Complex call longer than its REF, adding at least the length,
    /// with its high-quality ALT reads known and under the floor, and with REF reads, since a
    /// site where every read carries the insertion is the sample's own sequence.
    pub fn tests(&self, variant: &TumorOnlyVariant) -> bool {
        let hicnt = variant.hicnt_value();
        variant
            .inserted_length()
            .is_some_and(|length| length >= self.min_length)
            && !hicnt.is_missing()
            && hicnt < self.min_hicnt
            && variant.ref_forward.saturating_add(variant.ref_reverse) > 0
    }
}

/// The thresholds of the FILTER labels to apply; a label without a threshold is not applied.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FilterThresholds {
    /// The minimum mean distance from the variant to the nearer read end.
    pub near_read_end: Option<f32>,
    /// The minimum mean mapping quality of the ALT reads.
    pub low_mean_mapq: Option<f32>,
    /// The minimum homopolymer copies for a one-base indel to be labelled.
    pub homopolymer_indel: Option<f32>,
    /// The allele frequency below which a homopolymer indel is labelled; no limit when not given.
    pub homopolymer_indel_max_af: Option<f32>,
    /// The minimum repeat copies for a one-unit indel in a 2-6 bp tandem repeat to be labelled.
    pub tandem_repeat_indel: Option<f32>,
    /// The allele frequency below which a tandem repeat indel is labelled; no limit when not given.
    pub tandem_repeat_indel_max_af: Option<f32>,
    /// The maximum mean substitution mismatches per ALT read.
    pub high_mean_mismatches: Option<f32>,
    /// The allele frequency below which a call whose ALT reads share one read position is labelled.
    pub same_read_position: Option<f32>,
    /// The strand-bias Fisher p-value below which a call is labelled.
    pub strand_bias: Option<f32>,
    /// The folded strand odds ratio a labelled call must exceed, unless the table has an empty cell.
    pub strand_bias_min_odds_ratio: Option<f32>,
    /// The allele frequency below which a strand-biased call is labelled; no limit when not given.
    pub strand_bias_max_af: Option<f32>,
    /// The minimum allele frequency.
    pub low_af: Option<f32>,
    /// The minimum mean per-read variant base quality of the ALT reads.
    pub low_qmean: Option<f32>,
    /// The minimum read depth.
    pub low_dp: Option<i32>,
    /// The minimum number of high-quality ALT reads.
    pub low_hicnt: Option<i32>,
    /// Whether SNVs are recounted from the reads and held to the LOW_HICNT threshold; applied only
    /// with that threshold.
    pub foreign_clip_support: bool,
    /// The long insertions to label when they have few high-quality ALT reads and no read carries
    /// them whole.
    pub low_support_long_insertion: Option<LowSupportLongInsertion>,
}

impl FilterThresholds {
    /// Whether any FILTER label is applied.
    pub fn any(&self) -> bool {
        self.near_read_end.is_some()
            || self.low_mean_mapq.is_some()
            || self.homopolymer_indel.is_some()
            || self.tandem_repeat_indel.is_some()
            || self.high_mean_mismatches.is_some()
            || self.same_read_position.is_some()
            || self.strand_bias.is_some()
            || self.low_af.is_some()
            || self.low_qmean.is_some()
            || self.low_dp.is_some()
            || self.low_hicnt.is_some()
            || self.low_support_long_insertion.is_some()
    }

    /// The read floor FOREIGN_CLIP_SUPPORT holds recounted SNVs to: LOW_HICNT's, when both apply.
    pub fn foreign_clip_support_floor(&self) -> Option<i32> {
        self.low_hicnt.filter(|_| self.foreign_clip_support)
    }

    /// Whether an SNV call is recounted from the reads for FOREIGN_CLIP_SUPPORT: one whose known
    /// count of high-quality ALT reads is at the floor or above, and would fall under it without
    /// the reads VarDict reassigned to it, the only way soft clips are added to an SNV.
    pub fn recounts(&self, variant: &TumorOnlyVariant) -> bool {
        self.foreign_clip_support_floor().is_some_and(|min_hicnt| {
            let hicnt = variant.hicnt_value();
            variant.is_snv()
                && !hicnt.is_missing()
                && hicnt >= min_hicnt
                && i64::from(hicnt - min_hicnt) < variant.max_reassigned_reads()
        })
    }

    /// Return the ID and description of the FILTER header line of each applied label.
    pub fn header_filters(&self) -> Vec<(&'static str, String)> {
        let mut filters = Vec::new();
        if let Some(min_mean_dist) = self.near_read_end {
            filters.push((
                NEAR_READ_END,
                format!(r#"Mean distance from the variant to the nearer end of the ALT reads' aligned parts (FORMAT MEAN_DIST_TO_READ_END) is below {min_mean_dist}."#),
            ));
        }
        if let Some(min_mean_mapq) = self.low_mean_mapq {
            filters.push((
                LOW_MEAN_MAPQ,
                format!(r#"Mean mapping quality of the ALT reads (FORMAT MEAN_MAPQ) is below {min_mean_mapq}."#),
            ));
        }
        if let Some(min_copies) = self.homopolymer_indel {
            let af_limit = af_limit_clause(self.homopolymer_indel_max_af);
            filters.push((
                HOMOPOLYMER_INDEL,
                format!(r#"Call exactly one base longer or shorter than REF (an insertion, a deletion or a Complex call) in a homopolymer of at least {min_copies} copies (INFO REPEAT_UNIT_COPIES with REPEAT_UNIT_LEN 1){af_limit}, the signature of polymerase slippage."#),
            ));
        }
        if let Some(min_copies) = self.tandem_repeat_indel {
            let af_limit = af_limit_clause(self.tandem_repeat_indel_max_af);
            filters.push((
                TANDEM_REPEAT_INDEL,
                format!(r#"Call exactly one repeat unit longer or shorter than REF (an insertion, a deletion or a Complex call) in a 2-6 bp tandem repeat of at least {min_copies} copies (INFO REPEAT_UNIT_COPIES and REPEAT_UNIT_LEN){af_limit}, the signature of polymerase slippage."#),
            ));
        }
        if let Some(max_mean_mismatches) = self.high_mean_mismatches {
            filters.push((
                HIGH_MEAN_MISMATCHES,
                format!(r#"Mean substitution mismatches per ALT read (FORMAT MEAN_MISMATCHES) is above {max_mean_mismatches}, suggesting misaligned or paralogous reads."#),
            ));
        }
        if let Some(max_af) = self.same_read_position {
            filters.push((
                SAME_READ_POSITION,
                format!(r#"Every one of at least 2 ALT reads places the variant at the same distance from its read end (FORMAT ALT_READ_POS_VARIES 0) at an AF (FORMAT AF) below {max_af}, a sign of one fragment counted repeatedly or an error at a fixed sequencing cycle."#),
            ));
        }
        if let Some(max_p) = self.strand_bias {
            let odds_ratio = self
                .strand_bias_min_odds_ratio
                .map(|min_odds_ratio| format!(", the odds ratio of that table, folded so either direction counts, is above {min_odds_ratio} or the table has an empty cell"))
                .unwrap_or_default();
            let af_limit = self
                .strand_bias_max_af
                .map(|max_af| format!(", and AF (FORMAT AF) is below {max_af}"))
                .unwrap_or_default();
            filters.push((
                STRAND_BIAS,
                format!(r#"The ALT reads split across strands differently from the REF reads: the two-sided Fisher p-value of REF and ALT reads by strand (FORMAT STRAND_BIAS_FISHER_P) is below {max_p}{odds_ratio}{af_limit}."#),
            ));
        }
        if let Some(min_af) = self.low_af {
            filters.push((
                LOW_AF,
                format!(
                    r#"Fraction of reads carrying the ALT allele (FORMAT AF) is below {min_af}."#
                ),
            ));
        }
        if let Some(min_qmean) = self.low_qmean {
            filters.push((
                LOW_QMEAN,
                format!(r#"Mean per-read variant base quality of the ALT reads (FORMAT QMEAN) is below {min_qmean}."#),
            ));
        }
        if let Some(min_dp) = self.low_dp {
            filters.push((
                LOW_DP,
                format!(r#"Read depth (FORMAT DP) is below {min_dp}."#),
            ));
        }
        if let Some(min_hicnt) = self.low_hicnt {
            filters.push((
                LOW_HICNT,
                format!(r#"High-quality ALT reads (FORMAT HICNT) are fewer than {min_hicnt}."#),
            ));
        }
        if let Some(min_hicnt) = self.foreign_clip_support_floor() {
            filters.push((
                FOREIGN_CLIP_SUPPORT,
                format!(r#"SNV with at least {min_hicnt} high-quality ALT reads by VarDict (FORMAT HICNT), a count the reads it reassigned to the allele (FORMAT REALIGNED_FRAC_OF_DP) could have lifted to {min_hicnt}, but fewer than {min_hicnt} in the reads it was run on: reads carrying the ALT base aligned at the site at or above VarDict's -q, once its clipping of mismatches at read ends is replayed, or first in a soft clip at the site whose later bases match the reference. The shortfall is most often a soft clip VarDict added after comparing only its second base to the reference. Every count is left as VarDict wrote it."#),
            ));
        }
        if let Some(rule) = self.low_support_long_insertion {
            filters.push((
                LOW_SUPPORT_LONG_INSERTION,
                format!(r#"Insertion, or Complex call longer than its REF, adding at least {} bases to the reference, with fewer than {} high-quality ALT reads (FORMAT HICNT), at least one read VarDict counted as REF (FORMAT AD, left missing for Complex calls) and, when the reads VarDict ran on are given, no read carrying it whole: for an insertion, every added base in one insertion operation at a mean quality at or above VarDict's -q, at the call's position or one its insertion slides to; for a Complex call, its ALT spelled by the read's bases over its REF and the insertion operations among them. Every count is left as VarDict wrote it."#, rule.min_length, rule.min_hicnt),
            ));
        }
        filters
    }

    /// Return the FILTER header lines of the applied labels.
    pub fn header_lines(&self) -> Vec<String> {
        self.header_filters()
            .into_iter()
            .map(|(id, description)| format!(r#"##FILTER=<ID={id},Description="{description}">"#))
            .collect()
    }

    /// Return the FILTER labels a call fails by its FORMAT values, and by what the reads say of it
    /// when they were asked: none without reads of an ALT allele, and none from a value that is
    /// missing.
    pub fn labels(&self, variant: &TumorOnlyVariant, reads: ReadEvidence) -> Vec<&'static str> {
        let mut labels = Vec::new();
        if !variant.has_alt_reads() {
            return labels;
        }
        if self
            .near_read_end
            .is_some_and(|min_mean_dist| variant.mean_dist_to_read_end_value() < min_mean_dist)
        {
            labels.push(NEAR_READ_END);
        }
        if self
            .low_mean_mapq
            .is_some_and(|min_mean_mapq| variant.mean_mapq_value() < min_mean_mapq)
        {
            labels.push(LOW_MEAN_MAPQ);
        }
        if let (Some(copies), Some(unit_length)) = (
            variant.repeat_unit_copies_value(),
            variant.repeat_unit_len_value(),
        ) {
            let length_change = variant.ref_allele.len().abs_diff(variant.alt_allele.len());
            let is_one_unit_indel = unit_length > 0 && length_change == unit_length as usize;
            let af = variant.af_value();
            let below = |max_af: Option<f32>| max_af.is_none_or(|max_af| af < max_af);
            if is_one_unit_indel
                && unit_length == 1
                && self
                    .homopolymer_indel
                    .is_some_and(|min_copies| copies >= min_copies)
                && below(self.homopolymer_indel_max_af)
            {
                labels.push(HOMOPOLYMER_INDEL);
            }
            if is_one_unit_indel
                && unit_length > 1
                && self
                    .tandem_repeat_indel
                    .is_some_and(|min_copies| copies >= min_copies)
                && below(self.tandem_repeat_indel_max_af)
            {
                labels.push(TANDEM_REPEAT_INDEL);
            }
        }
        if self
            .high_mean_mismatches
            .is_some_and(|max_mean_mismatches| {
                variant.mean_mismatches_value() > max_mean_mismatches
            })
        {
            labels.push(HIGH_MEAN_MISMATCHES);
        }
        if self.same_read_position.is_some_and(|max_af| {
            variant.alt_depth >= 2
                && variant.alt_read_pos_varies_value() == 0
                && variant.af_value() < max_af
        }) {
            labels.push(SAME_READ_POSITION);
        }
        if self.strand_bias.is_some_and(|max_p| {
            self.strand_bias_max_af
                .is_none_or(|max_af| variant.af_value() < max_af)
                && self
                    .strand_bias_min_odds_ratio
                    .is_none_or(|min_odds_ratio| strand_odds_ratio_exceeds(variant, min_odds_ratio))
                && variant.strand_bias_fisher_p_value() < max_p
        }) {
            labels.push(STRAND_BIAS);
        }
        if self
            .low_af
            .is_some_and(|min_af| variant.af_value() < min_af)
        {
            labels.push(LOW_AF);
        }
        if self
            .low_qmean
            .is_some_and(|min_qmean| variant.qmean_value() < min_qmean)
        {
            labels.push(LOW_QMEAN);
        }
        if self.low_dp.is_some_and(|min_dp| variant.depth < min_dp) {
            labels.push(LOW_DP);
        }
        if self.low_hicnt.is_some_and(|min_hicnt| {
            let hicnt = variant.hicnt_value();
            !hicnt.is_missing() && hicnt < min_hicnt
        }) {
            labels.push(LOW_HICNT);
        }
        if let (Some(min_hicnt), Some(recount)) = (self.foreign_clip_support_floor(), reads.recount)
            && self.recounts(variant)
            && recount.supporting() < min_hicnt
        {
            labels.push(FOREIGN_CLIP_SUPPORT);
        }
        if self
            .low_support_long_insertion
            .is_some_and(|rule| rule.tests(variant))
            && reads.spanning.is_none_or(|spanning| spanning == 0)
        {
            labels.push(LOW_SUPPORT_LONG_INSERTION);
        }
        labels
    }
}

/// Whether the REF and ALT strand table has an empty cell or a folded odds ratio above the minimum.
fn strand_odds_ratio_exceeds(variant: &TumorOnlyVariant, min_odds_ratio: f32) -> bool {
    let cells = [
        variant.ref_forward,
        variant.ref_reverse,
        variant.alt_forward,
        variant.alt_reverse,
    ];
    if cells.contains(&0) {
        return true;
    }
    let odds_ratio = (variant.ref_forward as f64 * variant.alt_reverse as f64)
        / (variant.ref_reverse as f64 * variant.alt_forward as f64);
    odds_ratio.max(1.0 / odds_ratio) > min_odds_ratio as f64
}

/// Return the AF clause of a repeat label's description, empty when it has no AF limit.
fn af_limit_clause(max_af: Option<f32>) -> String {
    max_af
        .map(|max_af| format!(" at an AF (FORMAT AF) below {max_af}"))
        .unwrap_or_default()
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::rstest;

    use super::*;
    use crate::reads::Recount;

    fn variant(
        ref_allele: &'static str,
        alt_allele: &'static str,
        mean_dist: f32,
    ) -> TumorOnlyVariant<'static> {
        TumorOnlyVariant {
            ref_allele,
            alt_allele,
            mean_position_in_read: mean_dist,
            mean_mapping_quality: 60.0,
            alt_depth: 1,
            depth: 100,
            ..Default::default()
        }
    }

    #[test]
    fn test_a_call_without_alt_reads_gets_no_labels() {
        let filters = FilterThresholds {
            near_read_end: Some(8.0),
            low_mean_mapq: Some(10.0),
            high_mean_mismatches: Some(5.25),
            same_read_position: Some(0.35),
            strand_bias: Some(0.01),
            low_af: Some(0.01),
            low_qmean: Some(22.5),
            low_dp: Some(50),
            ..Default::default()
        };
        let no_alt_reads = TumorOnlyVariant {
            alt_depth: 0,
            depth: 30,
            mean_mapping_quality: 0.0,
            base_quality_mean: 0.0,
            mean_mismatches_in_reads: 9.0,
            ref_forward: 500,
            ref_reverse: 500,
            alt_reverse: 10,
            ..variant("G", "A", 0.0)
        };
        assert!(
            filters
                .labels(&no_alt_reads, ReadEvidence::default())
                .is_empty()
        );
    }

    fn mapq_variant(
        ref_allele: &'static str,
        alt_allele: &'static str,
        mean_mapq: f32,
    ) -> TumorOnlyVariant<'static> {
        TumorOnlyVariant {
            mean_mapping_quality: mean_mapq,
            ..variant(ref_allele, alt_allele, 30.0)
        }
    }

    #[test]
    fn test_no_label_is_applied_by_default() {
        let filters = FilterThresholds::default();
        assert!(!filters.any());
        assert!(filters.header_lines().is_empty());
        assert!(
            filters
                .labels(&variant("G", "A", 1.0), ReadEvidence::default())
                .is_empty()
        );
    }

    #[rstest]
    #[case("G", "A", 3.0, vec![NEAR_READ_END])]
    #[case("G", "A", 8.0, vec![])]
    #[case("G", "G", 3.0, vec![])]
    fn test_near_read_end(
        #[case] ref_allele: &'static str,
        #[case] alt_allele: &'static str,
        #[case] mean_dist: f32,
        #[case] expected: Vec<&str>,
    ) {
        let filters = FilterThresholds {
            near_read_end: Some(8.0),
            ..Default::default()
        };
        assert!(filters.any());
        assert_eq!(
            filters.labels(
                &variant(ref_allele, alt_allele, mean_dist),
                ReadEvidence::default()
            ),
            expected
        );
    }

    #[test]
    fn test_near_read_end_header_line_states_the_threshold() {
        let filters = FilterThresholds {
            near_read_end: Some(7.5),
            ..Default::default()
        };
        assert_eq!(
            filters.header_lines(),
            vec![
                r#"##FILTER=<ID=NEAR_READ_END,Description="Mean distance from the variant to the nearer end of the ALT reads' aligned parts (FORMAT MEAN_DIST_TO_READ_END) is below 7.5.">"#
            ]
        );
    }

    #[rstest]
    #[case("G", "A", 5.0, vec![LOW_MEAN_MAPQ])]
    #[case("G", "A", 10.0, vec![])]
    #[case("G", "G", 5.0, vec![])]
    fn test_low_mean_mapq(
        #[case] ref_allele: &'static str,
        #[case] alt_allele: &'static str,
        #[case] mean_mapq: f32,
        #[case] expected: Vec<&str>,
    ) {
        let filters = FilterThresholds {
            low_mean_mapq: Some(10.0),
            ..Default::default()
        };
        assert!(filters.any());
        assert_eq!(
            filters.labels(
                &mapq_variant(ref_allele, alt_allele, mean_mapq),
                ReadEvidence::default()
            ),
            expected
        );
    }

    #[test]
    fn test_low_mean_mapq_header_line_states_the_threshold() {
        let filters = FilterThresholds {
            low_mean_mapq: Some(10.0),
            ..Default::default()
        };
        assert_eq!(
            filters.header_lines(),
            vec![
                r#"##FILTER=<ID=LOW_MEAN_MAPQ,Description="Mean mapping quality of the ALT reads (FORMAT MEAN_MAPQ) is below 10.">"#
            ]
        );
    }

    fn repeat_variant(
        ref_allele: &'static str,
        alt_allele: &'static str,
        unit_length: i32,
        copies: f32,
        alt_depth: i32,
    ) -> TumorOnlyVariant<'static> {
        TumorOnlyVariant {
            microsatellite: copies,
            microsatellite_length: unit_length,
            alt_depth,
            depth: 100,
            ..variant(ref_allele, alt_allele, 30.0)
        }
    }

    #[rstest]
    #[case("GA", "G", 1, 13.0, 10, vec![HOMOPOLYMER_INDEL])]
    #[case("G", "GA", 1, 13.0, 10, vec![HOMOPOLYMER_INDEL])]
    #[case("GA", "G", 1, 12.0, 10, vec![])]
    #[case("GA", "G", 1, 13.0, 28, vec![])]
    #[case("GAA", "G", 1, 13.0, 10, vec![])]
    #[case("G", "A", 1, 20.0, 10, vec![])]
    #[case("GCA", "G", 2, 13.0, 10, vec![])]
    #[case("G", "G", 1, 13.0, 10, vec![])]
    fn test_homopolymer_indel(
        #[case] ref_allele: &'static str,
        #[case] alt_allele: &'static str,
        #[case] unit_length: i32,
        #[case] copies: f32,
        #[case] alt_depth: i32,
        #[case] expected: Vec<&str>,
    ) {
        let filters = FilterThresholds {
            homopolymer_indel: Some(13.0),
            homopolymer_indel_max_af: Some(0.275),
            ..Default::default()
        };
        assert!(filters.any());
        let variant = repeat_variant(ref_allele, alt_allele, unit_length, copies, alt_depth);
        assert_eq!(filters.labels(&variant, ReadEvidence::default()), expected);
    }

    #[rstest]
    #[case("GCA", "G", 2, 13.0, 10, vec![TANDEM_REPEAT_INDEL])]
    #[case("G", "GCAG", 3, 13.0, 10, vec![TANDEM_REPEAT_INDEL])]
    #[case("GCA", "G", 2, 12.5, 10, vec![])]
    #[case("GCA", "G", 2, 13.0, 20, vec![])]
    #[case("GCACA", "G", 2, 13.0, 10, vec![])]
    #[case("GA", "G", 1, 13.0, 10, vec![])]
    fn test_tandem_repeat_indel(
        #[case] ref_allele: &'static str,
        #[case] alt_allele: &'static str,
        #[case] unit_length: i32,
        #[case] copies: f32,
        #[case] alt_depth: i32,
        #[case] expected: Vec<&str>,
    ) {
        let filters = FilterThresholds {
            tandem_repeat_indel: Some(13.0),
            tandem_repeat_indel_max_af: Some(0.2),
            ..Default::default()
        };
        assert!(filters.any());
        let variant = repeat_variant(ref_allele, alt_allele, unit_length, copies, alt_depth);
        assert_eq!(filters.labels(&variant, ReadEvidence::default()), expected);
    }

    #[test]
    fn test_repeat_indels_have_no_af_limit_without_a_max_af() {
        let filters = FilterThresholds {
            homopolymer_indel: Some(13.0),
            tandem_repeat_indel: Some(13.0),
            ..Default::default()
        };
        assert_eq!(
            filters.labels(
                &repeat_variant("GA", "G", 1, 13.0, 90),
                ReadEvidence::default()
            ),
            vec![HOMOPOLYMER_INDEL]
        );
        assert_eq!(
            filters.labels(
                &repeat_variant("GCA", "G", 2, 13.0, 90),
                ReadEvidence::default()
            ),
            vec![TANDEM_REPEAT_INDEL]
        );
    }

    #[test]
    fn test_repeat_indel_header_lines_state_the_thresholds() {
        let filters = FilterThresholds {
            homopolymer_indel: Some(13.0),
            homopolymer_indel_max_af: Some(0.275),
            tandem_repeat_indel: Some(13.0),
            ..Default::default()
        };
        assert_eq!(
            filters.header_lines(),
            vec![
                r#"##FILTER=<ID=HOMOPOLYMER_INDEL,Description="Call exactly one base longer or shorter than REF (an insertion, a deletion or a Complex call) in a homopolymer of at least 13 copies (INFO REPEAT_UNIT_COPIES with REPEAT_UNIT_LEN 1) at an AF (FORMAT AF) below 0.275, the signature of polymerase slippage.">"#,
                r#"##FILTER=<ID=TANDEM_REPEAT_INDEL,Description="Call exactly one repeat unit longer or shorter than REF (an insertion, a deletion or a Complex call) in a 2-6 bp tandem repeat of at least 13 copies (INFO REPEAT_UNIT_COPIES and REPEAT_UNIT_LEN), the signature of polymerase slippage.">"#,
            ]
        );
    }

    #[rstest]
    #[case("G", "A", 5.5, vec![HIGH_MEAN_MISMATCHES])]
    #[case("G", "A", 5.25, vec![])]
    #[case("G", "G", 5.5, vec![])]
    fn test_high_mean_mismatches(
        #[case] ref_allele: &'static str,
        #[case] alt_allele: &'static str,
        #[case] mean_mismatches: f32,
        #[case] expected: Vec<&str>,
    ) {
        let filters = FilterThresholds {
            high_mean_mismatches: Some(5.25),
            ..Default::default()
        };
        assert!(filters.any());
        let variant = TumorOnlyVariant {
            mean_mismatches_in_reads: mean_mismatches,
            ..variant(ref_allele, alt_allele, 30.0)
        };
        assert_eq!(filters.labels(&variant, ReadEvidence::default()), expected);
    }

    #[test]
    fn test_high_mean_mismatches_header_line_states_the_threshold() {
        let filters = FilterThresholds {
            high_mean_mismatches: Some(5.25),
            ..Default::default()
        };
        assert_eq!(
            filters.header_lines(),
            vec![
                r#"##FILTER=<ID=HIGH_MEAN_MISMATCHES,Description="Mean substitution mismatches per ALT read (FORMAT MEAN_MISMATCHES) is above 5.25, suggesting misaligned or paralogous reads.">"#
            ]
        );
    }

    #[rstest]
    #[case("G", "A", 0.0, 3, vec![SAME_READ_POSITION])]
    #[case("G", "A", 1.0, 3, vec![])]
    #[case("G", "A", 0.0, 1, vec![])]
    #[case("G", "A", 0.0, 40, vec![])]
    #[case("G", "G", 0.0, 3, vec![])]
    fn test_same_read_position(
        #[case] ref_allele: &'static str,
        #[case] alt_allele: &'static str,
        #[case] positions_vary: f32,
        #[case] alt_depth: i32,
        #[case] expected: Vec<&str>,
    ) {
        let filters = FilterThresholds {
            same_read_position: Some(0.35),
            ..Default::default()
        };
        assert!(filters.any());
        let variant = TumorOnlyVariant {
            stdev_position_in_read: positions_vary,
            alt_depth,
            depth: 100,
            ..variant(ref_allele, alt_allele, 30.0)
        };
        assert_eq!(filters.labels(&variant, ReadEvidence::default()), expected);
    }

    #[test]
    fn test_same_read_position_header_line_states_the_threshold() {
        let filters = FilterThresholds {
            same_read_position: Some(0.35),
            ..Default::default()
        };
        assert_eq!(
            filters.header_lines(),
            vec![
                r#"##FILTER=<ID=SAME_READ_POSITION,Description="Every one of at least 2 ALT reads places the variant at the same distance from its read end (FORMAT ALT_READ_POS_VARIES 0) at an AF (FORMAT AF) below 0.35, a sign of one fragment counted repeatedly or an error at a fixed sequencing cycle.">"#
            ]
        );
    }

    fn strand_variant(
        ref_strands: (i32, i32),
        alt_strands: (i32, i32),
    ) -> TumorOnlyVariant<'static> {
        TumorOnlyVariant {
            ref_forward: ref_strands.0,
            ref_reverse: ref_strands.1,
            alt_forward: alt_strands.0,
            alt_reverse: alt_strands.1,
            alt_depth: alt_strands.0 + alt_strands.1,
            depth: 1000,
            ..variant("G", "A", 30.0)
        }
    }

    #[rstest]
    #[case::an_empty_cell((500, 500), (0, 10), vec![STRAND_BIAS])]
    #[case::odds_ratio_above_the_minimum((500, 500), (2, 60), vec![STRAND_BIAS])]
    #[case::odds_ratio_folded((500, 500), (60, 2), vec![STRAND_BIAS])]
    #[case::odds_ratio_below_the_minimum((500, 500), (40, 120), vec![])]
    #[case::p_value_above_the_maximum((500, 500), (0, 3), vec![])]
    #[case::af_above_the_maximum((500, 500), (0, 300), vec![])]
    fn test_strand_bias(
        #[case] ref_strands: (i32, i32),
        #[case] alt_strands: (i32, i32),
        #[case] expected: Vec<&str>,
    ) {
        let filters = FilterThresholds {
            strand_bias: Some(0.01),
            strand_bias_min_odds_ratio: Some(5.0),
            strand_bias_max_af: Some(0.25),
            ..Default::default()
        };
        assert!(filters.any());
        assert_eq!(
            filters.labels(
                &strand_variant(ref_strands, alt_strands),
                ReadEvidence::default()
            ),
            expected
        );
    }

    #[test]
    fn test_strand_bias_without_optional_limits() {
        let filters = FilterThresholds {
            strand_bias: Some(0.01),
            ..Default::default()
        };
        let modest_skew_at_high_af = strand_variant((500, 500), (200, 400));
        assert_eq!(
            filters.labels(&modest_skew_at_high_af, ReadEvidence::default()),
            vec![STRAND_BIAS]
        );
    }

    #[test]
    fn test_strand_bias_is_not_applied_without_an_alt() {
        let filters = FilterThresholds {
            strand_bias: Some(0.01),
            ..Default::default()
        };
        let reference_row = TumorOnlyVariant {
            alt_allele: "G",
            ..strand_variant((500, 500), (0, 10))
        };
        assert!(
            filters
                .labels(&reference_row, ReadEvidence::default())
                .is_empty()
        );
    }

    #[test]
    fn test_strand_bias_header_line_states_the_thresholds() {
        let filters = FilterThresholds {
            strand_bias: Some(0.01),
            strand_bias_min_odds_ratio: Some(5.0),
            strand_bias_max_af: Some(0.25),
            ..Default::default()
        };
        assert_eq!(
            filters.header_lines(),
            vec![
                r#"##FILTER=<ID=STRAND_BIAS,Description="The ALT reads split across strands differently from the REF reads: the two-sided Fisher p-value of REF and ALT reads by strand (FORMAT STRAND_BIAS_FISHER_P) is below 0.01, the odds ratio of that table, folded so either direction counts, is above 5 or the table has an empty cell, and AF (FORMAT AF) is below 0.25.">"#
            ]
        );
        let p_only = FilterThresholds {
            strand_bias: Some(0.01),
            ..Default::default()
        };
        assert_eq!(
            p_only.header_lines(),
            vec![
                r#"##FILTER=<ID=STRAND_BIAS,Description="The ALT reads split across strands differently from the REF reads: the two-sided Fisher p-value of REF and ALT reads by strand (FORMAT STRAND_BIAS_FISHER_P) is below 0.01.">"#
            ]
        );
    }

    #[rstest]
    #[case("G", "A", 1, vec![LOW_AF])]
    #[case("G", "A", 2, vec![])]
    #[case("G", "A", 5, vec![])]
    #[case("G", "G", 1, vec![])]
    fn test_low_af(
        #[case] ref_allele: &'static str,
        #[case] alt_allele: &'static str,
        #[case] alt_depth: i32,
        #[case] expected: Vec<&str>,
    ) {
        let filters = FilterThresholds {
            low_af: Some(0.002),
            ..Default::default()
        };
        assert!(filters.any());
        let variant = TumorOnlyVariant {
            alt_depth,
            depth: 1000,
            ..variant(ref_allele, alt_allele, 30.0)
        };
        assert_eq!(filters.labels(&variant, ReadEvidence::default()), expected);
    }

    #[test]
    fn test_low_af_header_line_states_the_threshold() {
        let filters = FilterThresholds {
            low_af: Some(0.001),
            ..Default::default()
        };
        assert_eq!(
            filters.header_lines(),
            vec![
                r#"##FILTER=<ID=LOW_AF,Description="Fraction of reads carrying the ALT allele (FORMAT AF) is below 0.001.">"#
            ]
        );
    }

    #[rstest]
    #[case("G", "A", 20.0, vec![LOW_QMEAN])]
    #[case("G", "A", 30.0, vec![])]
    #[case("G", "G", 20.0, vec![])]
    fn test_low_qmean(
        #[case] ref_allele: &'static str,
        #[case] alt_allele: &'static str,
        #[case] qmean: f32,
        #[case] expected: Vec<&str>,
    ) {
        let filters = FilterThresholds {
            low_qmean: Some(30.0),
            ..Default::default()
        };
        assert!(filters.any());
        let variant = TumorOnlyVariant {
            base_quality_mean: qmean,
            ..variant(ref_allele, alt_allele, 30.0)
        };
        assert_eq!(filters.labels(&variant, ReadEvidence::default()), expected);
    }

    #[test]
    fn test_low_qmean_header_line_states_the_threshold() {
        let filters = FilterThresholds {
            low_qmean: Some(30.0),
            ..Default::default()
        };
        assert_eq!(
            filters.header_lines(),
            vec![
                r#"##FILTER=<ID=LOW_QMEAN,Description="Mean per-read variant base quality of the ALT reads (FORMAT QMEAN) is below 30.">"#
            ]
        );
    }

    #[rstest]
    #[case("G", "A", 2, vec![LOW_DP])]
    #[case("G", "A", 3, vec![])]
    #[case("G", "G", 2, vec![])]
    fn test_low_dp(
        #[case] ref_allele: &'static str,
        #[case] alt_allele: &'static str,
        #[case] depth: i32,
        #[case] expected: Vec<&str>,
    ) {
        let filters = FilterThresholds {
            low_dp: Some(3),
            ..Default::default()
        };
        assert!(filters.any());
        let variant = TumorOnlyVariant {
            depth,
            alt_depth: 1,
            ..variant(ref_allele, alt_allele, 30.0)
        };
        assert_eq!(filters.labels(&variant, ReadEvidence::default()), expected);
    }

    #[test]
    fn test_low_dp_header_line_states_the_threshold() {
        let filters = FilterThresholds {
            low_dp: Some(3),
            ..Default::default()
        };
        assert_eq!(
            filters.header_lines(),
            vec![r#"##FILTER=<ID=LOW_DP,Description="Read depth (FORMAT DP) is below 3.">"#]
        );
    }

    #[rstest]
    #[case("G", "A", 1, vec![LOW_HICNT])]
    #[case("G", "A", 2, vec![])]
    #[case("G", "A", i32::missing(), vec![])]
    #[case("G", "G", 0, vec![])]
    fn test_low_hicnt(
        #[case] ref_allele: &'static str,
        #[case] alt_allele: &'static str,
        #[case] hicnt: i32,
        #[case] expected: Vec<&str>,
    ) {
        let filters = FilterThresholds {
            low_hicnt: Some(2),
            ..Default::default()
        };
        assert!(filters.any());
        let variant = TumorOnlyVariant {
            alt_depth: 3,
            high_quality_variant_reads: hicnt,
            ..variant(ref_allele, alt_allele, 30.0)
        };
        assert_eq!(filters.labels(&variant, ReadEvidence::default()), expected);
    }

    #[test]
    fn test_low_hicnt_header_line_states_the_threshold() {
        let filters = FilterThresholds {
            low_hicnt: Some(2),
            ..Default::default()
        };
        assert_eq!(
            filters.header_lines(),
            vec![
                r#"##FILTER=<ID=LOW_HICNT,Description="High-quality ALT reads (FORMAT HICNT) are fewer than 2.">"#
            ]
        );
    }

    /// An SNV with this many high-quality ALT reads, as VarDict counted them, or none known, of
    /// which VarDict reassigned this many to it, at a depth of 1,000.
    fn snv_reassigned(hicnt: i32, reassigned: i32) -> TumorOnlyVariant<'static> {
        TumorOnlyVariant {
            depth: 1_000,
            alt_depth: hicnt.max(1),
            high_quality_variant_reads: hicnt,
            af_adjusted: reassigned as f32 / 1_000.0,
            ..variant("G", "A", 30.0)
        }
    }

    /// An SNV with this many high-quality ALT reads, one of them reassigned to it by VarDict.
    fn snv_with_hicnt(hicnt: i32) -> TumorOnlyVariant<'static> {
        snv_reassigned(hicnt, 1)
    }

    fn recount(aligned: i32, clean_clips: i32, foreign_clips: i32) -> Recount {
        Recount {
            aligned,
            clean_clips,
            foreign_clips,
        }
    }

    #[rstest]
    #[case::one_aligned_read_and_a_foreign_clip(2, Some(recount(1, 0, 1)), vec![FOREIGN_CLIP_SUPPORT])]
    #[case::one_aligned_read_and_a_clean_clip(2, Some(recount(1, 1, 0)), vec![])]
    #[case::two_aligned_reads(2, Some(recount(2, 0, 1)), vec![])]
    #[case::a_shortfall_no_clip_explains(2, Some(recount(1, 0, 0)), vec![FOREIGN_CLIP_SUPPORT])]
    #[case::no_read_at_all(2, Some(recount(0, 0, 0)), vec![FOREIGN_CLIP_SUPPORT])]
    #[case::vardict_counted_under_the_floor(1, Some(recount(0, 0, 1)), vec![LOW_HICNT])]
    #[case::hicnt_unknown(i32::missing(), Some(recount(1, 0, 1)), vec![])]
    #[case::not_recounted(2, None, vec![])]
    fn test_foreign_clip_support(
        #[case] hicnt: i32,
        #[case] recount: Option<Recount>,
        #[case] expected: Vec<&str>,
    ) {
        let filters = FilterThresholds {
            low_hicnt: Some(2),
            foreign_clip_support: true,
            ..Default::default()
        };
        assert!(filters.any());
        assert_eq!(
            filters.labels(
                &snv_with_hicnt(hicnt),
                ReadEvidence {
                    recount,
                    ..Default::default()
                }
            ),
            expected
        );
    }

    #[rstest]
    #[case::two_of_three_at_a_floor_of_three(3, 3, Some(recount(2, 0, 1)), vec![FOREIGN_CLIP_SUPPORT])]
    #[case::two_of_three_at_a_floor_of_two(2, 3, Some(recount(2, 0, 1)), vec![])]
    fn test_foreign_clip_support_at_the_low_hicnt_floor(
        #[case] min_hicnt: i32,
        #[case] hicnt: i32,
        #[case] recount: Option<Recount>,
        #[case] expected: Vec<&str>,
    ) {
        let filters = FilterThresholds {
            low_hicnt: Some(min_hicnt),
            foreign_clip_support: true,
            ..Default::default()
        };
        assert_eq!(filters.foreign_clip_support_floor(), Some(min_hicnt));
        assert_eq!(
            filters.labels(
                &snv_with_hicnt(hicnt),
                ReadEvidence {
                    recount,
                    ..Default::default()
                }
            ),
            expected
        );
    }

    #[test]
    fn test_foreign_clip_support_is_not_applied_without_the_low_hicnt_floor() {
        let filters = FilterThresholds {
            foreign_clip_support: true,
            ..Default::default()
        };
        assert!(!filters.any());
        assert_eq!(filters.foreign_clip_support_floor(), None);
        assert!(!filters.recounts(&snv_with_hicnt(2)));
        assert!(filters.header_lines().is_empty());
        assert_eq!(
            filters.labels(
                &snv_with_hicnt(2),
                ReadEvidence {
                    recount: Some(recount(0, 0, 1)),
                    ..Default::default()
                }
            ),
            Vec::<&str>::new()
        );
    }

    #[rstest]
    #[case::an_snv_at_the_floor("G", "A", 2, 1, true)]
    #[case::an_snv_at_the_floor_with_no_read_reassigned("G", "A", 2, 0, false)]
    #[case::an_snv_its_reassigned_reads_could_have_lifted_to_the_floor("G", "A", 4, 3, true)]
    #[case::an_snv_over_the_floor_without_its_reassigned_reads("G", "A", 5, 2, false)]
    #[case::an_snv_under_the_floor("G", "A", 1, 1, false)]
    #[case::an_snv_with_no_count_known("G", "A", i32::missing(), 1, false)]
    #[case::an_mnv("GG", "AT", 2, 1, false)]
    #[case::a_deletion("GA", "G", 2, 1, false)]
    #[case::a_reference_row("G", "G", 2, 1, false)]
    #[case::a_symbolic_allele("G", "<DEL>", 2, 1, false)]
    fn test_only_snvs_reassigned_reads_could_have_lifted_to_the_floor_are_recounted(
        #[case] ref_allele: &'static str,
        #[case] alt_allele: &'static str,
        #[case] hicnt: i32,
        #[case] reassigned: i32,
        #[case] expected: bool,
    ) {
        let filters = FilterThresholds {
            low_hicnt: Some(2),
            foreign_clip_support: true,
            ..Default::default()
        };
        let call = TumorOnlyVariant {
            ref_allele,
            alt_allele,
            ..snv_reassigned(hicnt, reassigned)
        };
        assert_eq!(filters.recounts(&call), expected);
        assert!(!FilterThresholds::default().recounts(&call));
    }

    #[test]
    fn test_foreign_clip_support_header_line_states_the_threshold() {
        let filters = FilterThresholds {
            low_hicnt: Some(2),
            foreign_clip_support: true,
            ..Default::default()
        };
        assert_eq!(
            filters.header_lines(),
            vec![
                r#"##FILTER=<ID=LOW_HICNT,Description="High-quality ALT reads (FORMAT HICNT) are fewer than 2.">"#,
                r#"##FILTER=<ID=FOREIGN_CLIP_SUPPORT,Description="SNV with at least 2 high-quality ALT reads by VarDict (FORMAT HICNT), a count the reads it reassigned to the allele (FORMAT REALIGNED_FRAC_OF_DP) could have lifted to 2, but fewer than 2 in the reads it was run on: reads carrying the ALT base aligned at the site at or above VarDict's -q, once its clipping of mismatches at read ends is replayed, or first in a soft clip at the site whose later bases match the reference. The shortfall is most often a soft clip VarDict added after comparing only its second base to the reference. Every count is left as VarDict wrote it.">"#,
            ]
        );
    }

    /// An insertion of this many bases with these high-quality ALT reads, as VarDict counted them,
    /// or none known, and these REF reads.
    fn insertion(inserted: usize, hicnt: i32, ref_reads: i32) -> TumorOnlyVariant<'static> {
        let alt_allele: &'static str = Box::leak(
            format!("G{}", "ACGT".repeat(inserted.div_ceil(4)))[..inserted + 1]
                .to_string()
                .into_boxed_str(),
        );
        TumorOnlyVariant {
            alt_depth: hicnt.max(1),
            high_quality_variant_reads: hicnt,
            ref_forward: ref_reads,
            ref_reverse: 0,
            ..variant("G", alt_allele, 30.0)
        }
    }

    fn low_support() -> FilterThresholds {
        FilterThresholds {
            low_support_long_insertion: Some(LowSupportLongInsertion {
                min_length: 50,
                min_hicnt: 8,
            }),
            ..Default::default()
        }
    }

    #[rstest]
    #[case::few_reads_and_no_reads_given(insertion(60, 3, 10), None, vec![LOW_SUPPORT_LONG_INSERTION])]
    #[case::few_reads_and_no_read_carries_it(insertion(60, 3, 10), Some(0), vec![LOW_SUPPORT_LONG_INSERTION])]
    #[case::a_read_carries_it(insertion(60, 3, 10), Some(1), vec![])]
    #[case::at_the_length(insertion(50, 7, 10), None, vec![LOW_SUPPORT_LONG_INSERTION])]
    #[case::under_the_length(insertion(49, 3, 10), None, vec![])]
    #[case::at_the_floor(insertion(60, 8, 10), None, vec![])]
    #[case::no_ref_reads(insertion(60, 3, 0), None, vec![])]
    #[case::hicnt_unknown(insertion(60, i32::missing(), 10), None, vec![])]
    #[case::a_deletion(TumorOnlyVariant { ref_allele: insertion(60, 3, 10).alt_allele, ..insertion(60, 3, 10) }, None, vec![])]
    #[case::a_complex_call_adding_the_length(TumorOnlyVariant { ref_allele: "GA", variant_type: "Complex", ..insertion(51, 3, 10) }, None, vec![LOW_SUPPORT_LONG_INSERTION])]
    #[case::a_complex_call_adding_less(TumorOnlyVariant { ref_allele: "GA", variant_type: "Complex", ..insertion(50, 3, 10) }, None, vec![])]
    #[case::a_complex_call_a_read_spells(TumorOnlyVariant { ref_allele: "GA", variant_type: "Complex", ..insertion(51, 3, 10) }, Some(1), vec![])]
    #[case::a_symbolic_allele(TumorOnlyVariant { alt_allele: "<DUP>", ..insertion(60, 3, 10) }, None, vec![])]
    fn test_low_support_long_insertion(
        #[case] call: TumorOnlyVariant<'static>,
        #[case] spanning: Option<i32>,
        #[case] expected: Vec<&str>,
    ) {
        let filters = low_support();
        assert!(filters.any());
        let reads = ReadEvidence {
            spanning,
            ..Default::default()
        };
        assert_eq!(filters.labels(&call, reads), expected);
        assert_eq!(
            filters.low_support_long_insertion.unwrap().tests(&call),
            spanning.is_some() || !expected.is_empty()
        );
        assert!(FilterThresholds::default().labels(&call, reads).is_empty());
    }

    #[test]
    fn test_low_support_long_insertion_beside_low_hicnt() {
        let filters = FilterThresholds {
            low_hicnt: Some(2),
            ..low_support()
        };
        assert_eq!(
            filters.labels(&insertion(60, 1, 10), ReadEvidence::default()),
            vec![LOW_HICNT, LOW_SUPPORT_LONG_INSERTION]
        );
    }

    #[test]
    fn test_low_support_long_insertion_header_line_states_the_thresholds() {
        assert_eq!(
            low_support().header_lines(),
            vec![
                r#"##FILTER=<ID=LOW_SUPPORT_LONG_INSERTION,Description="Insertion, or Complex call longer than its REF, adding at least 50 bases to the reference, with fewer than 8 high-quality ALT reads (FORMAT HICNT), at least one read VarDict counted as REF (FORMAT AD, left missing for Complex calls) and, when the reads VarDict ran on are given, no read carrying it whole: for an insertion, every added base in one insertion operation at a mean quality at or above VarDict's -q, at the call's position or one its insertion slides to; for a Complex call, its ALT spelled by the read's bases over its REF and the insertion operations among them. Every count is left as VarDict wrote it.">"#
            ]
        );
    }

    #[test]
    fn test_labels_are_listed_in_a_fixed_order() {
        let filters = FilterThresholds {
            near_read_end: Some(8.0),
            low_mean_mapq: Some(10.0),
            ..Default::default()
        };
        let failing_both = TumorOnlyVariant {
            mean_mapping_quality: 5.0,
            ..variant("G", "A", 3.0)
        };
        assert_eq!(
            filters.labels(&failing_both, ReadEvidence::default()),
            vec![NEAR_READ_END, LOW_MEAN_MAPQ]
        );
    }
}
