//! FILTER labels for VarDict calls.

use crate::record::TumorOnlyVariant;

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
    }

    /// Return the FILTER header lines of the applied labels.
    pub fn header_lines(&self) -> Vec<String> {
        let mut lines = Vec::new();
        if let Some(min_mean_dist) = self.near_read_end {
            lines.push(format!(
                r#"##FILTER=<ID={NEAR_READ_END},Description="Mean distance from the variant to the nearer end of the ALT reads' aligned parts (FORMAT MEAN_DIST_TO_READ_END) is below {min_mean_dist}.">"#
            ));
        }
        if let Some(min_mean_mapq) = self.low_mean_mapq {
            lines.push(format!(
                r#"##FILTER=<ID={LOW_MEAN_MAPQ},Description="Mean mapping quality of the ALT reads (FORMAT MEAN_MAPQ) is below {min_mean_mapq}.">"#
            ));
        }
        if let Some(min_copies) = self.homopolymer_indel {
            let af_limit = af_limit_clause(self.homopolymer_indel_max_af);
            lines.push(format!(
                r#"##FILTER=<ID={HOMOPOLYMER_INDEL},Description="Insertion or deletion of exactly one base in a homopolymer of at least {min_copies} copies (INFO REPEAT_UNIT_COPIES with REPEAT_UNIT_LEN 1){af_limit}, the signature of polymerase slippage.">"#
            ));
        }
        if let Some(min_copies) = self.tandem_repeat_indel {
            let af_limit = af_limit_clause(self.tandem_repeat_indel_max_af);
            lines.push(format!(
                r#"##FILTER=<ID={TANDEM_REPEAT_INDEL},Description="Insertion or deletion of exactly one repeat unit in a 2-6 bp tandem repeat of at least {min_copies} copies (INFO REPEAT_UNIT_COPIES and REPEAT_UNIT_LEN){af_limit}, the signature of polymerase slippage.">"#
            ));
        }
        if let Some(max_mean_mismatches) = self.high_mean_mismatches {
            lines.push(format!(
                r#"##FILTER=<ID={HIGH_MEAN_MISMATCHES},Description="Mean substitution mismatches per ALT read (FORMAT MEAN_MISMATCHES) is above {max_mean_mismatches}, suggesting misaligned or paralogous reads.">"#
            ));
        }
        if let Some(max_af) = self.same_read_position {
            lines.push(format!(
                r#"##FILTER=<ID={SAME_READ_POSITION},Description="Every one of at least 2 ALT reads places the variant at the same distance from its read end (FORMAT ALT_READ_POS_VARIES 0) at an AF (FORMAT AF) below {max_af}, a sign of one fragment counted repeatedly or an error at a fixed sequencing cycle.">"#
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
            lines.push(format!(
                r#"##FILTER=<ID={STRAND_BIAS},Description="The ALT reads split across strands differently from the REF reads: the two-sided Fisher p-value of REF and ALT reads by strand (FORMAT STRAND_BIAS_FISHER_P) is below {max_p}{odds_ratio}{af_limit}.">"#
            ));
        }
        if let Some(min_af) = self.low_af {
            lines.push(format!(
                r#"##FILTER=<ID={LOW_AF},Description="Fraction of reads carrying the ALT allele (FORMAT AF) is below {min_af}.">"#
            ));
        }
        if let Some(min_qmean) = self.low_qmean {
            lines.push(format!(
                r#"##FILTER=<ID={LOW_QMEAN},Description="Mean per-read variant base quality of the ALT reads (FORMAT QMEAN) is below {min_qmean}.">"#
            ));
        }
        if let Some(min_dp) = self.low_dp {
            lines.push(format!(
                r#"##FILTER=<ID={LOW_DP},Description="Read depth (FORMAT DP) is below {min_dp}.">"#
            ));
        }
        lines
    }

    /// Return the FILTER labels a call fails by its FORMAT values: none without an ALT allele, and
    /// none from a value that is missing.
    pub fn labels(&self, variant: &TumorOnlyVariant) -> Vec<&'static str> {
        let mut labels = Vec::new();
        if variant.ref_allele == variant.alt_allele {
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
            variant.strand_bias_fisher_p_value() < max_p
                && self
                    .strand_bias_min_odds_ratio
                    .is_none_or(|min_odds_ratio| strand_odds_ratio_exceeds(variant, min_odds_ratio))
                && self
                    .strand_bias_max_af
                    .is_none_or(|max_af| variant.af_value() < max_af)
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
    fn test_a_call_without_alt_reads_gets_only_count_labels() {
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
            strand_bias_p_value: 0.0,
            ..variant("G", "A", 0.0)
        };
        assert_eq!(filters.labels(&no_alt_reads), vec![LOW_AF, LOW_DP]);
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
        assert!(filters.labels(&variant("G", "A", 1.0)).is_empty());
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
            filters.labels(&variant(ref_allele, alt_allele, mean_dist)),
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
            filters.labels(&mapq_variant(ref_allele, alt_allele, mean_mapq)),
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
        assert_eq!(filters.labels(&variant), expected);
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
        assert_eq!(filters.labels(&variant), expected);
    }

    #[test]
    fn test_repeat_indels_have_no_af_limit_without_a_max_af() {
        let filters = FilterThresholds {
            homopolymer_indel: Some(13.0),
            tandem_repeat_indel: Some(13.0),
            ..Default::default()
        };
        assert_eq!(
            filters.labels(&repeat_variant("GA", "G", 1, 13.0, 90)),
            vec![HOMOPOLYMER_INDEL]
        );
        assert_eq!(
            filters.labels(&repeat_variant("GCA", "G", 2, 13.0, 90)),
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
                r#"##FILTER=<ID=HOMOPOLYMER_INDEL,Description="Insertion or deletion of exactly one base in a homopolymer of at least 13 copies (INFO REPEAT_UNIT_COPIES with REPEAT_UNIT_LEN 1) at an AF (FORMAT AF) below 0.275, the signature of polymerase slippage.">"#,
                r#"##FILTER=<ID=TANDEM_REPEAT_INDEL,Description="Insertion or deletion of exactly one repeat unit in a 2-6 bp tandem repeat of at least 13 copies (INFO REPEAT_UNIT_COPIES and REPEAT_UNIT_LEN), the signature of polymerase slippage.">"#,
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
        assert_eq!(filters.labels(&variant), expected);
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
        assert_eq!(filters.labels(&variant), expected);
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
        p_value: f32,
        ref_strands: (i32, i32),
        alt_strands: (i32, i32),
    ) -> TumorOnlyVariant<'static> {
        TumorOnlyVariant {
            strand_bias_p_value: p_value,
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
    #[case(0.001, (500, 500), (0, 10), vec![STRAND_BIAS])]
    #[case(0.001, (500, 500), (2, 60), vec![STRAND_BIAS])]
    #[case(0.001, (500, 500), (60, 2), vec![STRAND_BIAS])]
    #[case(0.001, (500, 500), (20, 40), vec![])]
    #[case(0.5, (500, 500), (0, 10), vec![])]
    #[case(0.001, (500, 500), (0, 300), vec![])]
    fn test_strand_bias(
        #[case] p_value: f32,
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
            filters.labels(&strand_variant(p_value, ref_strands, alt_strands)),
            expected
        );
    }

    #[test]
    fn test_strand_bias_without_optional_limits() {
        let filters = FilterThresholds {
            strand_bias: Some(0.01),
            ..Default::default()
        };
        let modest_skew_at_high_af = strand_variant(0.001, (500, 500), (200, 400));
        assert_eq!(filters.labels(&modest_skew_at_high_af), vec![STRAND_BIAS]);
    }

    #[test]
    fn test_strand_bias_is_not_applied_without_an_alt() {
        let filters = FilterThresholds {
            strand_bias: Some(0.01),
            ..Default::default()
        };
        let reference_row = TumorOnlyVariant {
            alt_allele: "G",
            ..strand_variant(0.001, (500, 500), (0, 10))
        };
        assert!(filters.labels(&reference_row).is_empty());
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
        assert_eq!(filters.labels(&variant), expected);
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
        assert_eq!(filters.labels(&variant), expected);
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
        assert_eq!(filters.labels(&variant), expected);
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
            filters.labels(&failing_both),
            vec![NEAR_READ_END, LOW_MEAN_MAPQ]
        );
    }
}
