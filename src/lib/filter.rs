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
}

impl FilterThresholds {
    /// Whether any FILTER label is applied.
    pub fn any(&self) -> bool {
        self.near_read_end.is_some()
            || self.low_mean_mapq.is_some()
            || self.homopolymer_indel.is_some()
            || self.tandem_repeat_indel.is_some()
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
        lines
    }

    /// Return the FILTER labels a call fails; none for a record without an ALT allele.
    pub fn labels(&self, variant: &TumorOnlyVariant) -> Vec<&'static str> {
        let mut labels = Vec::new();
        if variant.ref_allele == variant.alt_allele {
            return labels;
        }
        if self
            .near_read_end
            .is_some_and(|min_mean_dist| variant.mean_position_in_read < min_mean_dist)
        {
            labels.push(NEAR_READ_END);
        }
        if self
            .low_mean_mapq
            .is_some_and(|min_mean_mapq| variant.mean_mapping_quality < min_mean_mapq)
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
        labels
    }
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
            ..Default::default()
        }
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
