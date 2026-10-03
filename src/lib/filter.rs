//! FILTER labels for VarDict calls.

use crate::record::TumorOnlyVariant;

/// The FILTER label for calls whose ALT reads carry the variant near their read ends.
pub const NEAR_READ_END: &str = "NEAR_READ_END";

/// The thresholds of the FILTER labels to apply; a label without a threshold is not applied.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FilterThresholds {
    /// The minimum mean distance from the variant to the nearer read end.
    pub near_read_end: Option<f32>,
}

impl FilterThresholds {
    /// Whether any FILTER label is applied.
    pub fn any(&self) -> bool {
        self.near_read_end.is_some()
    }

    /// Return the FILTER header lines of the applied labels.
    pub fn header_lines(&self) -> Vec<String> {
        let mut lines = Vec::new();
        if let Some(min_mean_dist) = self.near_read_end {
            lines.push(format!(
                r#"##FILTER=<ID={NEAR_READ_END},Description="Mean distance from the variant to the nearer end of the ALT reads' aligned parts (FORMAT MEAN_DIST_TO_READ_END) is below {min_mean_dist}.">"#
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
        labels
    }
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
            ..Default::default()
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
        };
        assert_eq!(
            filters.header_lines(),
            vec![
                r#"##FILTER=<ID=NEAR_READ_END,Description="Mean distance from the variant to the nearer end of the ALT reads' aligned parts (FORMAT MEAN_DIST_TO_READ_END) is below 7.5.">"#
            ]
        );
    }
}
