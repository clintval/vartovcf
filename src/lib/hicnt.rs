//! Recovers VarDict's count of ALT reads at or above its `-q` (its hicnt), which VarDict prints in
//! tumor-only output but not for a tumor-normal pair, from the SN it prints for each sample.

/// Returns the one hicnt that reproduces the SN VarDict printed for an allele with `alt_depth` ALT
/// reads, or `None` when no count or more than one does.
///
/// VarDict computes SN as `hicnt / locnt`, or `hicnt / 0.5` when `locnt` is 0, where `locnt` is the
/// ALT reads below `-q`, and prints the double to 3 decimals with Java's `DecimalFormat`. For a fixed
/// total of reads that ratio rises with hicnt, so the counts reproducing SN form one run, found by
/// bisection.
/// The total is usually the ALT depth, but when VarDict realigns an insertion it scales the ALT
/// depth, hicnt and locnt separately, rounding each toward zero, which can leave the ALT depth one
/// above `hicnt + locnt`; both totals are tried, and a count is returned only when they agree.
/// More than one count reproduces SN only when a thousand or more ALT reads fall below `-q`.
pub fn hicnt_from_sn(alt_depth: u32, sn_thousandths: u64) -> Option<u32> {
    let mut found = None;
    for total in [Some(alt_depth), alt_depth.checked_sub(1)]
        .into_iter()
        .flatten()
    {
        match counts_reproducing(total, sn_thousandths) {
            Fit::None => {}
            Fit::One(hicnt) if found.is_none_or(|other| other == hicnt) => found = Some(hicnt),
            Fit::One(_) | Fit::Several => return None,
        }
    }
    found
}

/// Parses a non-negative decimal with at most 3 decimals, the way VarDict prints SN (`20`, `0.333`,
/// or `2.` once it strips the zeros from `2.000`), into exact thousandths.
pub fn parse_thousandths(text: &str) -> Option<u64> {
    let (whole, fraction) = text.split_once('.').unwrap_or((text, ""));
    let digits = |part: &str| part.bytes().all(|byte| byte.is_ascii_digit());
    if whole.is_empty() || fraction.len() > 3 || !digits(whole) || !digits(fraction) {
        return None;
    }
    let fraction = fraction
        .bytes()
        .chain(std::iter::repeat(b'0'))
        .take(3)
        .fold(0, |sum, byte| sum * 10 + u64::from(byte - b'0'));
    whole
        .parse::<u64>()
        .ok()?
        .checked_mul(1000)?
        .checked_add(fraction)
}

/// The hicnt values that reproduce an SN for one total of ALT reads.
#[derive(Debug, PartialEq)]
enum Fit {
    None,
    One(u32),
    Several,
}

/// Finds the hicnt values from 0 to `total` whose SN, with the rest of the total as locnt, rounds to
/// `sn_thousandths`.
fn counts_reproducing(total: u32, sn_thousandths: u64) -> Fit {
    let sn = |hicnt: u32| sn_thousandths_of(hicnt, total - hicnt);
    let (mut low, mut high) = (0, total);
    while low < high {
        let middle = low + (high - low) / 2;
        if sn(middle) < sn_thousandths {
            low = middle + 1;
        } else {
            high = middle;
        }
    }
    if sn(low) != sn_thousandths {
        Fit::None
    } else if low < total && sn(low + 1) == sn_thousandths {
        Fit::Several
    } else {
        Fit::One(low)
    }
}

/// VarDict's SN for these counts in thousandths: the double `hicnt / locnt`, or `hicnt / 0.5` when
/// `locnt` is 0, rounded to 3 decimals as Java's `DecimalFormat` rounds it, which is half-even on
/// its exact binary value except that the double nearest 0.0005, just above it, prints as 0.
fn sn_thousandths_of(hicnt: u32, locnt: u32) -> u64 {
    let divisor = if locnt == 0 { 0.5 } else { f64::from(locnt) };
    let ratio = f64::from(hicnt) / divisor;
    if ratio == 0.0005 {
        return 0;
    }
    thousandths_half_even(ratio)
}

/// Rounds a finite, non-negative double half-even to a whole number of thousandths, settling ties
/// on its exact binary value rather than on a decimal approximation of it.
fn thousandths_half_even(value: f64) -> u64 {
    let bits = value.to_bits();
    let biased_exponent = (bits >> 52) as i32;
    if biased_exponent == 0 {
        return 0;
    }
    let shift = 1075 - biased_exponent;
    if shift <= 0 {
        return (value as u64).saturating_mul(1000);
    }
    if shift >= 64 {
        return 0;
    }
    let scaled = u128::from(bits & ((1 << 52) - 1) | 1 << 52) * 1000;
    let quotient = scaled >> shift;
    let remainder = scaled & ((1 << shift) - 1);
    let half = 1 << (shift - 1);
    let round_up = remainder > half || (remainder == half && quotient % 2 == 1);
    (quotient + u128::from(round_up)) as u64
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case(10, 20_000, Some(10))]
    #[case(1, 2_000, Some(1))]
    #[case(3, 0, Some(0))]
    #[case(22, 4_500, Some(18))]
    #[case(22, 10_000, Some(20))]
    #[case(8, 1_333, Some(4))]
    #[case(8, 333, Some(2))]
    #[case(10, 1_234, None)]
    #[case(0, 0, Some(0))]
    #[case(0, 60_000, None)]
    fn test_hicnt_from_sn(
        #[case] alt_depth: u32,
        #[case] sn_thousandths: u64,
        #[case] expected: Option<u32>,
    ) {
        assert_eq!(hicnt_from_sn(alt_depth, sn_thousandths), expected);
    }

    #[test]
    fn test_hicnt_from_sn_when_locnt_is_zero() {
        for alt_depth in 1..=500 {
            let sn = sn_thousandths_of(alt_depth, 0);
            assert_eq!(sn, 2_000 * u64::from(alt_depth));
            assert_eq!(hicnt_from_sn(alt_depth, sn), Some(alt_depth));
        }
    }

    #[test]
    fn test_hicnt_from_sn_reproduces_every_count_below_a_thousand_low_quality_reads() {
        for alt_depth in 0..=1_000 {
            for hicnt in 0..=alt_depth {
                let sn = sn_thousandths_of(hicnt, alt_depth - hicnt);
                assert_eq!(
                    hicnt_from_sn(alt_depth, sn),
                    Some(hicnt),
                    "{hicnt} of {alt_depth}"
                );
            }
        }
    }

    #[test]
    fn test_hicnt_from_sn_when_the_alt_depth_is_one_above_hicnt_plus_locnt() {
        for total in 0..=1_000 {
            for hicnt in 0..=total {
                let sn = sn_thousandths_of(hicnt, total - hicnt);
                assert_eq!(
                    hicnt_from_sn(total + 1, sn),
                    Some(hicnt),
                    "{hicnt} of {total} + 1"
                );
            }
        }
    }

    #[test]
    fn test_hicnt_from_sn_is_missing_when_counts_round_alike() {
        assert_eq!(sn_thousandths_of(1, 1_999), 1);
        assert_eq!(sn_thousandths_of(2, 1_998), 1);
        assert_eq!(hicnt_from_sn(2_000, 1), None);
        assert_eq!(counts_reproducing(2_000, 1), Fit::Several);
    }

    #[test]
    fn test_hicnt_from_sn_is_missing_when_the_two_totals_disagree() {
        assert_eq!(sn_thousandths_of(501, 1_001), 500);
        assert_eq!(counts_reproducing(1_502, 500), Fit::One(501));
        assert_eq!(counts_reproducing(1_501, 500), Fit::One(500));
        assert_eq!(hicnt_from_sn(1_502, 500), None);
    }

    #[rstest]
    #[case(10, 0, 20_000)]
    #[case(0, 5, 0)]
    #[case(1, 3, 333)]
    #[case(1, 16, 62)]
    #[case(3, 16, 188)]
    #[case(1, 1_999, 1)]
    #[case(1, 2_000, 0)]
    #[case(1, 2_001, 0)]
    #[case(1_999, 2_000, 1_000)]
    fn test_sn_thousandths_of(#[case] hicnt: u32, #[case] locnt: u32, #[case] expected: u64) {
        assert_eq!(sn_thousandths_of(hicnt, locnt), expected);
    }

    #[rstest]
    #[case(0.0, 0)]
    #[case(20.0, 20_000)]
    #[case(1.0 / 3.0, 333)]
    #[case(2.0 / 3.0, 667)]
    #[case(0.0625, 62)]
    #[case(0.1875, 188)]
    #[case(1.0 / 2_000.0, 1)]
    #[case(1.0 / 1_999.0, 1)]
    #[case(1.0 / 2_001.0, 0)]
    #[case(2_147_483_647.0, 2_147_483_647_000)]
    #[case(1e-300, 0)]
    fn test_thousandths_half_even(#[case] value: f64, #[case] expected: u64) {
        assert_eq!(thousandths_half_even(value), expected);
    }

    #[rstest]
    #[case("0", Some(0))]
    #[case("20", Some(20_000))]
    #[case("0.333", Some(333))]
    #[case("1.5", Some(1_500))]
    #[case("2.", Some(2_000))]
    #[case("0.", Some(0))]
    #[case("12345.679", Some(12_345_679))]
    #[case("0.3333", None)]
    #[case(".5", None)]
    #[case("-1", None)]
    #[case("1e3", None)]
    #[case("", None)]
    #[case("NaN", None)]
    fn test_parse_thousandths(#[case] text: &str, #[case] expected: Option<u64>) {
        assert_eq!(parse_thousandths(text), expected);
    }
}
