//! Fisher's exact test on a 2x2 table of read counts.

/// Return the one-sided Fisher exact p-value that the first row's share of the first column is
/// greater than the second row's, for the table `[[a, b], [c, d]]`: R's
/// `fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")`.
pub fn fisher_exact_greater(a: u64, b: u64, c: u64, d: u64) -> f64 {
    let (row, col, total) = ((a + b) as f64, (a + c) as f64, (a + b + c + d) as f64);
    let mode = ((row + 1.0) * (col + 1.0) / (total + 2.0)).floor();
    let ln_pmf =
        |x: f64| ln_choose(col, x) + ln_choose(total - col, row - x) - ln_choose(total, row);
    let other = d as f64 - a as f64;
    if a as f64 > mode {
        let upper = (a..=(a + b).min(a + c)).map(|x| x as f64);
        let ln_ratio = |x: f64| ((col - x) * (row - x)).ln() - ((x + 1.0) * (other + x + 1.0)).ln();
        return falling_sum(ln_pmf(a as f64), upper, ln_ratio).min(1.0);
    }
    let lowest = a.saturating_sub(d);
    if a == lowest {
        return 1.0;
    }
    let lower = (lowest..a).rev().map(|x| x as f64);
    let ln_ratio = |x: f64| (x * (other + x)).ln() - ((col - x + 1.0) * (row - x + 1.0)).ln();
    (1.0 - falling_sum(ln_pmf(a as f64 - 1.0), lower, ln_ratio)).max(0.0)
}

/// Sums terms that only fall, given the log of the first and the log ratio of each term to the
/// next, stopping once they drop below e^-40 of the first.
fn falling_sum(ln_first: f64, xs: impl Iterator<Item = f64>, ln_ratio: impl Fn(f64) -> f64) -> f64 {
    let (mut ln_term, mut sum) = (ln_first, 0.0);
    for x in xs {
        if ln_term < ln_first - 40.0 {
            break;
        }
        sum += (ln_term - ln_first).exp();
        ln_term += ln_ratio(x);
    }
    (ln_first + sum.ln()).exp()
}

/// The natural log of the binomial coefficient `n` choose `k`.
fn ln_choose(n: f64, k: f64) -> f64 {
    ln_gamma(n + 1.0) - ln_gamma(k + 1.0) - ln_gamma(n - k + 1.0)
}

/// The natural log of the gamma function for `x >= 1`, by the Lanczos approximation (g = 7).
fn ln_gamma(x: f64) -> f64 {
    const COEFFICIENTS: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    let x = x - 1.0;
    let series = COEFFICIENTS[1..]
        .iter()
        .enumerate()
        .fold(COEFFICIENTS[0], |acc, (i, coefficient)| {
            acc + coefficient / (x + i as f64 + 1.0)
        });
    let t = x + 7.5;
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + series.ln()
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case(10, 20, 0, 30, 0.0003985065657050687)]
    #[case(10, 20, 1, 30, 0.0023583749216316376)]
    #[case(15, 15, 15, 15, 0.6017271810907904)]
    #[case(30, 0, 15, 15, 2.916066847916421e-06)]
    #[case(10, 20, 3, 27, 0.028733136666041024)]
    #[case(0, 30, 15, 15, 1.0)]
    #[case(1, 30, 15, 15, 0.9999992829343817)]
    #[case(3, 27, 15, 15, 0.999927636742836)]
    #[case(10, 20, 0, 0, 1.0)]
    #[case(0, 0, 15, 15, 1.0)]
    #[case(2, 28, 0, 30, 0.2457627118644069)]
    #[case(10, 20, 2, 48, 0.0006576798232705035)]
    #[case(400, 7704, 0, 6000, 7.500022210244021e-99)]
    #[case(1, 8103, 0, 6000, 0.5745887691435053)]
    #[case(0, 0, 0, 0, 1.0)]
    #[case(5, 0, 0, 10, 0.00033300033300033305)]
    #[case(1, 29999, 15000, 0, 1.0)]
    #[case(5, 95, 50, 50, 0.9999999999999948)]
    #[case(100, 900, 500, 500, 1.0)]
    #[case(10, 990, 300, 700, 1.0)]
    #[case(2, 8, 5, 5, 0.9713622291021672)]
    #[case(3, 2, 4, 1, 0.9166666666666666)]
    #[case(250, 250, 260, 240, 0.7567260484183103)]
    #[case(1, 0, 1, 0, 1.0)]
    #[case(0, 5, 0, 5, 1.0)]
    fn test_fisher_exact_greater_matches_r(
        #[case] a: u64,
        #[case] b: u64,
        #[case] c: u64,
        #[case] d: u64,
        #[case] expected: f64,
    ) {
        let p = fisher_exact_greater(a, b, c, d);
        assert!(
            (p - expected).abs() <= expected * 1e-9,
            "fisher_exact_greater({a}, {b}, {c}, {d}) = {p}, expected {expected}"
        );
    }

    #[test]
    fn test_fisher_exact_greater_underflows_to_zero_like_r() {
        assert_eq!(fisher_exact_greater(4000, 4104, 10, 5990), 0.0);
    }
}
