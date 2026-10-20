//! Normalizes a call against the reference as `bcftools norm -f` does for a record with one ALT
//! allele (Tan, Abecasis and Kang 2015, "Unified representation of genetic variants").
//!
//! The REF allele must match the reference. A call's alleles lose their shared trailing bases, then
//! their shared leading bases down to one, and an insertion or deletion moves left while the base
//! before it equals its last base, keeping one base before it as its anchor. A substitution is only
//! trimmed. Every allele is written uppercase. Reference-only calls and symbolic structural variants
//! keep their alleles and position.

use crate::record::{TumorOnlyVariant, variant_class};
use crate::reference::Reference;

/// A call's position and alleles after normalization.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Normalized {
    /// The 1-based position of the first REF base.
    pub start: u64,
    /// The 1-based position of the last REF base, or VarDict's END for a symbolic allele.
    pub end: u64,
    /// The REF allele, uppercase.
    pub ref_allele: String,
    /// The ALT allele, uppercase unless symbolic, and equal to REF for a reference-only call.
    pub alt_allele: String,
    /// The bases an insertion or deletion can slide toward the 3' end, `None` for other calls.
    pub indel_3p_shift: Option<i32>,
}

impl Normalized {
    /// Returns a call with this position and these alleles in place of its own, and this 3' shift.
    pub fn apply<'b, 'c: 'b>(&'b self, call: &TumorOnlyVariant<'c>) -> TumorOnlyVariant<'b> {
        let call: TumorOnlyVariant<'b> = call.clone();
        TumorOnlyVariant {
            start: self.start,
            end: self.end,
            ref_allele: &self.ref_allele,
            alt_allele: &self.alt_allele,
            num_bases_3_prime_shift_for_deletions: self.indel_3p_shift.unwrap_or(0),
            ..call
        }
    }
}

/// Normalizes a call against the reference, refusing one whose REF allele the reference does not
/// have at its position.
pub fn normalize(reference: &mut Reference, call: &TumorOnlyVariant) -> Result<Normalized, String> {
    let (contig, start) = (call.contig, call.start);
    check_ref(reference, contig, start, call.ref_allele)?;
    let unchanged = Normalized {
        start,
        end: call.end,
        ref_allele: uppercase(call.ref_allele),
        alt_allele: uppercase(call.alt_allele),
        indel_3p_shift: None,
    };
    if call.ref_allele == call.alt_allele
        || !is_sequence(call.ref_allele)
        || !is_sequence(call.alt_allele)
    {
        return Ok(unchanged);
    }
    let ref_allele = unchanged.ref_allele.into_bytes();
    let alt_allele = unchanged.alt_allele.into_bytes();
    let (position, ref_allele, alt_allele) =
        left_align(reference, contig, start - 1, ref_allele, alt_allele)?;
    let ref_allele = String::from_utf8(ref_allele).expect("Alleles are ASCII letters!");
    let alt_allele = String::from_utf8(alt_allele).expect("Alleles are ASCII letters!");
    let indel_3p_shift = match variant_class(&ref_allele, &alt_allele) {
        Some("Insertion" | "Deletion") => Some(three_prime_shift(
            reference,
            contig,
            position,
            &ref_allele,
            &alt_allele,
        )?),
        _ => None,
    };
    Ok(Normalized {
        start: position + 1,
        end: position + ref_allele.len() as u64,
        ref_allele,
        alt_allele,
        indel_3p_shift,
    })
}

/// Refuses a REF allele the reference does not spell at the call's position, ignoring case.
fn check_ref(
    reference: &mut Reference,
    contig: &str,
    start: u64,
    ref_allele: &str,
) -> Result<(), String> {
    if start == 0 {
        return Err(format!(
            "Expected a 1-based position, found {contig}:{start}!"
        ));
    }
    let expected = reference.bases(contig, start - 1, start - 1 + ref_allele.len() as u64)?;
    if !expected.eq_ignore_ascii_case(ref_allele.as_bytes()) {
        return Err(format!(
            "The REF allele '{ref_allele}' at {contig}:{start} is '{}' in the reference: give --reference the FASTA VarDict was run with!",
            String::from_utf8_lossy(expected)
        ));
    }
    Ok(())
}

/// Whether an allele is a plain sequence of bases rather than symbolic.
fn is_sequence(allele: &str) -> bool {
    !allele.is_empty() && allele.bytes().all(|base| base.is_ascii_alphabetic())
}

/// Uppercases an allele's bases, leaving a symbolic allele as written.
fn uppercase(allele: &str) -> String {
    if is_sequence(allele) {
        allele.to_ascii_uppercase()
    } else {
        allele.to_string()
    }
}

/// Trims the shared trailing bases of two alleles at a 0-based position, extending both with the
/// reference base before them whenever one runs out, then trims their shared leading bases while
/// both keep more than one, as `bcftools norm` does; it stops at the contig's first base.
fn left_align(
    reference: &mut Reference,
    contig: &str,
    mut position: u64,
    mut ref_allele: Vec<u8>,
    mut alt_allele: Vec<u8>,
) -> Result<(u64, Vec<u8>, Vec<u8>), String> {
    while ref_allele.last() == alt_allele.last() {
        if ref_allele.len().min(alt_allele.len()) <= 1 && position == 0 {
            break;
        }
        ref_allele.pop();
        alt_allele.pop();
        if ref_allele.is_empty() || alt_allele.is_empty() {
            position -= 1;
            let base = reference.base(contig, position)?;
            ref_allele.insert(0, base);
            alt_allele.insert(0, base);
        }
    }
    let shared = ref_allele
        .iter()
        .zip(&alt_allele)
        .take(ref_allele.len().min(alt_allele.len()) - 1)
        .take_while(|(r, a)| r == a)
        .count();
    ref_allele.drain(..shared);
    alt_allele.drain(..shared);
    Ok((position + shared as u64, ref_allele, alt_allele))
}

/// The bases a left-aligned insertion or deletion at a 0-based position can slide toward the 3'
/// end and still describe the same haplotype, counted as VarDict counts its Shift3: the inserted
/// or deleted bases followed by the reference after them, compared base by base with that reference.
fn three_prime_shift(
    reference: &mut Reference,
    contig: &str,
    position: u64,
    ref_allele: &str,
    alt_allele: &str,
) -> Result<i32, String> {
    let length = reference.length(contig)?;
    let after = position + ref_allele.len() as u64;
    let inserted = alt_allele.as_bytes().get(1..).unwrap_or_default();
    let deleted = ref_allele.len() as u64 - 1;
    let mut shift: u64 = 0;
    while after + shift < length {
        let moved = if deleted > 0 {
            reference.base(contig, after + shift - deleted)?
        } else if let Some(&base) = inserted.get(shift as usize) {
            base
        } else {
            reference.base(contig, after + shift - inserted.len() as u64)?
        };
        if moved != reference.base(contig, after + shift)? {
            break;
        }
        shift += 1;
    }
    Ok(i32::try_from(shift).unwrap_or(i32::MAX))
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::rstest;

    use super::*;
    use crate::reference::tests::fasta;

    /// Normalizes a call on contig `c` of a reference holding `sequence`.
    fn norm(
        sequence: &str,
        start: u64,
        ref_allele: &str,
        alt_allele: &str,
    ) -> Result<Normalized, String> {
        let (_dir, path) = fasta(&[("c", sequence)]);
        let mut reference = Reference::from_path(&path).unwrap();
        let call = TumorOnlyVariant {
            contig: "c",
            start,
            end: start + ref_allele.len() as u64 - 1,
            ref_allele,
            alt_allele,
            ..Default::default()
        };
        normalize(&mut reference, &call)
    }

    /// The 1-based position and alleles of a normalized call.
    fn site(normalized: Normalized) -> (u64, String, String) {
        (
            normalized.start,
            normalized.ref_allele,
            normalized.alt_allele,
        )
    }

    #[rstest]
    #[case::deletion_in_a_homopolymer("GCAAAAT", 5, "AA", "A", (2, "CA", "C"))]
    #[case::insertion_in_a_homopolymer("GCAAAAT", 6, "A", "AA", (2, "C", "CA"))]
    #[case::deletion_in_a_tandem_repeat("GTCACACAG", 6, "ACA", "A", (2, "TCA", "T"))]
    #[case::insertion_in_a_tandem_repeat("GTCACACAG", 7, "C", "CAC", (2, "T", "TCA"))]
    #[case::deletion_of_a_repeat_unit_rotated("GTCACACAG", 3, "CAC", "C", (2, "TCA", "T"))]
    #[case::left_aligned_deletion_is_unchanged("GCAAAAT", 2, "CA", "C", (2, "CA", "C"))]
    #[case::deletion_without_a_left_anchor("GCAAAAT", 6, "AT", "T", (2, "CA", "C"))]
    #[case::insertion_without_a_left_anchor("GCAAAAT", 7, "T", "AT", (2, "C", "CA"))]
    #[case::complex_with_shared_flanks("GCATGAT", 2, "CATG", "CTTG", (3, "A", "T"))]
    #[case::complex_indel_with_a_shared_suffix("GCATGAT", 2, "CATG", "CTGG", (3, "AT", "TG"))]
    #[case::mnv_with_shared_flanks("GCATGAT", 2, "CATGA", "CGCGA", (3, "AT", "GC"))]
    #[case::mnv_trimmed_to_an_snv("GCATGAT", 3, "AT", "CT", (3, "A", "C"))]
    #[case::snv_is_unchanged("GCATGAT", 4, "T", "G", (4, "T", "G"))]
    #[case::deletion_at_the_contig_start("AAAACG", 1, "AA", "A", (1, "AA", "A"))]
    #[case::insertion_at_the_contig_start("AAAACG", 3, "A", "AA", (1, "A", "AA"))]
    #[case::deletion_stops_at_an_n("GNCACAGT", 4, "ACA", "A", (2, "NCA", "N"))]
    #[case::deletion_of_ns("GNNNNT", 3, "NN", "N", (1, "GN", "G"))]
    #[case::soft_masked_padding_is_uppercase("GcaaaaT", 5, "AA", "A", (2, "CA", "C"))]
    #[case::lowercase_alleles_are_uppercased("GCAAAAT", 5, "aa", "a", (2, "CA", "C"))]
    #[case::reference_only_is_unchanged("GCAAAAT", 3, "a", "a", (3, "A", "A"))]
    #[case::symbolic_is_unchanged("GCAAAAT", 2, "C", "<DEL>", (2, "C", "<DEL>"))]
    fn test_normalize(
        #[case] sequence: &str,
        #[case] start: u64,
        #[case] ref_allele: &str,
        #[case] alt_allele: &str,
        #[case] expected: (u64, &str, &str),
    ) {
        let normalized = norm(sequence, start, ref_allele, alt_allele).unwrap();
        let expected = (expected.0, expected.1.to_string(), expected.2.to_string());
        assert_eq!(site(normalized), expected);
    }

    #[rstest]
    #[case::deletion_in_a_homopolymer("GCAAAAT", 5, "AA", "A", Some(3))]
    #[case::insertion_in_a_homopolymer("GCAAAAT", 6, "A", "AA", Some(4))]
    #[case::deletion_in_a_tandem_repeat("GTCACACAG", 6, "ACA", "A", Some(4))]
    #[case::insertion_in_a_tandem_repeat("GTCACACAG", 7, "C", "CAC", Some(6))]
    #[case::deletion_that_cannot_slide("GCATGAT", 2, "CA", "C", Some(0))]
    #[case::deletion_sliding_to_the_contig_end("GCAAAA", 2, "CA", "C", Some(3))]
    #[case::insertion_sliding_to_the_contig_end("GCAAAA", 2, "C", "CA", Some(4))]
    #[case::substitution("GCATGAT", 3, "AT", "CG", None)]
    #[case::complex("GCATGAT", 3, "AT", "CGG", None)]
    fn test_normalize_counts_the_3p_shift(
        #[case] sequence: &str,
        #[case] start: u64,
        #[case] ref_allele: &str,
        #[case] alt_allele: &str,
        #[case] expected: Option<i32>,
    ) {
        let normalized = norm(sequence, start, ref_allele, alt_allele).unwrap();
        assert_eq!(normalized.indel_3p_shift, expected);
    }

    #[test]
    fn test_normalize_keeps_the_end_of_a_symbolic_allele() {
        let (_dir, path) = fasta(&[("c", "GCAAAAT")]);
        let mut reference = Reference::from_path(&path).unwrap();
        let call = TumorOnlyVariant {
            contig: "c",
            start: 2,
            end: 6,
            ref_allele: "C",
            alt_allele: "<DEL>",
            ..Default::default()
        };
        assert_eq!(normalize(&mut reference, &call).unwrap().end, 6);
    }

    #[rstest]
    #[case::wrong_base(
        "GCAAAAT",
        2,
        "G",
        "A",
        "The REF allele 'G' at c:2 is 'C' in the reference: give --reference the FASTA VarDict was run with!"
    )]
    #[case::wrong_symbolic_base(
        "GCAAAAT",
        2,
        "A",
        "<DEL>",
        "The REF allele 'A' at c:2 is 'C' in the reference: give --reference the FASTA VarDict was run with!"
    )]
    #[case::past_the_end(
        "GCAAAAT",
        7,
        "TA",
        "T",
        "Position c:8 is beyond the end of c, which has 7 bases in the reference: give --reference the FASTA VarDict was run with!"
    )]
    fn test_normalize_refuses_a_ref_allele_the_reference_lacks(
        #[case] sequence: &str,
        #[case] start: u64,
        #[case] ref_allele: &str,
        #[case] alt_allele: &str,
        #[case] expected: &str,
    ) {
        assert_eq!(
            norm(sequence, start, ref_allele, alt_allele).unwrap_err(),
            expected
        );
    }

    #[test]
    fn test_apply_puts_the_normalized_site_on_a_call() {
        let call = TumorOnlyVariant {
            contig: "c",
            start: 5,
            end: 6,
            ref_allele: "AA",
            alt_allele: "A",
            alt_depth: 7,
            num_bases_3_prime_shift_for_deletions: 0,
            ..Default::default()
        };
        let normalized = Normalized {
            start: 2,
            end: 3,
            ref_allele: "CA".to_string(),
            alt_allele: "C".to_string(),
            indel_3p_shift: Some(3),
        };
        let applied = normalized.apply(&call);
        assert_eq!(
            (
                applied.start,
                applied.end,
                applied.ref_allele,
                applied.alt_allele
            ),
            (2, 3, "CA", "C")
        );
        assert_eq!(applied.num_bases_3_prime_shift_for_deletions, 3);
        assert_eq!(applied.alt_depth, 7);
    }
}
