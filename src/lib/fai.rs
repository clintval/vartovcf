//! A module for serialization-deserialization of FASTA FAI index records.
use std::error;
use std::fmt::Debug;
use std::path::Path;

use anyhow::Result;
use csv::ReaderBuilder;
use rust_htslib::bcf::Header;
use serde::Deserialize;

/// A FASTA FAI index record.
#[derive(Debug, Deserialize)]
pub struct FaiRecord<'a> {
    /// Name of this reference sequence.
    pub name: &'a str,
    /// Total length of this reference sequence, in bases.
    pub length: u64,
    /// Offset within the FASTA file of this sequence's first base.
    pub offset: u64,
    /// The number of bases on each line.
    pub linebases: u64,
    /// The number of bytes in each line, including the newline.
    pub linewidth: u64,
}

impl<'a> FaiRecord<'a> {
    /// Convert this FAI record into a contig header line for a VCF file.
    fn to_vcf_contig_record(&self) -> String {
        format!("##contig=<ID={},length={}>", self.name, self.length)
    }
}

/// Read a FASTA FAI file and create a list of VCF contig header records.
pub fn vcf_contig_header_records<I>(fai: I) -> Result<Vec<String>, Box<dyn error::Error>>
where
    I: AsRef<Path> + Debug,
{
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(fai)?; // Find FAI file based on reference file location.

    let mut carry = csv::StringRecord::new();
    let mut records: Vec<String> = Vec::new();

    while reader.read_record(&mut carry)? {
        let rec: FaiRecord = carry.deserialize(None)?;
        records.push(rec.to_vcf_contig_record());
    }

    Ok(records)
}

/// Add a collection of VCF contig header records to the supplied VCF header.
pub fn contigs_to_vcf_header(contigs: &[String], mut header: Header) -> Header {
    for contig in contigs.iter() {
        header.push_record(contig.as_bytes());
    }
    header
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use anyhow::Result;
    use tempfile::NamedTempFile;

    use super::*;

    #[test]
    fn test_fai_to_vcf_contig_record() {
        let rec = FaiRecord {
            name: "chr1",
            length: 248956422,
            offset: 112,
            linebases: 70,
            linewidth: 71,
        };
        let expected = "##contig=<ID=chr1,length=248956422>";
        assert_eq!(rec.to_vcf_contig_record(), expected);
    }

    #[test]
    fn test_vcf_contig_header_records() -> Result<(), Box<dyn std::error::Error>> {
        let mut file = NamedTempFile::new()?;
        writeln!(file, "chr1\t248956422\t112\t70\t71")?;
        writeln!(file, "chr2\t242193529\t252513167\t70\t71")?;
        writeln!(file, "chr3\t198295559\t498166716\t70\t71")?;
        writeln!(file, "chr4\t190214555\t699295181\t70\t71")?;
        writeln!(file, "chr5\t181538259\t892227221\t70\t71")?;
        let actual = vcf_contig_header_records(file.path())?;
        let expected = vec![
            "##contig=<ID=chr1,length=248956422>",
            "##contig=<ID=chr2,length=242193529>",
            "##contig=<ID=chr3,length=198295559>",
            "##contig=<ID=chr4,length=190214555>",
            "##contig=<ID=chr5,length=181538259>",
        ];
        for (left, right) in actual.iter().zip(expected.iter()) {
            assert_eq!(&left, &right)
        }
        Ok(())
    }
}
