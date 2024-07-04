//! A module for serialization-deserialization of FASTA FAI index records.
use std::default::Default;
use std::error;
use std::ffi::OsString;
use std::fmt::Debug;
use std::path::{Path, PathBuf};

use anyhow::Result;
use csv::ReaderBuilder;
use rust_htslib::bcf::Header;
use serde::{Deserialize, Serialize};

/// A FASTA FAI index record.
#[derive(Debug, Default, Deserialize, Serialize)]
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

/// Given a path to a FASTA file, resolve its sibling FAI index file.
pub fn fai_file<I: AsRef<Path> + Debug>(fasta: I) -> PathBuf {
    let mut fasta = OsString::from(fasta.as_ref().to_path_buf());
    fasta.push(".fai");
    PathBuf::from(fasta)
}

/// Add a collection of VCF contig header records to the supplied VCF header.
pub fn contigs_to_vcf_header(contigs: &[String], header: &mut Header) {
    for contig in contigs.iter() {
        header.push_record(contig.as_bytes());
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
        .from_path(&fai)
        .unwrap_or_else(|_| panic!("Could not open an FAI reader for file path: {:?}", &fai));

    let mut carry = csv::StringRecord::new();
    let mut records: Vec<String> = Vec::new();

    while reader
        .read_record(&mut carry)
        .expect("Failed to read the FAI record!")
    {
        let rec: FaiRecord = carry
            .deserialize(None)
            .expect("Failed to deserialize the FAI record!");
        records.push(rec.to_vcf_contig_record());
    }

    Ok(records)
}

/// Add the reference contigs from an indexed FASTA file to the VCF header.
pub fn fasta_contigs_to_vcf_header<I>(fasta: I, header: &mut Header)
where
    I: AsRef<Path> + Debug,
{
    let fai = fai_file(&fasta);
    let contigs = vcf_contig_header_records(fai).expect("Could not read the FAI index records!");
    contigs_to_vcf_header(&contigs, header);
}

/// Add a VCF header record containing the filepath of the reference FASTA.
pub fn fasta_path_to_vcf_header<I>(
    fasta: I,
    header: &mut Header,
) -> Result<(), Box<dyn error::Error>>
where
    I: AsRef<Path> + Debug,
{
    match fasta.as_ref().to_str() {
        Some(path) => header.push_record(format!("##reference={}", path).as_bytes()),
        None => return Err("Could not create a string out of the FASTA path".into()),
    };
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use anyhow::Result;
    use pretty_assertions::assert_eq;
    use rstest::*;
    use rust_htslib::bcf::{Format, HeaderRecord, Read};
    use rust_htslib::bcf::{Reader as VcfReader, Writer as VcfWriter};
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
    fn test_fai_file() {
        let fasta = PathBuf::from("/references/reference.fasta");
        assert_eq!(
            fai_file(fasta),
            PathBuf::from("/references/reference.fasta.fai")
        );
        let fasta = PathBuf::from("/references/reference.fa");
        assert_eq!(
            fai_file(fasta),
            PathBuf::from("/references/reference.fa.fai")
        );
    }

    #[fixture]
    fn fai() -> NamedTempFile {
        let file = NamedTempFile::new().expect("Cannot create temporary file!");
        writeln!(&file, "chr1\t248956422\t112\t70\t71").expect("Could not write bytes!");
        writeln!(&file, "chr2\t242193529\t252513167\t70\t71").expect("Could not write bytes!");
        writeln!(&file, "chr3\t198295559\t498166716\t70\t71").expect("Could not write bytes!");
        writeln!(&file, "chr4\t190214555\t699295181\t70\t71").expect("Could not write bytes!");
        writeln!(&file, "chr5\t181538259\t892227221\t70\t71").expect("Could not write bytes!");
        file
    }

    #[fixture]
    fn contig_header_lines() -> Vec<String> {
        vec![
            "##contig=<ID=chr1,length=248956422>".into(),
            "##contig=<ID=chr2,length=242193529>".into(),
            "##contig=<ID=chr3,length=198295559>".into(),
            "##contig=<ID=chr4,length=190214555>".into(),
            "##contig=<ID=chr5,length=181538259>".into(),
        ]
    }

    #[rstest]
    fn test_vcf_contig_header_records(
        fai: NamedTempFile,
        contig_header_lines: Vec<String>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let actual = vcf_contig_header_records(&fai.path())
            .expect("Could not parse contig header records from file!");
        for (left, right) in actual.iter().zip(contig_header_lines.iter()) {
            assert_eq!(&left, &right)
        }
        Ok(())
    }

    #[rstest]
    fn test_contigs_to_vcf_header(
        fai: NamedTempFile,
        contig_header_lines: Vec<String>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let actual = vcf_contig_header_records(&fai.path())
            .expect("Could not parse contig header records from file!");

        let mut header = Header::default();
        contigs_to_vcf_header(&actual, &mut header);

        let file = NamedTempFile::new().expect("Cannot create temporary file!");
        let _ = VcfWriter::from_path(&file.path(), &mut header, true, Format::Vcf).unwrap();
        let reader = VcfReader::from_path(&file.path()).expect("Error opening tempfile!");
        let records = reader.header().header_records();

        fn header_record_matches_contig(record: &HeaderRecord, line: &str) -> bool {
            match record {
                HeaderRecord::Contig { key, values } => {
                    line == &format!(
                        "##{}=<ID={},length={}>",
                        key, values["ID"], values["length"]
                    )
                }
                _ => false,
            }
        }

        let value = contig_header_lines.iter().all(|contig| {
            records
                .iter()
                .any(|record| header_record_matches_contig(record, contig))
        });
        assert!(value);
        Ok(())
    }

    #[test]
    fn test_fasta_path_to_vcf_header_exists() {
        let mut header = Header::default();
        fasta_path_to_vcf_header(&"/references/hg19.fa", &mut header)
            .expect("Could not add the FASTA file path to the VCF header!");

        let file = NamedTempFile::new().expect("Cannot create temporary file!");
        let _ = VcfWriter::from_path(&file.path(), &mut header, true, Format::Vcf).unwrap();
        let reader = VcfReader::from_path(&file.path()).expect("Error opening tempfile!");
        let records = reader.header().header_records();

        let test = records.iter().any(|rec| match rec {
            HeaderRecord::Generic { key, value } => {
                key == "reference" && value == "/references/hg19.fa"
            }
            _ => false,
        });
        assert!(
            test,
            "Could not find the reference record in the VCF header"
        )
    }
}
