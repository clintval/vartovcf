//! A module for serialization-deserialization of FASTA FAI index records.
use std::error;
use std::fmt::Debug;
use std::path::Path;

use csv::ReaderBuilder;
use rust_htslib::bcf::Header;
use serde::Deserialize;

/// A FASTA FAI index record.
#[derive(Debug, Deserialize)]
pub struct FaiRecord<'a> {
    /// Name of this reference sequence.
    pub name: &'a str,
    /// Total length of this reference sequence, in bases.
    pub length: &'a str,
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
pub fn vcf_contig_header_records<I>(input: I) -> Result<Vec<String>, Box<dyn error::Error>>
where
    I: AsRef<Path> + Debug,
{
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(input)?; // Find FAI file based on reference file location.

    let mut carry = csv::StringRecord::new();
    let mut records: Vec<String> = Vec::new();
    while reader.read_record(&mut carry)? {
        let fai: FaiRecord = carry.deserialize(None)?;
        records.push(fai.to_vcf_contig_record());
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
