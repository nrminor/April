use std::{borrow::Cow, path::Path, rc::Rc};

use anyhow::Result;
use needletail::{parse_fastx_file, FastxReader, Sequence};

pub enum KmerOptions {
    K32 = 32,
    K52 = 52,
}

/// Read the reference sequence in a separate function in preparationg
/// for a future ansyc implementation
fn buffer_ref(target_ref: &Path) -> Result<Box<dyn FastxReader>> {
    let ref_reader = parse_fastx_file(target_ref)?;

    Ok(ref_reader)
}

/// Collect all non-complement k-mers of length `kmer_len` for the
/// reference records in the provided FASTA. Takes a FASTA record
/// with normalized sequences (i.e., with bases 'A', 'T', 'G', or 'C'),
/// along with their reverse complements, and returns a vector of UTF-
/// 8 characters for each base. In the future this will be replaced with
/// bit-kmers for lower memory usage.
fn collect_kmers<'a>(
    normalized: &'a Cow<'a, [u8]>,
    reverse_comp: &'a [u8],
    kmer_len: Rc<u8>,
) -> Result<Vec<&'a [u8]>> {
    let kmer_db: Vec<&[u8]> = normalized
        .canonical_kmers(*kmer_len, reverse_comp)
        .filter_map(
            |(_, kmer, whether_comp)| {
                if whether_comp {
                    None
                } else {
                    Some(kmer)
                }
            },
        )
        .collect();

    match kmer_db.is_empty() {
        true => Err(anyhow::Error::msg(
            "Unable to compute k-mers for the provided reference FASTA.",
        )),
        false => Ok(kmer_db),
    }
}

/// Oversees the conversion of all records in the reference FASTA pointed to
/// with `target_ref` into k-mers of length `kmer_len`. Returns a Result type
/// containing a database of all k-mers with the same orientation as the reference.
pub fn ref_to_kmers(target_ref: &Path, kmer_len: Rc<u8>) -> Result<Vec<Rc<[u8]>>> {
    let mut ref_reader: Box<dyn FastxReader> = buffer_ref(target_ref)?;
    let mut kmer_db: Vec<Rc<[u8]>> = Vec::new();

    while let Some(record) = ref_reader.next() {
        let fasta_entry = record?;
        let normalized = fasta_entry.normalize(true);
        let reverse_comp = fasta_entry.reverse_complement();
        let record_kmers = collect_kmers(&normalized, &reverse_comp, kmer_len.clone())?;

        kmer_db.extend(
            record_kmers
                .into_iter()
                .map(|kmer| Rc::from(kmer.to_vec().into_boxed_slice())),
        );
    }

    Ok(kmer_db)
}
