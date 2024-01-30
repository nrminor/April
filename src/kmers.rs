use std::{borrow::Cow, path::Path, rc::Rc};

use anyhow::Result;
use needletail::{parse_fastx_file, FastxReader, Sequence};

/// Read the reference sequence in a separate function in preparationg
/// for a future ansyc implementation
fn buffer_ref(target_ref: &Path) -> Result<Box<dyn FastxReader>> {
    let ref_reader = parse_fastx_file(target_ref)?;

    Ok(ref_reader)
}

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

///TODO
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
