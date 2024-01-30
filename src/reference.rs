use anyhow::Result;
use needletail::{parse_fastx_file, parser::SequenceRecord, sequence::minimizer, Sequence};
use std::path::Path;
use std::rc::Rc;

use crate::kmers::ref_to_kmers;

/// TODO
fn score_infiltration(
    seqrec: &SequenceRecord<'_>,
    kmer_len: Rc<u8>,
    test_kmer: &[u8; 32],
) -> Result<usize> {
    // determine the sequence length (potentially for use later?)
    let record_len = Rc::from(seqrec.num_bases());
    let kmer_space = (*record_len / (*kmer_len as usize)) * 2;
    eprintln!("{} possible kmers for the provided reference", kmer_space);

    // double check that the sequence is uncorrupted
    let norm_seq = seqrec.normalize(true);

    // retrieve the reverse complement of the fasta record
    let rc = norm_seq.reverse_complement();

    // make vectors to contain kmer "hits" and their locations, along with
    // a counter for all kmers
    let mut hits: Vec<&[u8]> = Vec::new();
    let mut hit_locs: Vec<usize> = Vec::new();
    let mut counter = 0;

    // figure out which kmers are hits
    for (i, (pos, kmer, _)) in norm_seq.canonical_kmers(*kmer_len, &rc).enumerate() {
        counter = i;
        let _mini = minimizer(kmer, kmer.len());
        if kmer
            .windows((*test_kmer).len())
            .any(|window| window == *test_kmer)
        {
            eprintln!("Kmer at position {} contains the test kmer:", pos);
            eprintln!("{}", std::str::from_utf8(kmer)?);
            hits.push(kmer);
            hit_locs.push(pos);
        }
    }

    eprintln!(
        "{} of {} kmers matched reference kmers.",
        hit_locs.len(),
        counter
    );

    let score = hits.len() / counter;

    Ok(score)
}

/// TODO
pub fn score_all_records(input_path: &Path, ref_kmers: Rc<[&[u8; 32]; 1]>) -> Result<()> {
    // TODO: Replace with actual slice of all minimizer kmers from reference
    let test_kmer: &[u8; 32] = ref_kmers.first().unwrap();
    let kmer_len = Rc::from(test_kmer.len() as u8);
    let _ref_kmers = ref_to_kmers(Path::new("test_dataset/PP141354.1.fasta"), kmer_len.clone());

    // open fasta reader
    let mut reader = parse_fastx_file(input_path)?;
    while let Some(record) = reader.next() {
        // unwrap the record
        let seqrec = record?;

        let scoring_attempt = score_infiltration(&seqrec, kmer_len.clone(), test_kmer);
        match scoring_attempt {
            Ok(score) => {
                eprintln!(
                    "Sequence record {:?} has reference infiltration score of {}",
                    &seqrec.id(),
                    score,
                );
            }
            Err(message) => eprintln!(
                "Scoring record {:?} for reference infiltration failed:\n{:?}",
                &seqrec.id(),
                message
            ),
        }
    }

    Ok(())
}
