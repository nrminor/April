use anyhow::Result;
use needletail::{parse_fastx_file, parser::SequenceRecord, Sequence};
use std::borrow::Cow;
use std::path::Path;
use std::rc::Rc;
use textdistance::str::hamming;

use crate::kmers::ref_to_kmers;

/// .
///
/// # Errors
///
/// This function will return an error if .
fn seq_to_windows<'a>(normalized: &'a Cow<'a, [u8]>, kmer_len: Rc<u8>) -> Result<Vec<&'a [u8]>> {
    let windows: Vec<&[u8]> = normalized.chunks(*kmer_len as usize).collect();

    match windows.is_empty() {
        true => Err(anyhow::Error::msg(
            "Unable to compute k-mers for the provided reference FASTA.",
        )),
        false => Ok(windows),
    }
}

/// .
///
/// # Errors
///
/// This function will return an error if .
fn window_min_score(window: &[u8], ref_kmers: &[Rc<[u8]>]) -> Result<Option<usize>> {
    let scoring_attempt: Result<Vec<usize>> = ref_kmers
        .iter()
        .map(|kmer| {
            let kmer_str = std::str::from_utf8(kmer)?;
            let window_str = std::str::from_utf8(window)?;
            Ok(hamming(kmer_str, window_str))
        })
        .collect();

    match scoring_attempt {
        Ok(dist_scores) => Ok(dist_scores.into_iter().min()),
        Err(message) => Err(anyhow::Error::msg(format!(
            "Scoring windows with reference k-mers failed with the following message:\n{:?}",
            message
        ))),
    }
}

/// .
///
/// # Errors
///
/// This function will return an error if .
fn score_infiltration(
    seqrec: &SequenceRecord<'_>,
    kmer_len: Rc<u8>,
    ref_kmers: &[Rc<[u8]>],
) -> Result<f64> {
    // determine the sequence length (potentially for use later?)
    let record_len = Rc::from(seqrec.num_bases());
    eprintln!("{:?}", record_len);

    // double check that the sequence is uncorrupted
    let norm_seq = seqrec.normalize(true);

    // generate windows for comparison with ref kmers
    let windows = seq_to_windows(&norm_seq, kmer_len.clone())?;
    eprintln!(
        "Number of windows for record of length {} with a K of {}: {}",
        record_len,
        *kmer_len,
        windows.len()
    );

    // collect minimum distances per window
    let scores: Vec<usize> = windows
        .iter()
        .map(|window| {
            let scoring_attempt = window_min_score(window, ref_kmers);
            match scoring_attempt {
                Ok(Some(score)) => score,
                _ => 0_usize,
            }
        })
        .collect();

    // compute a stand-in score that will be replaced with a statistic
    let nil_count = &scores.iter().filter(|score| score == &&0_usize).count();
    let score = (*nil_count as f64) / (scores.len() as f64);

    Ok(score)
}

/// .
///
/// # Panics
///
/// Panics if .
///
/// # Errors
///
/// This function will return an error if .
pub fn score_all_records(input_path: &Path, ref_path: &Path, kmer_len: Rc<u8>) -> Result<()> {
    let ref_kmers = ref_to_kmers(ref_path, kmer_len.clone())?;

    // open fasta reader
    let mut reader = parse_fastx_file(input_path)?;
    while let Some(record) = reader.next() {
        // unwrap the record
        let seqrec = record?;

        let scoring_attempt = score_infiltration(&seqrec, kmer_len.clone(), &ref_kmers);
        match scoring_attempt {
            Ok(score) => {
                eprintln!(
                    "Sequence record {:?} has reference infiltration score of {}",
                    String::from_utf8(seqrec.id().to_vec()).unwrap(),
                    score,
                );
            }
            Err(message) => eprintln!(
                "Scoring record {:?} for reference infiltration failed:\n{:?}",
                String::from_utf8(seqrec.id().to_vec()).unwrap(),
                message
            ),
        }
    }

    Ok(())
}
