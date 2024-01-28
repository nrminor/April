use anyhow::Result;
use needletail::{parse_fastx_file, Sequence};
use std::path::PathBuf;
use std::rc::Rc;

pub fn compute_match_intervals(input_path: &PathBuf, ref_kmers: Rc<[&[u8; 32]; 1]>) -> Result<()> {
    // TODO: Replace with actual slice of all minimizer kmers from reference
    let test_kmer = *ref_kmers.first().unwrap();

    // open fasta reader
    let mut reader = parse_fastx_file(input_path)?;
    while let Some(record) = reader.next() {
        // unwrap the record
        let seqrec = record?;
        let kmer_len = test_kmer.len() as u8;

        // determine the sequence length (potentially for use later?)
        let record_len = Rc::from(seqrec.num_bases());
        let kmer_space = (*record_len / (kmer_len as usize)) * 2;
        println!("{} possible kmers for the provided reference", kmer_space);

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
        for (i, (pos, kmer, _)) in norm_seq.canonical_kmers(kmer_len, &rc).enumerate() {
            counter = i;
            if kmer
                .windows((*test_kmer).len())
                .any(|window| window == *test_kmer)
            {
                println!("Kmer at position {} contains the test kmer:", pos);
                println!("{}", std::str::from_utf8(kmer)?);
                hits.push(kmer);
                hit_locs.push(pos);
            }
        }

        println!(
            "{} of {} kmers matched reference kmers.",
            hit_locs.len(),
            counter
        );

        // create a companion vector of each previous kmer's locations
        let mut prev_locs = vec![0_usize];
        prev_locs.append(&mut hit_locs.clone());
        prev_locs = prev_locs[0..hit_locs.len()].to_vec();

        // compute intervals between kmer hits
        let mut intervals = Vec::new();
        for (current, previous) in hit_locs.iter().zip(prev_locs.iter()) {
            if previous == &0_usize {
                continue;
            }
            let interval = current - previous;
            intervals.push(interval)
        }

        println!("Interval-ing successful:\n{:?}", intervals)
    }

    Ok(())
}
