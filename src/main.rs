use anyhow::Result;
use needletail::{parse_fastx_file, Sequence};
use std::rc::Rc;

fn main() -> Result<()> {
    //  variables that will be outsourced to the command line interface
    let filename = "test_dataset/PP141354.1.fasta";
    let test_kmer = Rc::from(b"AAAATTAATTTT");
    let kmer_len = 32;
    let _amplicon_interval = 300;

    // open fasta reader
    let mut reader = parse_fastx_file(filename)?;
    while let Some(record) = reader.next() {
        // unwrap the record
        let seqrec = record?;

        // determine the sequence length (potentially for use later?)
        let _record_len = Rc::from(seqrec.num_bases());

        // double check that the sequence is uncorrupted
        let norm_seq = seqrec.normalize(true);

        // retrieve the reverse complement of the fasta record
        let rc = norm_seq.reverse_complement();

        // make a vector to contain kmer "hits" and their locations
        let mut hits: Vec<&[u8]> = Vec::new();
        let mut hit_locs: Vec<usize> = Vec::new();

        // figure out which kmers are hits
        for (pos, kmer, _) in norm_seq.canonical_kmers(kmer_len, &rc) {
            if kmer
                .windows((*test_kmer).len())
                .any(|window| window == (*test_kmer))
            {
                println!("Kmer at position {} contains the test kmer:", pos);
                println!("{}", std::str::from_utf8(kmer)?);
                hits.push(kmer);
                hit_locs.push(pos);
            }
        }

        // count the hits
        let hit_count = hits.len();
        println!("{}", hit_count);
        println!("{:?}", hit_locs);

        // get hit intervals
        let mut _hit_intervals: Vec<usize> = Vec::with_capacity(hits.len());
    }

    Ok(())
}
