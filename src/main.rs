use anyhow::Result;
use needletail::{parse_fastx_file, Sequence};
use std::rc::Rc;

fn main() -> Result<()> {
    //  variables that will be outsourced to the command line interface
    let filename = "sandbox/PP141354.1.fasta";
    let test_kmer = Rc::from(b"AAAATTAATTTT");
    let kmer_len = 32;

    let mut reader = parse_fastx_file(filename)?;
    while let Some(record) = reader.next() {
        let seqrec = record?;

        let _record_len = Rc::from(seqrec.num_bases());

        let norm_seq = seqrec.normalize(true);

        let rc = norm_seq.reverse_complement();

        for (i, kmer, _) in norm_seq.canonical_kmers(kmer_len, &rc) {
            if kmer
                .windows((*test_kmer).len())
                .any(|window| window == (*test_kmer))
            {
                println!("Kmer at position {} contains the test kmer:", i);
                println!("{}", std::str::from_utf8(kmer)?);
            }
        }
    }

    Ok(())
}
