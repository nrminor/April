pub mod cli;
pub mod intervals;
pub mod kmers;
pub mod reference;

use anyhow::Result;
use std::{path::PathBuf, rc::Rc};

fn main() -> Result<()> {
    //  variables that will be outsourced to the command line interface
    let filename = PathBuf::from("test_dataset/PP141354.1.fasta");
    let ref_kmers: Rc<[&[u8; 32]; 1]> = Rc::from([b"AAGAAATGCTGGACAACAGGGCAACCTTACAA"]);
    let _amplicon_interval = 300;

    // TODO: these are simply demo calls. They will be replace with clap
    // subcommands soon
    reference::score_all_records(&filename, ref_kmers.clone())?;
    intervals::compute_match_intervals(&filename, ref_kmers.clone())?;

    Ok(())
}
