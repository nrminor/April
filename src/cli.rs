use clap::{Parser, Subcommand};

/// Shout out to https://patorjk.com/software/taag/ for the ASCII art.
const INFO: &str = r"

 ░▒▓██████▓▒░  ░▒▓███████▓▒░  ░▒▓███████▓▒░  ░▒▓█▓▒░ ░▒▓█▓▒░        
░▒▓█▓▒░░▒▓█▓▒░ ░▒▓█▓▒░░▒▓█▓▒░ ░▒▓█▓▒░░▒▓█▓▒░ ░▒▓█▓▒░ ░▒▓█▓▒░        
░▒▓█▓▒░░▒▓█▓▒░ ░▒▓█▓▒░░▒▓█▓▒░ ░▒▓█▓▒░░▒▓█▓▒░ ░▒▓█▓▒░ ░▒▓█▓▒░        
░▒▓████████▓▒░ ░▒▓███████▓▒░  ░▒▓███████▓▒░  ░▒▓█▓▒░ ░▒▓█▓▒░        
░▒▓█▓▒░░▒▓█▓▒░ ░▒▓█▓▒░        ░▒▓█▓▒░░▒▓█▓▒░ ░▒▓█▓▒░ ░▒▓█▓▒░        
░▒▓█▓▒░░▒▓█▓▒░ ░▒▓█▓▒░        ░▒▓█▓▒░░▒▓█▓▒░ ░▒▓█▓▒░ ░▒▓█▓▒░        
░▒▓█▓▒░░▒▓█▓▒░ ░▒▓█▓▒░        ░▒▓█▓▒░░▒▓█▓▒░ ░▒▓█▓▒░ ░▒▓████████▓▒░ 


April: Help consensus sequences spring up and stand out
=======================================================

Quality control tool for identifying consensus sequences that over-index
with reference sequence k-mers or primer sequences. April helps noteworthy
genomes spring up and stand out from large sequence datasets, where quality-
control requirements may not guarantee that accessions are without common
bioinformatic artifacts. Specifically, April searchers for either 1) Primer
sequences that were not trimmed from the assembly prior to consensus-calling,
or 2) Reference sequences that may inadvertantly have been inserted into the
consensus sequence while it was being generated.

";

#[derive(Parser)]
#[clap(name = "alpine")]
#[clap(about = INFO)]
#[clap(version = "v0.1.2")]
struct Cli {
    #[command(flatten)]
    verbose: clap_verbosity_flag::Verbosity,

    #[arg(short, long, default_value_t = 3, required = false)]
    threads: u8,

    #[arg(short, long, default_value = "april_report.txt", required = false)]
    report_name: String,

    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    #[clap(
        about = "Find and eliminate FASTA records that are over-enriched with primer sequences."
    )]
    FindPrimers {
        /// FASTA format file of sequences. Reads from stdin if no path is provided.
        #[arg(short, long, required = true)]
        input_fasta: Option<String>,

        /// FASTA of primer sequences. Can be produced from a primer BED file with `bedtools getfasta`.
        #[arg(short, long, required = true)]
        primer_fasta: String,

        /// Output file path. Writes to stdout if no path is provided.
        #[arg(short, long, required = false)]
        output_fasta: Option<String>,
    },

    #[clap(
        about = "Find and eliminate FASTA records that are over-enriched with reference k-mers."
    )]
    FindRef {
        /// FASTA format file of sequences. Reads from stdin if no path is provided.
        #[arg(short, long, required = true)]
        input_fasta: Option<String>,

        /// FASTA file containing the reference sequence(s).
        #[arg(short, long, required = true)]
        ref_fasta: String,

        /// Desired kmer size for the reference sequence scan.
        #[arg(short, long, required = false, default_value_t = 32)]
        kmer_size: i8,

        /// Output file path. Writes to stdout if no path is provided.
        #[arg(short, long, required = false)]
        output_fasta: Option<String>,
    },
}
