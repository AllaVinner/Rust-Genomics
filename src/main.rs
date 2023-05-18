use clap::{Parser, Subcommand};
mod fasta;

use std::fs;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Adds files to myapp
    NumRecords{file_path: String},
    RecordLengths{file_path: String},
}

fn main() {
    let cli = Cli::parse();

    // You can check for the existence of subcommands, and if found use their
    // matches just as you would the top level cmd
    match &cli.command {
        Commands::NumRecords {file_path} => {
            let file = fs::read_to_string(file_path).unwrap().replace("\r", "");
            let num_records = fasta::num_records(&file);
            println!("Num records was {:?}", num_records);
        }
        Commands::RecordLengths {file_path} => {
            let file = fs::read_to_string(file_path).unwrap().replace("\r", "");
            let lengths = fasta::record_lengths(&file);
            println!("Num records was {:?}", lengths);
        }
    }
}