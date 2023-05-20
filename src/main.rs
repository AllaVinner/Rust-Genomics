#![feature(iter_array_chunks)]
#![feature(array_chunks)]
use std::fs;
mod virtual_fasta;


fn main() {
    let file_path = "./data/dna-example.fasta";
    let file = fs::read_to_string(file_path).unwrap().replace("\r", "");
    virtual_fasta::main(&file);
}