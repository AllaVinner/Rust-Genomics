#![feature(iter_array_chunks)]
#![feature(array_chunks)]
use std::fs;
mod basic_bio;


fn main() {
    let file_path = "./data/dna2.fasta";
    let file = fs::read_to_string(file_path).unwrap().replace("\r", "");
    basic_bio::main(&file);
}