use std::str::{Chars, Split};
use std::iter::{FilterMap, Map, Skip};


#[derive(Debug, Copy, Clone, PartialEq)]
pub enum DNABase {
    A,
    C,
    G,
    T
}


fn char_to_dna_base(c: char) -> Option<DNABase> {
    match c {
        'A' | 'a' => Some(DNABase::A),
        'C' | 'c' => Some(DNABase::C),
        'G' | 'g' => Some(DNABase::G),
        'T' | 't' => Some(DNABase::T),
        '\n' => None,
        other => panic!("Got unexpected caracter in raw data. {:?}", other)
    }
}

pub type VirtualSequence<'a> = FilterMap<Chars<'a>, fn(char) -> Option<DNABase>>;
pub struct VirtualRead<'a> {
    id: &'a str,
    desc: & 'a str,
    seq: VirtualSequence<'a>
}

pub fn str_to_virtual_seq<'a>(input: &'a str) -> VirtualSequence<'a> {
    input.chars().filter_map(char_to_dna_base)
}

pub fn str_to_virtual_read<'a>(input: &'a str) -> VirtualRead<'a> {
    let (header, data) = input.split_once("\n").unwrap();
    let (id, desc) = header.split_once(" ").unwrap();
    VirtualRead {
        id: id,
        desc: desc,
        seq: str_to_virtual_seq(data)
    }
}


pub type VirtualFasta<'a> = Map<Skip<Split<'a, &'a str>>, fn(&'a str) -> VirtualRead<'a>>;
pub fn str_to_virtual_fasta<'a>(input: &'a str) -> VirtualFasta<'a> {
    input.split(">").skip(1).map(str_to_virtual_read)
}


pub fn main(fasta: &str) {
    let c = str_to_virtual_fasta(fasta).count();
    println!("Number of records: {c}");
    let (longest_read, read_len) = str_to_virtual_fasta(fasta)
        .map(|read| (read.id, read.seq.count()))
        .reduce(|(max_id, max_len), (id, len)| {
            if len > max_len {
                return (id, len)
            }
            return (max_id, max_len);
        }).unwrap();
    println!("Longest read: {longest_read}, with len {read_len}.");
    
}
