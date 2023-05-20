

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum DNABase {
    A,
    C,
    G,
    T
}

type DNASequence = Vec<DNABase>;
pub struct Read {
    id: String,
    desc: String,
    seq: DNASequence
}

pub struct Fasta {
    reads: Vec<Read>
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

pub fn str_to_seq(input: &str) -> DNASequence {
    input.chars().map(char_to_dna_base).collect()
}

pub fn str_to_read(input: str) -> Read{
    let (a, b) =  input.split_once(">").unwrap();
    if a == "" {
        input = b;
    }

    let (header, data) = input.split_once("\n").unwrap();
    let (id, desc) = header.split_once(" ").unwrap();
    Read {
        id: id.to_string(),
        desc: desc.to_string(),
        seq: str_to_seq(data)
    }
}

pub fn main(fasta: &str) {
    let c = str_to_fasta(fasta).len();
    println!("Number of records: {c}");

    let (longest_read, read_len) = str_to_virtual_fasta(fasta)
        .iter()
        .map(|read| (read.id, read.seq.len()))
        .reduce(|(max_id, max_len), (id, len)| {
            if len > max_len {
                return (id, len)
            }
            return (max_id, max_len);
        }).unwrap();
    println!("Longest read: {longest_read}, with len {read_len}.");
    
}





