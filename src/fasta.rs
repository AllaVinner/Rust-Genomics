

const STOP_CODONS: [[DNABase; 3]; 3] = [[DNABase::T, DNABase::A, DNABase::A],
                                        [DNABase::T, DNABase::A, DNABase::G],
                                        [DNABase::T, DNABase::G, DNABase::A]];                  

const START_CODONS: [[DNABase; 3]; 1] = [[DNABase::A, DNABase::T, DNABase::G]];
const CODON_LENGTH: usize = 3;

type DNASequence = Vec<DNABase>;
type ORF<'a> = &'a [DNABase];

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum DNABase {
    A,
    C,
    G,
    T
}

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
    input.chars().filter_map(char_to_dna_base).collect()
}

pub fn str_to_read(input: &str) -> Read{
    let start = match input.split_once(">") {
        None => input,
        Some((a,b)) => if a == "" {
            b
        } else {
            input
        },
    };

    let (header, data) = start.split_once("\n").unwrap();
    let (id, desc) = header.split_once(" ").unwrap();
    Read {
        id: id.to_string(),
        desc: desc.to_string(),
        seq: str_to_seq(data)
    }
}

pub fn str_to_fasta(input: &str) -> Vec<Read> {
    input.split(">").skip(1).map(str_to_read).collect()
}

pub fn seq_to_orfs<'a>(seq: &'a DNASequence) -> Vec<ORF<'a>> {
    seq.array_chunks::<CODON_LENGTH>()
   .enumerate()
   .fold((Vec::<usize>::new(), Vec::<ORF>::new()), 
         |(mut buffer, mut orfs), (codon_i, codon)| {
         
             if START_CODONS.contains(codon) {
                 buffer.push(codon_i);
             }
             
             if STOP_CODONS.contains(codon) {
                 for start_i in &buffer{
                    orfs.push(&seq[CODON_LENGTH*start_i..CODON_LENGTH*(codon_i+1)]);
                 }
                 buffer.clear();
             }
             return (buffer, orfs)
         }
         ).1
}



   
         

pub fn main(fasta: &str) {
    let c = str_to_fasta(fasta).len();
    println!("Number of records: {c}");

    let (longest_read, read_len) = str_to_fasta(fasta)
        .iter()
        .map(|read| (read.id.clone(), read.seq.len()))
        .reduce(|(max_id, max_len), (id, len)| {
            if len > max_len {
                return (id, len)
            }
            return (max_id, max_len);
        }).unwrap();
    println!("Longest read: {longest_read}, with len {read_len}.");

}





