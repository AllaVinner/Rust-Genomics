use std::cmp::min;



const STOP_CODONS: [[&DNABase; 3]; 3] = [[&DNABase::T, &DNABase::A, &DNABase::A],
                                        [&DNABase::T, &DNABase::A, &DNABase::G],
                                        [&DNABase::T, &DNABase::G, &DNABase::A]];                  

const START_CODONS: [[&DNABase; 3]; 1] = [[&DNABase::A, &DNABase::T, &DNABase::G]];
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

pub fn seq_to_orfs<'a>(seq: &'a DNASequence, offset: usize) -> Vec<ORF<'a>> {
    seq.iter()
        .skip(offset)
        .array_chunks::<CODON_LENGTH>()
        .enumerate()
        .fold((Vec::<usize>::new(), Vec::<ORF>::new()), 
            |(mut buffer, mut orfs), (codon_i, codon)| {
            
                if START_CODONS.contains(&codon) {
                    buffer.push(codon_i);
                }
                
                if STOP_CODONS.contains(&codon) {
                    for start_i in &buffer{
                        orfs.push(&seq[CODON_LENGTH*start_i..CODON_LENGTH*(codon_i+1)]);
                    }
                    buffer.clear();
                }
                return (buffer, orfs)
            }
        ).1
}


fn find_n_mers(seq: &DNASequence, n: usize) -> Vec<Vec<usize>> {
    let mut n_mers_indices: Vec<Vec<usize>> = Vec::new();
    for base_i in 0..(seq.len()-n+1) {
        let base_slice = &seq[base_i..(base_i+n)];
        let is_base_slice_checked = n_mers_indices.iter()
                           .map(|repeats| repeats.get(0).unwrap())
                           .map(|&start_i| &seq[start_i..(start_i+n)])
                           .any(|compare_slice| base_slice == compare_slice);
                           
        if is_base_slice_checked {
            continue;
        }
        
        let mut base_repeats = Vec::new();
        for compare_i in base_i..(seq.len()-n+1) {
            if base_slice == &seq[compare_i..(compare_i+n)] {
                base_repeats.push(compare_i);
            }
        }
         
        if base_repeats.len() > 1 {
            n_mers_indices.push(base_repeats)
        }
    }  
    return n_mers_indices 
}


fn max_of_tupple<A, B>(tup1: (A, B), tup2: (A, B)) -> (A, B) 
where B: PartialOrd {
    if tup2.1 > tup1.1 {
        return tup2
    }
    return tup1;
}   

fn min_of_tupple<A, B>(tup1: (A, B), tup2: (A, B)) -> (A, B) 
where B: PartialOrd {
    if tup2.1 < tup1.1 {
        return tup2
    }
    return tup1;
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
    
    let (shortest_read, shortest_read_len) = str_to_fasta(fasta)
        .iter()
        .map(|read| (read.id.clone(), read.seq.len()))
        .reduce(min_of_tupple)
        .unwrap();
    println!("Shortest read: {shortest_read}, with len {shortest_read_len}.");


    let (longest_ord_read, orf_len) = str_to_fasta(fasta)
        .iter()
        .map(|read| (read.id.clone(), seq_to_orfs(&read.seq,0)))
        .map(|(id, orfs)| (id, orfs.iter().map(|orf| orf.len()).max().unwrap_or(0)))
        .reduce(|(max_id, max_len), (id, len)| {
            if len > max_len {
                return (id, len)
            }
            return (max_id, max_len);
        }).unwrap();
    println!("Longest orf read: {longest_ord_read}, with len {orf_len}.");

    let (longest_ord_read, orf_len) = str_to_fasta(fasta)
    .iter()
    .filter(|read| read.id == "gi|142022655|gb|EQ086233.1|16")
    .map(|read| (read.id.clone(), seq_to_orfs(&read.seq,2)))
    .map(|(id, orfs)| (id, orfs.iter().map(|orf| orf.len()).max().unwrap_or(0)))
    .reduce(|(max_id, max_len), (id, len)| {
        if len > max_len {
            return (id, len)
        }
        return (max_id, max_len);
    }).unwrap();
    println!("Longest orf read in gi|142022655|gb|EQ086233.1|16: {longest_ord_read}, with len {orf_len}.");  
    
    let ((longest_nmer_read, seq_start), nmer_len) = str_to_fasta(fasta)
        .iter()
        .map(|read| (read.id.clone(), find_n_mers(&read.seq, 10)))
        .map(|(id, nmers)| (id, nmers.iter().map(|nmer| (*nmer.get(0).unwrap(), nmer.len()))
                                                                     .reduce(max_of_tupple)
                                                                     .unwrap_or((0, 0))
                                                                    ))
        .map(|(id, (start_i, max_len))| ((id,start_i), max_len))
        .reduce(max_of_tupple)
        .unwrap();
    println!("Longest nmer read: {longest_nmer_read}, with len {nmer_len},  starting at {seq_start}.");
    


}





