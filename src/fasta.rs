




const STOP_CODONS: [[&DNABase; 3]; 3] = [[&DNABase::T, &DNABase::A, &DNABase::A],
                                        [&DNABase::T, &DNABase::A, &DNABase::G],
                                        [&DNABase::T, &DNABase::G, &DNABase::A]];
                     

const START_CODONS: [[&DNABase; 3]; 1] = [[&DNABase::A, &DNABase::T, &DNABase::G]];

const CODON_LENGTH: usize = 3;

//TODO: Move this at some point
#[derive(Debug, Copy, Clone, PartialEq)]
pub enum DNABase {
    A,
    C,
    G,
    T
}

#[derive(Debug)]
pub struct FastaSequence {
    identifier: String,
    description: String,
    data: Vec<DNABase>
}


pub fn parse_fasta_sequence(input: &str) -> FastaSequence{

    // Assumptions
    // Has one new line
    let (header, str_data) = input.split_once("\n").unwrap();
    let (identifier, description) = if header.contains(" ") {
        header.split_once(" ").unwrap()
    } else {
        (header, "")
    };

    
    let mut data: Vec<DNABase> = str_data
        .trim_end()
        .chars().filter_map(|c| match c {
            'A' | 'a' => Some(DNABase::A),
            'C' | 'c' => Some(DNABase::C),
            'G' | 'g' => Some(DNABase::G),
            'T' | 't' => Some(DNABase::T),
            '\n' => None,
            _ => panic!("Data contains non actg data."),
        }).collect();

    FastaSequence{
        identifier: identifier.to_string(),
        description: description.to_string(),
        data: data
    }
}

pub fn parse_fasta(fasta_input: &str) -> Vec<FastaSequence>{
    fasta_input.split(">")
        .skip_while(|seq| seq == &"")
        .map(|seq| parse_fasta_sequence(seq))
        .collect()
}

fn get_orfs(seq: &Vec<DNABase>, offset: usize) -> Vec<&[DNABase]> { 
    let (_, orfs) = seq
        .iter()
        .skip(offset)
        .array_chunks::<CODON_LENGTH>()
        .enumerate()
        .fold((Vec::<usize>::new(), Vec::<&[DNABase]>::new()), 
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
            );
    return orfs
}




pub fn num_records(fasta_input: &str) -> i32 {
    fasta_input.split(">").count() as i32 - 1
}

pub fn record_lengths(fasta_input: &str) -> Vec<(String, usize)> {
    parse_fasta(fasta_input).iter().map(|f| (f.identifier.clone(), f.data.len())).collect()
}


pub fn record_max_orf(fasta_input: &str) -> Vec<(String, usize)> {
    parse_fasta(fasta_input).iter().map(|f| (f.identifier.clone(), f.data.len())).collect()
}






