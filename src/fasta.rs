//TODO: Move this at some point
#[derive(Debug, Copy, Clone)]
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

pub fn num_records(fasta_input: &str) -> i32 {
    fasta_input.split(">").count() as i32 - 1
}

pub fn record_lengths(fasta_input: &str) -> Vec<(String, usize)> {
    parse_fasta(fasta_input).iter().map(|f| (f.identifier.clone(), f.data.len())).collect()
}




