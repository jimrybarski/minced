use minced::parse;
use std::fs::File;
use std::io::{BufReader, Read};

fn main() {
    let file = File::open("examples/minced.txt").unwrap();
    let mut reader = BufReader::new(file);
    let mut input = String::new();
    reader.read_to_string(&mut input).unwrap();
    let contigs = parse(&input).unwrap();
    for contig in contigs {
        println!("{} has {} arrays", contig.accession, contig.arrays.len());
    }
}
