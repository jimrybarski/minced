# minced-parser
A Rust parser for the [MinCED](https://github.com/ctSkennerton/minced) CRISPR array annotation tool.

### Installation

Add the following to Cargo.toml:

`minced-parser = 1.0.0`

### Usage

```rust
use std::fs::File;
use std::io::{BufReader, Read};

fn main() {
    let file = File::open("minced.txt").unwrap();
    let mut reader = BufReader::new(file);
    let mut input = String::new();
    reader.read_to_string(&mut input).unwrap();
    let contigs = minced_parser::parse(&input).unwrap();
    for contig in contigs {
        println!("{} has {} arrays", contig.accession, contig.arrays.len());
    }
}
```
