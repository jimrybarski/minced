use minced::parse;
use std::include_str;

fn main() {
    let data = include_str!("example.minced");
    let arrays = parse(data).unwrap();
    for (n, array) in arrays.iter().enumerate() {
        let count = array.repeat_spacers.len();
        println!("Array {n} has {count} repeat spacers");
    }
}
