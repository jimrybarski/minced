#![deny(warnings, missing_docs)]
//! Parses the output produced by MinCED (<https://github.com/ctSkennerton/minced>), a CRISPR array
//! annotation tool.
//!
//! ## Example
//!
//! ```rust
//! use minced_parser::parse;
//! use std::fs::File;
//! use std::io::{BufReader, Read};
//!
//! let file = File::open("examples/minced.txt").unwrap();
//! let mut reader = BufReader::new(file);
//! let mut input = String::new();
//! reader.read_to_string(&mut input).unwrap();
//! let contigs = parse(&input).unwrap();
//! for contig in contigs {
//!     println!("{} has {} arrays", contig.accession, contig.arrays.len());
//! }
//! ```

use nom::{
    branch::alt,
    bytes::complete::{tag, take_until},
    character::complete::{alpha1, char, digit1, line_ending, multispace1, not_line_ending},
    error::Error,
    multi::{many0, many1},
    sequence::{pair, tuple},
    Err, IResult,
};

#[derive(Debug, PartialEq)]
/// A single repeat and spacer.
pub struct RepeatSpacer<'a> {
    /// Sequence of the repeat.
    pub repeat: &'a str,
    /// Sequence of the spacer.
    pub spacer: &'a str,
    /// Zero-indexed inclusive start coordinate.
    pub start: usize,
    /// Zero-indexed exclusive end coordinate.
    pub end: usize,
    /// Zero-indexed inclusive start coordinate of the spacer.
    pub spacer_start: usize,
    /// Zero-indexed exclusive end coordinate of the spacer.
    pub spacer_end: usize,
    /// Zero-indexed inclusive start coordinate of the repeat.
    pub repeat_start: usize,
    /// Zero-indexed exclusive end coordinate of the repeat.
    pub repeat_end: usize,
}

#[derive(Debug, PartialEq)]
/// A single repeat, without a spacer. This is the last repeat in the CRISPR array.
pub struct RepeatOnly<'a> {
    /// Sequence of the repeat.
    pub repeat: &'a str,
    /// Zero-indexed inclusive start coordinate.
    pub start: usize,
    /// Zero-indexed exclusive end coordinate.
    pub end: usize,
}

/// Represents one component of a CRISPR array.
#[derive(Debug, PartialEq)]
pub enum Repeat<'a> {
    /// A repeat with a spacer
    WithSpacer(RepeatSpacer<'a>),
    /// A repeat without a spacer (the last repeat in the array)
    WithoutSpacer(RepeatOnly<'a>),
}

#[derive(Debug, PartialEq)]
/// A single CRISPR array.
pub struct Array<'a> {
    /// The nth CRISPR array in this genome/contig.
    pub order: usize,
    /// Zero-indexed inclusive start coordinate.
    pub start: usize,
    /// Zero-indexed exclusive end coordinate.
    pub end: usize,
    /// All of the repeat-spacer pairs in this CRISPR array.
    pub repeat_spacers: Vec<Repeat<'a>>,
}

#[derive(Debug)]
/// Represents all of the CRISPR arrays in a single contig or genome.
pub struct Contig<'a> {
    /// Accession of the contig/genome.
    pub accession: &'a str,
    /// Length of the contig/genome in base pairs.
    pub bp: usize,
    /// The CRISPR arrays in this contig/genome.
    pub arrays: Vec<Array<'a>>,
}

/// Parses the output of minCED for a single contig/genome.
pub fn parse(input: &str) -> Result<Vec<Contig>, Err<Error<&str>>> {
    let result = many0(parse_contig_arrays)(input);
    match result {
        Ok((_, contigs)) => Ok(contigs),
        Err(e) => Err(e),
    }
}

/// Parses the accession and arrays for a single contig/genome
fn parse_contig_arrays(input: &str) -> IResult<&str, Contig> {
    let result = tuple((
        parse_accession_line,
        skip_empty_line,
        many1(parse_array),
        parse_footer,
    ))(input);
    match result {
        Ok((remainder, ((accession, bp), _, arrays, _))) => Ok((
            remainder,
            Contig {
                accession,
                bp,
                arrays,
            },
        )),
        Err(e) => Err(e),
    }
}

/// Parses a single CRISPR array.
fn parse_array(input: &str) -> IResult<&str, Array> {
    let result = tuple((
        skip_empty_line,
        parse_crispr_order_and_coordinates,
        skip_empty_line,
        skip_one_line,
        skip_one_line,
        many1(parse_repeat_spacer_line),
        skip_one_line,
        skip_one_line,
    ))(input);
    match result {
        Ok((remainder, (_, (order, start, end), _, _, _, repeat_spacers, _, _))) => Ok((
            remainder,
            Array {
                order,
                start,
                end,
                repeat_spacers,
            },
        )),
        Err(e) => Err(e),
    }
}

/// Skips a line with text.
fn skip_one_line(input: &str) -> IResult<&str, ()> {
    let result = pair(not_line_ending, line_ending)(input);
    match result {
        Ok((remaining, _)) => Ok((remaining, ())),
        Err(e) => Err(e),
    }
}

/// Skips an empty line.
fn skip_empty_line(input: &str) -> IResult<&str, ()> {
    let result = line_ending(input);
    match result {
        Ok((remaining, _)) => Ok((remaining, ())),
        Err(e) => Err(e),
    }
}

/// Skips the four lines at the end of each contig.
fn parse_footer(input: &str) -> IResult<&str, ()> {
    let result = tuple((
        skip_empty_line,
        skip_one_line,
        skip_empty_line,
        skip_empty_line,
    ))(input);
    match result {
        Ok((remainder, _)) => Ok((remainder, ())),
        Err(e) => Err(e),
    }
}

/// Parses the order (i.e. the nth CRISPR array found for a given run of minCED) and start/end
/// coordinates of the array.
fn parse_crispr_order_and_coordinates(input: &str) -> IResult<&str, (usize, usize, usize)> {
    let result = tuple((
        tag("CRISPR"),
        char(' '),
        digit1,
        multispace1,
        tag("Range:"),
        char(' '),
        digit1,
        tag(" - "),
        digit1,
    ))(input);
    match result {
        Ok((remaining, (_, _, raw_order, _, _, _, start, _, end))) => Ok((
            remaining,
            (
                raw_order.parse::<usize>().unwrap() - 1,
                start.parse::<usize>().unwrap() - 1,
                end.parse::<usize>().unwrap(),
            ),
        )),
        Err(e) => Err(e),
    }
}

/// Parses the contig/genome accession and length
fn parse_accession_line(input: &str) -> IResult<&str, (&str, usize)> {
    let result = tuple((
        tag("Sequence '"),
        take_until("'"),
        tag("'"),
        char(' '),
        tag("("),
        take_until(" "),
        tag(" bp)"),
    ))(input);
    match result {
        Ok((remainder, (_, accession, _, _, _, bp, _))) => {
            Ok((remainder, (accession, bp.parse::<usize>().unwrap())))
        }
        Err(e) => Err(e),
    }
}

/// Parses a single repeat/spacer line
fn parse_repeat_spacer_line(input: &str) -> IResult<&str, Repeat> {
    alt((parse_repeat_with_spacer, parse_repeat_only))(input)
}

/// Parses a repeat entry that has no spacer. This is always the final repeat in the array.
fn parse_repeat_only(input: &str) -> IResult<&str, Repeat> {
    let result = tuple((digit1, multispace1, alpha1, multispace1))(input);
    match result {
        Ok((remaining, (raw_start, _, repeat, _))) => {
            let start = raw_start.parse::<usize>().unwrap() - 1;
            Ok((
                remaining,
                Repeat::WithoutSpacer(RepeatOnly {
                    repeat,
                    start,
                    end: start + repeat.len(),
                }),
            ))
        }
        Err(e) => Err(e),
    }
}

/// Parses a repeat and spacer entry.
fn parse_repeat_with_spacer(input: &str) -> IResult<&str, Repeat> {
    let result = tuple((
        digit1,
        multispace1,
        alpha1,
        multispace1,
        alpha1,
        not_line_ending,
        line_ending,
    ))(input);
    match result {
        Ok((remaining, (raw_start, _, repeat, _, spacer, _, _))) => {
            let start = raw_start.parse::<usize>().unwrap() - 1;
            Ok((
                remaining,
                Repeat::WithSpacer(RepeatSpacer {
                    repeat,
                    spacer,
                    start,
                    end: start + repeat.len() + spacer.len(),
                    repeat_start: start,
                    repeat_end: start + repeat.len(),
                    spacer_start: start + repeat.len(),
                    spacer_end: start + repeat.len() + spacer.len(),
                }),
            ))
        }
        Err(e) => Err(e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_array() {
        let input = "\nCRISPR 1   Range: 10648 - 10814
POSITION        REPEAT                          SPACER
--------        -----------------------------   ----------------------------------------
10648           CAAGTGCACCAACCAATCTCACCACCTCA   GGGGGTGCACTTAAAGGGGGTGCACTTGTCTCAAGTGCACCAAGAA  [ 29, 46 ]
10723           CAAGTGCACCAACCAATCTCACCACCTCA   CCATCTCACCACCTCTCAGGGGGTGCAGTTGTCT      [ 29, 34 ]
10786           CAAGTGCACCAACCAATCTCACCACCTCA
--------        -----------------------------   ----------------------------------------
Repeats: 3      Average Length: 29              Average Length: 40\n";
        let expected = Array {
            order: 0,
            start: 10647,
            end: 10814,
            repeat_spacers: vec![
                Repeat::WithSpacer(RepeatSpacer {
                    start: 10647,
                    end: 10722,
                    repeat_start: 10647,
                    repeat_end: 10676,
                    spacer_start: 10676,
                    spacer_end: 10722,
                    repeat: "CAAGTGCACCAACCAATCTCACCACCTCA",
                    spacer: "GGGGGTGCACTTAAAGGGGGTGCACTTGTCTCAAGTGCACCAAGAA",
                }),
                Repeat::WithSpacer(RepeatSpacer {
                    start: 10722,
                    end: 10785,
                    repeat_start: 10722,
                    repeat_end: 10751,
                    spacer_start: 10751,
                    spacer_end: 10785,
                    repeat: "CAAGTGCACCAACCAATCTCACCACCTCA",
                    spacer: "CCATCTCACCACCTCTCAGGGGGTGCAGTTGTCT",
                }),
                Repeat::WithoutSpacer(RepeatOnly {
                    start: 10785,
                    end: 10814,
                    repeat: "CAAGTGCACCAACCAATCTCACCACCTCA",
                }),
            ],
        };
        let (_, actual) = parse_array(input).unwrap();
        assert_eq!(expected, actual);
    }

    #[test]
    fn test_parse_crispr_order_and_coordinates() {
        let input = r"CRISPR 1   Range: 1214 - 1776";
        let expected = (0, 1213, 1776);
        let (_, actual) = parse_crispr_order_and_coordinates(input).unwrap();
        assert_eq!(expected, actual);
    }

    #[test]
    fn test_parse_accession_line() {
        let input = "Sequence 'MGYG000166779_38' (12280 bp)";
        let expected = ("MGYG000166779_38", 12280);
        let (_, actual) = parse_accession_line(input).unwrap();
        assert_eq!(expected, actual);
    }

    #[test]
    fn test_parse_repeat_spacer() {
        let input = "10723           CAAGTGCACCAACCAATCTCACCACCTCA   CCATCTCACCACCTCTCAGGGGGTGCAGTTGTCT      [ 29, 34 ]\n";
        let expected = RepeatSpacer {
            repeat: "CAAGTGCACCAACCAATCTCACCACCTCA",
            spacer: "CCATCTCACCACCTCTCAGGGGGTGCAGTTGTCT",
            start: 10722,
            end: 10785,
            repeat_start: 10722,
            repeat_end: 10751,
            spacer_start: 10751,
            spacer_end: 10785,
        };
        let (_, actual) = parse_repeat_spacer_line(input).unwrap();
        match actual {
            Repeat::WithSpacer(act) => {
                assert_eq!(expected, act);
            }
            _ => {
                unreachable!()
            }
        }
    }

    #[test]
    fn test_parse_repeat_only_line() {
        let input = "10786		CAAGTGCACCAACCAATCTCACCACCTCA\n";
        let expected = RepeatOnly {
            repeat: "CAAGTGCACCAACCAATCTCACCACCTCA",
            start: 10785,
            end: 10814,
        };
        let (_, actual) = parse_repeat_spacer_line(input).unwrap();
        match actual {
            Repeat::WithoutSpacer(act) => {
                assert_eq!(expected, act);
            }
            _ => {
                unreachable!()
            }
        }
    }

    #[test]
    fn test_parse_contig_arrays() {
        let input = "Sequence 'MGYG000242676_4' (164254 bp)

CRISPR 3   Range: 60487 - 61025
POSITION	REPEAT				SPACER
--------	------------------------------------	--------------------------
60487		TTTAATAACCCTATATAATTTCTACTATTGTAGATA	TCTCCTTTGTAACTTCTTTGATTCGG	[ 36, 26 ]
60549		TTTAATAACCCTATATAATTTCTACTGTCGTAGATA	TTGTTCTTTTATATGTGTACATAGCTAGA	[ 36, 29 ]
60990		TTTAATAACCCTATATAATTTCTACTTTTTTGATTA
--------	------------------------------------	--------------------------
Repeats: 9	Average Length: 36		Average Length: 26

CRISPR 4   Range: 157550 - 157915
POSITION	REPEAT				SPACER
--------	------------------------------------	------------------------------
157550		GTTTTACTACCTTATAGATTTACACTATTCTCAAAC	GAGGGGTTGTCCTTCATGTACTCTTTACCT	[ 36, 30 ]
157748		GTTTTACTACCTTATAGATTTACACTATTCTCAAAC	GGGCTTATACTCTGACTTTCAACAAGTTAG	[ 36, 30 ]
157814		GTTTTACTACCTTATAGATTTACACTATTCTCAAAC	CCGATTTTTTCATTGCCAAAACGATATTTT	[ 36, 30 ]
157880		GTTTTACTACCTTATAGATTTACACTATTCTCAAAC
--------	------------------------------------	------------------------------
Repeats: 6	Average Length: 36		Average Length: 30

Time to find repeats: 22 ms


";
        let (_, contig) = parse_contig_arrays(input).unwrap();
        assert_eq!(contig.accession, "MGYG000242676_4");
        assert_eq!(contig.bp, 164254);
        assert_eq!(contig.arrays.len(), 2);
    }

    #[test]
    fn test_parse() {
        let input = "Sequence 'MGYG000166779_38' (12280 bp)

CRISPR 1   Range: 10648 - 10814
POSITION	REPEAT				SPACER
--------	-----------------------------	----------------------------------------
10648		CAAGTGCACCAACCAATCTCACCACCTCA	GGGGGTGCACTTAAAGGGGGTGCACTTGTCTCAAGTGCACCAAGAA	[ 29, 46 ]
10723		CAAGTGCACCAACCAATCTCACCACCTCA	CCATCTCACCACCTCTCAGGGGGTGCAGTTGTCT	[ 29, 34 ]
10786		CAAGTGCACCAACCAATCTCACCACCTCA	
--------	-----------------------------	----------------------------------------
Repeats: 3	Average Length: 29		Average Length: 40

Time to find repeats: 3 ms


Sequence 'MGYG000166779_43' (11302 bp)

CRISPR 2   Range: 4 - 1413
POSITION	REPEAT				SPACER
--------	------------------------------------	-----------------------------
4		GTTGTGGTTTGATGTAGGAATCAAAAGATATACAAC	ACGGGTGCACTTTCGATGTCGCACTTTTTG	[ 36, 30 ]
70		GTTGTGGTTTGATGTAGGAATCAAAAGATATACAAC	TATACATCATCGTACATATAAGCATACAG	[ 36, 29 ]
135		GTTGTGGTTTGATGTAGGAATCAAAAGATATACAAC	GAAAAATCAGAGCCCAAAGTACGAGTAAC	[ 36, 29 ]
200		GTTGTGGTTTGATGTAGGAATCAAAAGATATACAAC	CCAGTTCCCGAATTTGATGCTCTTGGCAT	[ 36, 29 ]
265		GTTGTGGTTTGATGTAGGAATCAAAAGATATACAAC	ACTTACAACAACAACAATAACAATAAATG	[ 36, 29 ]
330		GTTGTGGTTTGATGTAGGAATCAAAAGATATACAAC	ATACGTGTGCTCTATATACGCACCCATTGG	[ 36, 30 ]
396		GTTGTGGTTTGATGTAGGAATCAAAAGATATACAAC	GGAGCTCTTTCGATGTCGCACTTTCTGAAG	[ 36, 30 ]
462		GTTGTGGTTTGATGTAGGAATCAAAAGATATACAAC	CGTGCTCGCTTTGAATTTGTAGAACCCGA	[ 36, 29 ]
527		GTTGTGGTTTGATGTAGGAATCAAAAGATATACAAC	TCTCGACACTATTTCTAACGAGGAAATTAA	[ 36, 30 ]
593		GTTGTGGTTTGATGTAGGAATCAAAAGATATACAAC	GCGCTGAGAAGTTACCACCGACCGCTTGA	[ 36, 29 ]
658		GTTGTGGTTTGATGTAGGAATCAAAAGATATACAAC	AATACTAAACCAAGATTGCCAAAGGTCCA	[ 36, 29 ]
723		GTTGTGGTTTGATGTAGGAATCAAAAGATATACAAC	AGATGATCTACGCTCAATATTAGAAAAAC	[ 36, 29 ]
788		GTTGTGGTTTGATGTAGGAATCAAAAGATATACAAC	GTATCTGCGGAACAAGTACAGAGAACATGA	[ 36, 30 ]
854		GTTGTGGTTTGATGTAGGAATCAAAAGATATACAAC	CGAACCTAATACGGCTTTAGCCTTTTTGCA	[ 36, 30 ]
920		GTTGTGGTTTGATGTAGGAATCAAAAGATATACAAC	AATGAGTACCAAAAGCAAAGAACAAATCGA	[ 36, 30 ]
986		GTTGTGGTTTGATGTAGGAATCAAAAGATATACAAC	TATATTTTTGTGCGTTACCCGTCCGTGAGG	[ 36, 30 ]
1052		GTTGTGGTTTGATGTAGGAATCAAAAGATATACAAC	TTACTTACGACTATTACGACCAGGTGAAC	[ 36, 29 ]
1117		GTTGTGGTTTGATGTAGGAATCAAAAGATATACAAC	ATAATTATAATCGGAAATCAAGCGGATAA	[ 36, 29 ]
1182		GTTGTGGTTTGATGTAGGAATCAAAAGATATACAAC	TTTGAATTAGATTCGGCAACCTTAGCATT	[ 36, 29 ]
1247		GTTGTGGTTTGATGTAGGAATCAAAAGATATACAAC	TTTTCATCATATTCATAAGAATAGCGACC	[ 36, 29 ]
1312		GTTGTGGTTTGATGTAGGAATCAAAAGATATACAAC	CCATACGCTCCTTGGTGGTCTTGGTAAGGA	[ 36, 30 ]
1378		GTTGTGGTTTGATGTAGAAATCAAAAGACATACAAC	
--------	------------------------------------	-----------------------------
Repeats: 22	Average Length: 36		Average Length: 29

Time to find repeats: 3 ms


Sequence 'MGYG000242676_4' (164254 bp)

CRISPR 3   Range: 60487 - 61025
POSITION	REPEAT				SPACER
--------	------------------------------------	--------------------------
60487		TTTAATAACCCTATATAATTTCTACTATTGTAGATA	TCTCCTTTGTAACTTCTTTGATTCGG	[ 36, 26 ]
60549		TTTAATAACCCTATATAATTTCTACTGTCGTAGATA	TTGTTCTTTTATATGTGTACATAGCTAGA	[ 36, 29 ]
60614		TTTAATAACCCTATATAATTTCTACTATTGTAGATA	ACCTCCTTTGGATTTTCAGCAAATCAGG	[ 36, 28 ]
60678		TTTAATAACCCTATATAATTTCTACTATTTTAGATA	ATACTGCTTGTTCTGTAAAAATTTTG	[ 36, 26 ]
60740		TTTAATAACTCTATATAATTTCTACTATTGTAGATG	GAGTTCTCCAACCGTTTGCGGCAATA	[ 36, 26 ]
60802		TTTAATAACCCTATATAATTTCTACTATTGTAGATA	ACGGTTGAATCAATGAGAAATGTTGTG	[ 36, 27 ]
60865		TTTAATAACCCTATATAATTTCTACTATTGTAGATA	TGATATTGACGGTGACCTGATTAACCG	[ 36, 27 ]
60928		TTTAATAACCCTATATAATTTCTACTATTGTAGATA	TGTCAATCACATCTGTGACCGCAAGG	[ 36, 26 ]
60990		TTTAATAACCCTATATAATTTCTACTTTTTTGATTA	
--------	------------------------------------	--------------------------
Repeats: 9	Average Length: 36		Average Length: 26

CRISPR 4   Range: 157550 - 157915
POSITION	REPEAT				SPACER
--------	------------------------------------	------------------------------
157550		GTTTTACTACCTTATAGATTTACACTATTCTCAAAC	GAGGGGTTGTCCTTCATGTACTCTTTACCT	[ 36, 30 ]
157616		GTTTTACTACCTTATAGATTTACACTATTCTCAAAC	ATACAAATGCATTGCCGAGGACAGTGTTTT	[ 36, 30 ]
157682		GTTTTACTACCTTATAGATTTACACTATTCTCAAAC	TATACGGTTTGCCCGTGCAGTCTTGTACAA	[ 36, 30 ]
157748		GTTTTACTACCTTATAGATTTACACTATTCTCAAAC	GGGCTTATACTCTGACTTTCAACAAGTTAG	[ 36, 30 ]
157814		GTTTTACTACCTTATAGATTTACACTATTCTCAAAC	CCGATTTTTTCATTGCCAAAACGATATTTT	[ 36, 30 ]
157880		GTTTTACTACCTTATAGATTTACACTATTCTCAAAC	
--------	------------------------------------	------------------------------
Repeats: 6	Average Length: 36		Average Length: 30

Time to find repeats: 22 ms


Sequence 'MGYG000273829_14' (62198 bp)

CRISPR 5   Range: 15191 - 17205
POSITION	REPEAT				SPACER
--------	------------------------------------	-----------------------------
15191		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	ATCGCTGAACCTACAACAGACGCAAGAACA	[ 36, 30 ]
15257		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	GATATTGTCATACCTAAGTAAATAGGTGCG	[ 36, 30 ]
15323		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	TACAATCAGTCTATAACATTTGCAACTACG	[ 36, 30 ]
15389		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	TATTATAGACAGCAAGCAACTTGATGTAT	[ 36, 29 ]
15454		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	CTTGAATTTGGGGAGATGTTCTCAGCTGGT	[ 36, 30 ]
15520		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	AAAGTTTGCTGACAGGGACATTCAAAGCCG	[ 36, 30 ]
15586		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	AAACCTGTCTGTCCGATCTGCACCATATAT	[ 36, 30 ]
15652		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	ACCGCTATTGCGCTGCAGCATCCACAAGGA	[ 36, 30 ]
15718		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	GCATCTTCCTGCGCTCTCTCTGAAAACATG	[ 36, 30 ]
15784		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	CGAAGCCTAAAGCTCATTTCGCTTAGGCTT	[ 36, 30 ]
15850		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	AAACATGGTGTTGATGTCAAAGAGCTGTAT	[ 36, 30 ]
15916		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	TGGGCGAATATAAATTCCATCGGTGGCAAG	[ 36, 30 ]
15982		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	AATATTGGATGATGTGTATGGCATTTTACT	[ 36, 30 ]
16048		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	AGTGTATATGTGAACCCTGCTCCCAGTGCT	[ 36, 30 ]
16114		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	AAAGACCGGAGCAAAGATGTCCGGGAGCCG	[ 36, 30 ]
16180		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	TGAAAGTGGTGTAATTGTTATAACTCATTG	[ 36, 30 ]
16246		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	AGACAACAGGTGTGGAAGCATATGTCTTTA	[ 36, 30 ]
16312		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	TGCTGCATAGGTGTGTATTTTCTCATGTCG	[ 36, 30 ]
16378		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	GGTAATGGTGGTGGCGGTTATACCGCAACT	[ 36, 30 ]
16444		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	ATGGTCGGGGCTACATATTACGCCGCAGTA	[ 36, 30 ]
16510		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	CGTGAGGTCTCCGACCGTGAAAACAGTTCT	[ 36, 30 ]
16576		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	ACGAACTTAGTACCCTTTTCTGGGCGGCAT	[ 36, 30 ]
16642		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	CCGCAGGTGCTACCGCTGTTATACTCTGTT	[ 36, 30 ]
16708		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	CGTAAATCGTTGGCGAAACGCTACCAACTG	[ 36, 30 ]
16774		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	CCTCGGTCTGCTCTAACAGATCCCCCAAGT	[ 36, 30 ]
16840		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	ACAGAGAAAGAAAGAGAGATTAACGACTAC	[ 36, 30 ]
16906		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	TGAAACGGAGTGGACAGGTAAAGGAATGGG	[ 36, 30 ]
16972		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	TGCGGTCCCTTGGTTCCGTCAACAACATCA	[ 36, 30 ]
17038		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	TGTCCTATTCCCTTTTATGCTGCGTGTATA	[ 36, 30 ]
17104		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	AATACAAGCATAAAGAACGAACCGCAACGG	[ 36, 30 ]
17170		GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT	
--------	------------------------------------	-----------------------------
Repeats: 31	Average Length: 36		Average Length: 29

Time to find repeats: 9 ms


";
        let contigs = parse(input).unwrap();
        assert_eq!(contigs.len(), 4);
        let array_count: usize = contigs.iter().map(|c| c.arrays.len()).sum();
        assert_eq!(array_count, 5);
    }
}
