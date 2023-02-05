use nom::{
    branch::alt,
    bytes::complete::{tag, take_until},
    character::complete::{
        alpha1, char, digit1, line_ending, multispace0, multispace1, not_line_ending,
    },
    multi::{many0, many1},
    sequence::{pair, separated_pair, tuple},
    IResult,
};

// #[derive(Debug)]
// pub struct MincedError {}

#[derive(Debug, PartialEq)]
pub struct RepeatSpacer<'a> {
    /// sequence of the repeat
    repeat: &'a str,
    /// sequence of the spacer
    spacer: Option<&'a str>,
    /// zero-indexed inclusive start coordinate
    start: usize,
    /// zero-indexed exclusive end coordinate
    end: usize,
}

#[derive(Debug, PartialEq)]
pub struct Array<'a> {
    /// The nth CRISPR array in this genome/contig
    pub order: usize,
    /// zero-indexed inclusive start coordinate
    pub start: usize,
    /// zero-indexed exclusive end coordinate
    pub end: usize,
    /// all of the repeat-spacer pairs in this CRISPR array
    pub repeat_spacers: Vec<RepeatSpacer<'a>>,
}

// pub fn parse(input: &str) -> Result<Vec<Array>, MincedError> {
//     let result = many0(parse_array)(input);
//     match result {
//         Ok((_, arrays)) => Ok(arrays),
//         Err(_) => Err(MincedError {}),
//     }
// }

// fn parse_contig_arrays(input: &str) -> IResult<&str, Vec<Array>> {
//     let result = tuple((parse_accession_line, many1(parse_array), parse_footer))(input);
// }

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

fn skip_one_line(input: &str) -> IResult<&str, ()> {
    let result = pair(not_line_ending, line_ending)(input);
    match result {
        Ok((remaining, _)) => Ok((remaining, ())),
        Err(e) => Err(e),
    }
}

fn skip_empty_line(input: &str) -> IResult<&str, ()> {
    let result = line_ending(input);
    match result {
        Ok((remaining, _)) => Ok((remaining, ())),
        Err(e) => Err(e),
    }
}

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

// fn parse_crispr_range(input: &str) -> IResult<&str, (usize, usize)> {
//     let result = tuple((tag("Range:"), char(' '), digit1, tag(" - "), digit1))(input);
//     match result {
//         Ok((remaining, (_, _, start, _, end))) => Ok((
//             remaining,
//             (
//                 start.parse::<usize>().unwrap() - 1,
//                 end.parse::<usize>().unwrap(),
//             ),
//         )),
//         Err(e) => Err(e),
//     }
// }

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

// fn parse_array_header(input: &str) -> IResult<&str, (usize, (usize, usize))> {
//     separated_pair(parse_crispr_order, multispace1, parse_crispr_range)(input)
// }

// // fn parse_header(input: &str) -> IResult<&str, Header> {
// //     let result = separated_pair(parse_accession_line, multispace1)(input);
// //     match result {
// //         Ok((remaining, ((accession, bp), (order, (start, end))))) => Ok((
// //             remaining,
// //             Header {
// //                 accession,
// //                 bp,
// //                 order,
// //                 start,
// //                 end,
// //             },
// //         )),
// //         Err(e) => Err(e),
// //     }
// // }

fn parse_repeat_spacer_line(input: &str) -> IResult<&str, RepeatSpacer> {
    alt((parse_repeat_with_spacer, parse_repeat_only))(input)
}

fn parse_repeat_only(input: &str) -> IResult<&str, RepeatSpacer> {
    let result = tuple((digit1, multispace1, alpha1, multispace1))(input);
    match result {
        Ok((remaining, (raw_start, _, repeat, _))) => {
            let start = raw_start.parse::<usize>().unwrap() - 1;
            Ok((
                remaining,
                RepeatSpacer {
                    repeat,
                    spacer: None,
                    start,
                    end: start + repeat.len(),
                },
            ))
        }
        Err(e) => Err(e),
    }
}

fn parse_repeat_with_spacer(input: &str) -> IResult<&str, RepeatSpacer> {
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
            dbg!("{raw_start} {repeat} {spacer}");
            let start = raw_start.parse::<usize>().unwrap() - 1;
            Ok((
                remaining,
                RepeatSpacer {
                    repeat,
                    spacer: Some(spacer),
                    start,
                    end: start + repeat.len() + spacer.len(),
                },
            ))
        }
        Err(e) => Err(e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Here is a real file and the python coordinates for those sequences.
    // For the range, we should return: (start - 1, end)
    // For each repeat-spacer, we should return (start - 1, start - 1 + length(repeat+spacer)
    // Sequence 'MGYG000166779_38' (12280 bp)

    // CRISPR 1   Range: 10648 - 10814
    // POSITION	REPEAT				SPACER
    // --------	-----------------------------	----------------------------------------
    // 10648		CAAGTGCACCAACCAATCTCACCACCTCA	GGGGGTGCACTTAAAGGGGGTGCACTTGTCTCAAGTGCACCAAGAA	[ 29, 46 ]
    // 10723		CAAGTGCACCAACCAATCTCACCACCTCA   CCATCTCACCACCTCTCAGGGGGTGCAGTTGTCT	[ 29, 34 ]
    // 10786		CAAGTGCACCAACCAATCTCACCACCTCA
    // --------	-----------------------------	----------------------------------------
    // Repeats: 3	Average Length: 29		Average Length: 40
    //     In [18]: str(fna[10722:10785].seq)
    // Out[18]: 'CAAGTGCACCAACCAATCTCACCACCTCACCATCTCACCACCTCTCAGGGGGTGCAGTTGTCT'

    // In [19]: str(fna[10647:10814].seq)
    // Out[19]: 'CAAGTGCACCAACCAATCTCACCACCTCAGGGGGTGCACTTAAAGGGGGTGCACTTGTCTCAAGTGCACCAAGAACAAGTGCACCAACCAATCTCACCACCTCACCATCTCACCACCTCTCAGGGGGTGCAGTTGTCTCAAGTGCACCAACCAATCTCACCACCTCA'

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
                RepeatSpacer {
                    start: 10647,
                    end: 10722,
                    repeat: "CAAGTGCACCAACCAATCTCACCACCTCA",
                    spacer: Some("GGGGGTGCACTTAAAGGGGGTGCACTTGTCTCAAGTGCACCAAGAA"),
                },
                RepeatSpacer {
                    start: 10722,
                    end: 10785,
                    repeat: "CAAGTGCACCAACCAATCTCACCACCTCA",
                    spacer: Some("CCATCTCACCACCTCTCAGGGGGTGCAGTTGTCT"),
                },
                RepeatSpacer {
                    start: 10785,
                    end: 10814,
                    repeat: "CAAGTGCACCAACCAATCTCACCACCTCA",
                    spacer: None,
                },
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
            spacer: Some("CCATCTCACCACCTCTCAGGGGGTGCAGTTGTCT"),
            start: 10722,
            end: 10785,
        };
        let (_, actual) = parse_repeat_spacer_line(input).unwrap();
        assert_eq!(expected, actual);
    }

    #[test]
    fn test_parse_repeat_only_line() {
        let input = "10786		CAAGTGCACCAACCAATCTCACCACCTCA\n";
        let expected = RepeatSpacer {
            repeat: "CAAGTGCACCAACCAATCTCACCACCTCA",
            spacer: None,
            start: 10785,
            end: 10814,
        };
        let (_, actual) = parse_repeat_spacer_line(input).unwrap();
        assert_eq!(expected, actual);
    }

    //     #[test]
    //     fn test_parse_contig_arrays() {
    //         let input = r"Sequence 'MGYG000242676_4' (164254 bp)

    // CRISPR 3   Range: 60487 - 61025
    // POSITION	REPEAT				SPACER
    // --------	------------------------------------	--------------------------
    // 60487		TTTAATAACCCTATATAATTTCTACTATTGTAGATA	TCTCCTTTGTAACTTCTTTGATTCGG	[ 36, 26 ]
    // 60549		TTTAATAACCCTATATAATTTCTACTGTCGTAGATA	TTGTTCTTTTATATGTGTACATAGCTAGA	[ 36, 29 ]
    // 60990		TTTAATAACCCTATATAATTTCTACTTTTTTGATTA
    // --------	------------------------------------	--------------------------
    // Repeats: 9	Average Length: 36		Average Length: 26

    // CRISPR 4   Range: 157550 - 157915
    // POSITION	REPEAT				SPACER
    // --------	------------------------------------	------------------------------
    // 157550		GTTTTACTACCTTATAGATTTACACTATTCTCAAAC	GAGGGGTTGTCCTTCATGTACTCTTTACCT	[ 36, 30 ]
    // 157748		GTTTTACTACCTTATAGATTTACACTATTCTCAAAC	GGGCTTATACTCTGACTTTCAACAAGTTAG	[ 36, 30 ]
    // 157814		GTTTTACTACCTTATAGATTTACACTATTCTCAAAC	CCGATTTTTTCATTGCCAAAACGATATTTT	[ 36, 30 ]
    // 157880		GTTTTACTACCTTATAGATTTACACTATTCTCAAAC
    // --------	------------------------------------	------------------------------
    // Repeats: 6	Average Length: 36		Average Length: 30

    // Time to find repeats: 22 ms";
    //     }

    //     #[test]
    //     fn test_parse_array() {
    //         let input = r"Sequence 'GUT_GENOME226819_47' (4332 bp)

    // CRISPR 1   Range: 1214 - 1776
    // POSITION        REPEAT                          SPACER
    // --------        -----------------------------   -----------------------------------------------
    // 1214            TGGTTTAGTTTTAGATGATTCTACTTTTA   CTAATAACACTGCTAAGATTGGTGGTGCAATTTATAACTCTGCTGA  [ 29, 46 ]
    // 1289            TTGCTTTGTTGTAGGTAATTCTACTTTTG   CTAATAACACTGCTACAAGTAATGGTGGAGTAATCTTTAATTATGGTAT       [ 29, 49 ]
    // 1367            TGGTTTTGTTGTAGGTAATTCTACTTTTG   TGAATAATAGTGCTGCAGACGGTGCTGGTGCAATCCTGAATGGTGGCCG       [ 29, 49 ]
    // 1445            TGGTTTTGTTGTAGGTAATTCTACTTTTG   CTAATAACACTGCTACAAGTAAAGGTGGTGCCATTTATAATTATGGTAT       [ 29, 49 ]
    // 1523            TGGTTTTGTTGTAGGTAATTCTACTTTTG   CTAATAACACTGCAGAAGATGCCGGTGCAGTTTATAATGAGGGTGA  [ 29, 46 ]
    // 1598            TAACTCTGTTGTAGGTAATTCTACTTTTG   TTAATAATACTGCAACTTCAATAGGTGGAGCTATAATTAATAATGG  [ 29, 46 ]
    // 1673            TAAATTAGTAGTTGATAATTCTGCATTTG   AAGATAATGCTGCTAATTATTATGGTGGAGCTATCTTTAACTGGGA  [ 29, 46 ]
    // 1748            TGATTTACAAGTTACTAACTCTGCATTTG
    // --------        -----------------------------   -----------------------------------------------
    // Repeats: 8      Average Length: 29              Average Length: 47

    // Time to find repeats: 10 ms";
    //         let expected = Array {
    //             order: 1,
    //             start: 1213,
    //             end: 1776,
    //             accession: "GUT_GENOME226819_47",
    //             bp: 4332,
    //             repeat_spacers: vec![
    //                 RepeatSpacer {
    //                     repeat: "TGGTTTAGTTTTAGATGATTCTACTTTTA",
    //                     spacer: Some("CTAATAACACTGCTAAGATTGGTGGTGCAATTTATAACTCTGCTGA"),
    //                     start: 1213,
    //                     end: 1213 + 29 + 46 + 1,
    //                 },
    //                 RepeatSpacer {
    //                     repeat: "TTGCTTTGTTGTAGGTAATTCTACTTTTG",
    //                     spacer: Some("CTAATAACACTGCTACAAGTAATGGTGGAGTAATCTTTAATTATGGTAT"),
    //                     start: 1288,
    //                     end: 1288 + 29 + 49 + 1,
    //                 },
    //                 RepeatSpacer {
    //                     repeat: "TGGTTTTGTTGTAGGTAATTCTACTTTTG",
    //                     spacer: Some("TGAATAATAGTGCTGCAGACGGTGCTGGTGCAATCCTGAATGGTGGCCG"),
    //                     start: 1366,
    //                     end: 1366 + 29 + 49 + 1,
    //                 },
    //                 RepeatSpacer {
    //                     repeat: "TGGTTTTGTTGTAGGTAATTCTACTTTTG",
    //                     spacer: Some("CTAATAACACTGCTACAAGTAAAGGTGGTGCCATTTATAATTATGGTAT"),
    //                     start: 1444,
    //                     end: 1444 + 29 + 49 + 1,
    //                 },
    //                 RepeatSpacer {
    //                     repeat: "TGGTTTTGTTGTAGGTAATTCTACTTTTG",
    //                     spacer: Some("CTAATAACACTGCAGAAGATGCCGGTGCAGTTTATAATGAGGGTGA"),
    //                     start: 1522,
    //                     end: 1522 + 29 + 46 + 1,
    //                 },
    //                 RepeatSpacer {
    //                     repeat: "TAACTCTGTTGTAGGTAATTCTACTTTTG",
    //                     spacer: Some("TTAATAATACTGCAACTTCAATAGGTGGAGCTATAATTAATAATGG"),
    //                     start: 1597,
    //                     end: 1597 + 29 + 46 + 1,
    //                 },
    //                 RepeatSpacer {
    //                     repeat: "TAAATTAGTAGTTGATAATTCTGCATTTG",
    //                     spacer: Some("AAGATAATGCTGCTAATTATTATGGTGGAGCTATCTTTAACTGGGA"),
    //                     start: 1672,
    //                     end: 1672 + 29 + 46 + 1,
    //                 },
    //                 RepeatSpacer {
    //                     repeat: "TGATTTACAAGTTACTAACTCTGCATTTG",
    //                     spacer: None,
    //                     start: 1747,
    //                     end: 1747 + 29 + 1,
    //                 },
    //             ],
    //         };
    //         let (_, actual) = parse_array(input).unwrap();
    //         assert_eq!(expected, actual);
    //     }
}
