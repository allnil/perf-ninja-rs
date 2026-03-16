use std::cmp::max;

#[cfg(test)]
mod tests;

pub const SEQUENCE_SIZE: usize = 200; // The length of the generated sequences.
pub const SEQUENCE_COUNT: usize = 16; // The number of sequences to generate for both sequence collections.

pub type Sequence = [u8; SEQUENCE_SIZE];
pub type AlignResult = [i16; SEQUENCE_COUNT];

// The alignment algorithm which computes the alignment of the given sequence
// pairs.
pub fn compute_alignment(sequences1: &[Sequence], sequences2: &[Sequence]) -> AlignResult {
    let mut result: AlignResult = [0; SEQUENCE_COUNT];

    type Score = i16;
    type Column = [Score; SEQUENCE_SIZE + 1];

    let gap_open: Score = -11;
    let gap_extension: Score = -1;
    let match_score: Score = 6; // previously 'match'
    let mismatch: Score = -4;

    let mut score_column = [[0i16; SEQUENCE_COUNT]; SEQUENCE_SIZE + 1];
    let mut horizontal_gap_column = [[0i16; SEQUENCE_COUNT]; SEQUENCE_SIZE + 1];
    let mut last_vertical_gap = [0i16; SEQUENCE_COUNT];
    let mut last_diagonal_score = [0i16; SEQUENCE_COUNT];

    for seq in 0..SEQUENCE_COUNT {
        horizontal_gap_column[0][seq] = gap_open;
        last_vertical_gap[seq] = gap_open;

        for i in 1..=SEQUENCE_SIZE {
            score_column[i][seq] = last_vertical_gap[seq];
            horizontal_gap_column[i][seq] = last_vertical_gap[seq] + gap_open;
            last_vertical_gap[seq] += gap_extension;
        }
    }

    for col in 1..=SEQUENCE_SIZE {
        for seq in 0..SEQUENCE_COUNT {
            last_diagonal_score[seq] = score_column[0][seq];
            score_column[0][seq] = horizontal_gap_column[0][seq];
            last_vertical_gap[seq] = horizontal_gap_column[0][seq] + gap_open;
            horizontal_gap_column[0][seq] += gap_extension;
        }

        for row in 1..=SEQUENCE_SIZE {
            for seq in 0..SEQUENCE_COUNT {
                let diag = sequences2[seq][col - 1];

                // Compute next score from diagonal direction with match/mismatch.
                let eq = (sequences1[seq][row - 1] == diag) as i16;
                let diag_score = mismatch + eq * (match_score - mismatch);
                let mut best_cell_score = last_diagonal_score[seq] + diag_score;

                // Determine best score from diagonal, vertical, or horizontal
                // direction.
                best_cell_score = max(best_cell_score, last_vertical_gap[seq]);
                best_cell_score = max(best_cell_score, horizontal_gap_column[row][seq]);

                // Cache next diagonal value and store optimum in score_column.
                last_diagonal_score[seq] = score_column[row][seq];
                score_column[row][seq] = best_cell_score;
                // Compute the next values for vertical and horizontal gap.
                best_cell_score += gap_open;
                last_vertical_gap[seq] += gap_extension;
                horizontal_gap_column[row][seq] += gap_extension;

                // Store optimum between gap open and gap extension.
                last_vertical_gap[seq] = max(last_vertical_gap[seq], best_cell_score);
                horizontal_gap_column[row][seq] =
                    max(horizontal_gap_column[row][seq], best_cell_score);
            }
        }
    }

    // Report the best score.

    for seq in 0..SEQUENCE_COUNT {
        result[seq] = score_column[SEQUENCE_SIZE][seq];
    }

    result
}

// C++ version includes a branchless max which is ported here (copied from crate lokacore) but
// not used as it doesn't seem to affect branch-misses, and is slower.
//
// Clang-12 compiler generates branches for std::max, which are often mispredicted
// in this benchmark. That's the reason we provide branchless version of max function.
// fn max(a: i16, b: i16) -> i16 {
//     a ^ ((a ^ b) & ((a < b) as i16).wrapping_neg())
// }

// Initialises a pair of sequence collections given a fixed sequence size.
pub fn init() -> (Vec<Sequence>, Vec<Sequence>) {
    use rand::distributions::Uniform;
    use rand::prelude::*;
    let mut random_engine = thread_rng();

    // Simulate DNA alphabet with 4 symbols.
    let symbol_distribution = Uniform::from(0..4);

    let mut generate_sequences = || -> Vec<Sequence> {
        let mut sequences: Vec<Sequence> = vec![[0u8; SEQUENCE_SIZE]; SEQUENCE_COUNT];
        for sequence in &mut sequences {
            sequence.fill_with(|| symbol_distribution.sample(&mut random_engine));
        }
        sequences
    };

    (generate_sequences(), generate_sequences())
}
