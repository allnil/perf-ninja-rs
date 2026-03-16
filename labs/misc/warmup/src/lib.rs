#[cfg(test)]
mod tests;

pub fn solution(arr: &[i32], n: usize) -> i32 {
    let N = n as i32;
    (N * (N + 1)) / 2
}
