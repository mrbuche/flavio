// Start by doing what's necessary,
// then do what's possible,
// and suddenly you are doing the impossible.
// --- St. Francis of Assisi

// I can accept failure,
// everyone fails at something.
// But I can't accept not trying.
// --- Michael Jordan

fn main() {
    println!(
        "\x1b[1m{} {}\x1b[0m",
        env!("CARGO_PKG_NAME"),
        env!("CARGO_PKG_VERSION")
    );
    println!("\x1b[2m{}\x1b[0m\n", env!("CARGO_PKG_AUTHORS"));
    println!("\x1b[1;96mf l a v i o\x1b[0m");
    println!("\x1b[2;96mflavio welcomes you\x1b[0m\n");
}
