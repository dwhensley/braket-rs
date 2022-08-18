mod complex;
mod vector;

use complex::C64;
use vector::{Bra, Ket, Vector};

fn main() {
    let ket1: Vector<Ket, 3> =
        Vector::new([C64::new(1.0, 0.0), C64::new(0.5, 0.5), C64::new(0.75, 0.25)]);

    let ket2: Vector<Ket, 5> = Vector::new([
        C64::new(1.0, 0.0),
        C64::new(0.5, 0.5),
        C64::new(0.75, 0.25),
        C64::new(-0.55, 0.8),
        C64::new(0.1, 0.2),
    ]);

    let bra1 = ket1.to_bra();
    let bra2 = ket2.to_bra();

    let bra3: Vector<Bra, 3> =
        Vector::new([C64::new(0.8, 0.2), C64::new(-1.1, 0.9), C64::new(0.6, 0.1)]);

    println!("\nKet1: {}", &ket1);
    println!("\nKet2: {}", &ket2);
    println!("\nBra1: {}", &bra1);
    println!("\nBra2: {}", &bra2);
    println!("\nKet1 Normalized: {}", &ket1.to_normalized());
    println!(
        "\nInner Product of Ket1 and Bra3: {:?}",
        ket1.inner_prod(&bra3)
    );
}
