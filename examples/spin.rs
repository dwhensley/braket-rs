use braket::complex::C64;
use braket::operator::HermitianMatrix;
use braket::vector::{Bra, InnerProductDualSpace, Ket, Vector};

fn main() {
    let one_over_sqrt2 = 1. / 2.0_f64.sqrt();

    let u: Vector<Ket, 2> = Vector::from_arr([C64::new(1.0, 0.0), C64::new(0.0, 0.0)]);
    let d: Vector<Ket, 2> = Vector::from_arr([C64::new(0.0, 0.0), C64::new(1.0, 0.0)]);

    let l: Vector<Ket, 2> = one_over_sqrt2 * u - one_over_sqrt2 * d;
    let r: Vector<Ket, 2> = one_over_sqrt2 * u + one_over_sqrt2 * d;

    let i: Vector<Ket, 2> = one_over_sqrt2 * u + C64::i() * one_over_sqrt2 * d;
    let o: Vector<Ket, 2> = one_over_sqrt2 * u - C64::i() * one_over_sqrt2 * d;

    println!("\nUp spin state vector: {}", &u);
    println!("\nDown spin state vector: {}", &d);
    println!("\nLeft spin state vector: {}", &l);
    println!("\nRight spin state vector: {}", &r);
    println!("\nIn spin state vector: {}", &i);
    println!("\nOut spin state vector: {}", &o);

    println!("\n<o|u><u|o>: {}", (o.to_bra() * u) * (u.to_bra() * o));

    let sigma_z: HermitianMatrix<2> = HermitianMatrix::from_arr([
        [C64::new(1.0, 0.0), C64::new(0.0, 0.0)],
        [C64::new(0.0, 0.0), C64::new(-1.0, 0.0)],
    ])
    .unwrap();
    let sigma_x: HermitianMatrix<2> = HermitianMatrix::from_arr([
        [C64::new(0.0, 0.0), C64::new(1.0, 0.0)],
        [C64::new(1.0, 0.0), C64::new(0.0, 0.0)],
    ])
    .unwrap();
    let sigma_y: HermitianMatrix<2> = HermitianMatrix::from_arr([
        [C64::new(0.0, 0.0), C64::new(0.0, -1.0)],
        [C64::new(0.0, 1.0), C64::new(0.0, 0.0)],
    ])
    .unwrap();

    println!("\nSpin operator (z): {}", &sigma_z);
    println!("\nSpin operator (x): {}", &sigma_x);
    println!("\nSpin operator (y): {}", &sigma_y);

    println!("\nSz|u>: {}", sigma_z * u);
    println!("\nSz|d>: {}", sigma_z * d);
    println!("\nSz|r>: {}", sigma_z * r);
}
