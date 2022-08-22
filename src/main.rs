mod complex;
mod operator;
mod vector;

use complex::C64;
use operator::HermitianMatrix;
use vector::{Bra, InnerProductDualSpace, Ket, Vector};

fn main() {
    let u: Vector<Ket, 2> = Vector::from_arr([C64::new(1.0, 0.0), C64::new(0.0, 0.0)]);
    let d: Vector<Ket, 2> = Vector::from_arr([C64::new(0.0, 0.0), C64::new(1.0, 0.0)]);

    println!("\nUp spin state vector {}", &u);
    println!("\nDown spin state vector {}", &d);

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

    println!("\nSpin operator (z): {:?}", &sigma_z);
    println!("\nSpin operator (x): {:?}", &sigma_x);
    println!("\nSpin operator (y): {:?}", &sigma_y);

    println!("\nSz|u> = {}", &sigma_z * &u);
    println!("\nSz|d> = {}", &sigma_z * &d);
}
