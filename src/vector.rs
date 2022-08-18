use core::fmt;
use core::marker::PhantomData;

use crate::complex::C64;

/// Marker Trait for Bra and Ket type states.
pub trait BraKet {}
#[derive(Debug)]
pub enum Bra {}
#[derive(Debug)]
pub enum Ket {}
impl BraKet for Bra {}
impl BraKet for Ket {}

/// Vector in complex D-dimensional dual space (populated by bras and kets).
#[derive(Debug)]
pub struct Vector<S: BraKet, const D: usize> {
    inner: [C64; D],
    _s: PhantomData<S>,
}

impl<S: BraKet, const D: usize> Clone for Vector<S, D> {
    fn clone(&self) -> Self {
        Self {
            inner: self.inner,
            _s: PhantomData,
        }
    }
}

impl<S: BraKet, const D: usize> fmt::Display for Vector<S, D> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self.inner)
    }
}

impl<S: BraKet, const D: usize> Vector<S, D> {
    pub fn new(arr: [C64; D]) -> Self {
        Self {
            inner: arr,
            _s: PhantomData,
        }
    }
}

impl<const D: usize> Vector<Ket, D> {
    pub fn inner_prod(&self, other: &Vector<Bra, D>) -> C64 {
        self.inner
            .iter()
            .copied()
            .zip(other.inner.iter().copied())
            .fold(C64::new(0.0, 0.0), |acc, (a, b)| acc + a * b)
    }
    pub fn to_bra(&self) -> Vector<Bra, D> {
        let mut inner: [C64; D] = self.inner;
        inner.iter_mut().for_each(|c| *c = c.conj());
        Vector::<Bra, D> {
            inner,
            _s: PhantomData,
        }
    }
    pub fn to_normalized(&self) -> Self {
        let mut norm = self.clone();
        let ip = norm.inner_prod(&norm.to_bra());
        let magnitude = (ip.real() + ip.imag()).sqrt();
        norm.inner.iter_mut().for_each(|c| *c = *c / magnitude);
        norm
    }
}

impl<const D: usize> Vector<Bra, D> {
    pub fn inner_prod(&self, other: &Vector<Ket, D>) -> C64 {
        self.inner
            .iter()
            .copied()
            .zip(other.inner.iter().copied())
            .fold(C64::new(0.0, 0.0), |acc, (a, b)| acc + a * b)
    }
    pub fn to_ket(&self) -> Vector<Ket, D> {
        let mut inner: [C64; D] = self.inner;
        inner.iter_mut().for_each(|c| *c = c.conj());
        Vector::<Ket, D> {
            inner,
            _s: PhantomData,
        }
    }
    pub fn to_normalized(&self) -> Self {
        let mut norm = self.clone();
        let ip = norm.inner_prod(&norm.to_ket());
        let magnitude = (ip.real() + ip.imag()).sqrt();
        norm.inner.iter_mut().for_each(|c| *c = *c / magnitude);
        norm
    }
}

#[cfg(test)]
mod tests {
    use crate::complex::C64;
    use crate::vector::{Bra, Ket, Vector};

    #[test]
    fn test_bra_ket_round_trip() {
        let ket: Vector<Ket, 3> = Vector::new([
            C64::new(1.3, 0.01),
            C64::new(0.25, 0.75),
            C64::new(0.75, 0.),
        ]);
        let bra = ket.to_bra();
        let ket_reconstituted = bra.to_ket();
        for (vo, vr) in ket
            .inner
            .into_iter()
            .zip(ket_reconstituted.inner.into_iter())
        {
            let diff = vr - vo;
            assert!(diff.real().abs() < 0.0001 && diff.imag().abs() < 0.0001);
        }
    }

    #[test]
    fn test_normalize_ket() {
        let ket: Vector<Ket, 2> = Vector::new([C64::new(1.0, 0.0), C64::new(0.0, 1.0)]);
        let ket_norm = ket.to_normalized();
        let one_over_sqrt2 = 1.0 / f64::sqrt(2.0);
        let elem0_real = ket_norm.inner[0].real();
        let elem0_imag = ket_norm.inner[0].imag();
        assert!((one_over_sqrt2 - elem0_real).abs() < 0.0001);
        assert!(elem0_imag.abs() < 0.0001);

        let elem1_real = ket_norm.inner[1].real();
        let elem1_imag = ket_norm.inner[1].imag();
        assert!(elem1_real.abs() < 0.0001);
        assert!((one_over_sqrt2 - elem1_imag).abs() < 0.0001);
    }

    #[test]
    fn test_normalize_bra() {
        let bra: Vector<Bra, 2> = Vector::new([C64::new(1.0, 0.0), C64::new(0.0, 1.0)]);
        let bra_norm = bra.to_normalized();
        let one_over_sqrt2 = 1.0 / f64::sqrt(2.0);
        let elem0_real = bra_norm.inner[0].real();
        let elem0_imag = bra_norm.inner[0].imag();
        assert!((one_over_sqrt2 - elem0_real).abs() < 0.0001);
        assert!(elem0_imag.abs() < 0.0001);

        let elem1_real = bra_norm.inner[1].real();
        let elem1_imag = bra_norm.inner[1].imag();
        assert!(elem1_real.abs() < 0.0001);
        assert!((one_over_sqrt2 - elem1_imag).abs() < 0.0001);
    }
}
