use core::fmt;
use core::iter::IntoIterator;
use core::marker::PhantomData;
use core::ops::{Index, IndexMut, Mul};
use core::slice::SliceIndex;

use crate::complex::C64;
use crate::operator::HermitianMatrix;

/// A dual (complex-valued) inner product space.
pub trait InnerProductDualSpace {
    type Dual;
    type Scalar;

    fn to_dual(&self) -> Self::Dual;
    fn inner_product(&self, dual: &Self::Dual) -> Self::Scalar;
    fn normalize(&mut self);
}

/// Marker Trait for Bra and Ket type states.
pub trait BraKet {}
/// Marker for a "bra" in the bra-ket notation.
#[derive(Debug)]
pub enum Bra {}
/// Marker for a "ket" in the bra-ket notation.
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

impl<S: BraKet, const D: usize> Vector<S, D> {
    pub fn new() -> Self {
        Self {
            inner: [C64::default(); D],
            _s: PhantomData,
        }
    }
    pub fn from_arr(arr: [C64; D]) -> Self {
        Self {
            inner: arr,
            _s: PhantomData,
        }
    }
    pub fn iter(&self) -> core::slice::Iter<'_, C64> {
        self.into_iter()
    }
    pub fn iter_mut(&mut self) -> core::slice::IterMut<'_, C64> {
        self.into_iter()
    }
}

impl<S: BraKet, const D: usize> IntoIterator for Vector<S, D> {
    type Item = C64;
    type IntoIter = core::array::IntoIter<C64, D>;

    fn into_iter(self) -> Self::IntoIter {
        self.inner.into_iter()
    }
}

impl<'a, S: BraKet, const D: usize> IntoIterator for &'a Vector<S, D> {
    type Item = &'a C64;
    type IntoIter = core::slice::Iter<'a, C64>;

    fn into_iter(self) -> Self::IntoIter {
        self.inner.iter()
    }
}

impl<'a, S: BraKet, const D: usize> IntoIterator for &'a mut Vector<S, D> {
    type Item = &'a mut C64;
    type IntoIter = core::slice::IterMut<'a, C64>;

    fn into_iter(self) -> Self::IntoIter {
        self.inner.iter_mut()
    }
}

impl<S: BraKet, const D: usize, Idx> Index<Idx> for Vector<S, D>
where
    Idx: SliceIndex<[C64], Output = C64>,
{
    type Output = C64;

    #[inline(always)]
    fn index(&self, index: Idx) -> &Self::Output {
        self.inner.index(index)
    }
}

impl<S: BraKet, const D: usize, Idx> IndexMut<Idx> for Vector<S, D>
where
    Idx: SliceIndex<[C64], Output = C64>,
{
    #[inline(always)]
    fn index_mut(&mut self, index: Idx) -> &mut Self::Output {
        self.inner.index_mut(index)
    }
}

impl<S: BraKet, const D: usize> Clone for Vector<S, D> {
    fn clone(&self) -> Self {
        Self {
            inner: self.inner,
            _s: PhantomData,
        }
    }
}

impl<S: BraKet, const D: usize> Default for Vector<S, D> {
    fn default() -> Self {
        Self::new()
    }
}

impl<const D: usize> fmt::Display for Vector<Bra, D> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "(<|) [")?;
        for idx in 0..D - 1 {
            let v = self[idx];
            let join_op = if v.imag() >= 0.0 { "+" } else { "-" };
            write!(f, "{} {} {}i, ", v.real(), join_op, v.imag())?;
        }
        let v = self[D - 1];
        let join_op = if v.imag() >= 0.0 { "+" } else { "-" };
        write!(f, "{} {} {}i", v.real(), join_op, v.imag())?;
        write!(f, "]")
    }
}

impl<const D: usize> Vector<Bra, D> {
    pub fn to_ket(&self) -> Vector<Ket, D> {
        let mut inner: [C64; D] = self.inner;
        inner.iter_mut().for_each(|c| *c = c.conj());
        Vector::<Ket, D> {
            inner,
            _s: PhantomData,
        }
    }
}

impl<const D: usize> InnerProductDualSpace for Vector<Bra, D> {
    type Dual = Vector<Ket, D>;
    type Scalar = C64;

    fn to_dual(&self) -> Self::Dual {
        self.to_ket()
    }
    fn inner_product(&self, dual: &Self::Dual) -> Self::Scalar {
        self * dual
    }
    fn normalize(&mut self) {
        let ip = self.inner_product(&(&*self).to_ket());
        let magnitude = (ip.real() + ip.imag()).sqrt();
        self.iter_mut().for_each(|c| *c /= magnitude);
    }
}

impl<const D: usize> Mul<&HermitianMatrix<D>> for &Vector<Bra, D> {
    type Output = Vector<Bra, D>;

    fn mul(self, rhs: &HermitianMatrix<D>) -> Vector<Bra, D> {
        let mut out_bra: Vector<Bra, D> = Vector::default();
        for ridx in 0..D {
            let b = self[ridx];
            for cidx in 0..D {
                out_bra[cidx] = b * rhs.inner[ridx][cidx];
            }
        }
        out_bra
    }
}

impl<const D: usize> Mul<&Vector<Ket, D>> for &Vector<Bra, D> {
    type Output = C64;

    fn mul(self, rhs: &Vector<Ket, D>) -> C64 {
        self.iter()
            .copied()
            .zip(rhs.iter().copied())
            .fold(C64::new(0.0, 0.0), |acc, (a, b)| acc + a * b)
    }
}

impl<const D: usize> fmt::Display for Vector<Ket, D> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "(|>) [")?;
        for idx in 0..D - 1 {
            let v = self[idx];
            let join_op = if v.imag() >= 0.0 { "+" } else { "-" };
            write!(f, "{} {} {}i, ", v.real(), join_op, v.imag())?;
        }
        let v = self[D - 1];
        let join_op = if v.imag() >= 0.0 { "+" } else { "-" };
        write!(f, "{} {} {}i", v.real(), join_op, v.imag())?;
        write!(f, "]")
    }
}

impl<const D: usize> Vector<Ket, D> {
    pub fn to_bra(&self) -> Vector<Bra, D> {
        let mut inner: [C64; D] = self.inner;
        inner.iter_mut().for_each(|c| *c = c.conj());
        Vector::<Bra, D> {
            inner,
            _s: PhantomData,
        }
    }
}

impl<const D: usize> InnerProductDualSpace for Vector<Ket, D> {
    type Dual = Vector<Bra, D>;
    type Scalar = C64;

    fn to_dual(&self) -> Self::Dual {
        self.to_bra()
    }
    fn inner_product(&self, dual: &Self::Dual) -> Self::Scalar {
        dual * self
    }
    fn normalize(&mut self) {
        let ip = self.inner_product(&(&*self).to_bra());
        let magnitude = (ip.real() + ip.imag()).sqrt();
        self.iter_mut().for_each(|c| *c /= magnitude);
    }
}

#[cfg(test)]
mod tests {
    use crate::complex::C64;
    use crate::vector::{Bra, InnerProductDualSpace, Ket, Vector};

    #[test]
    fn test_bra_ket_round_trip() {
        let ket: Vector<Ket, 3> = Vector::from_arr([
            C64::new(1.3, 0.01),
            C64::new(0.25, 0.75),
            C64::new(0.75, 0.),
        ]);
        let bra = ket.to_bra();
        let ket_reconstituted = bra.to_ket();
        for (vo, vr) in ket.into_iter().zip(ket_reconstituted.into_iter()) {
            let diff = vr - vo;
            assert!(diff.real().abs() < 0.0001 && diff.imag().abs() < 0.0001);
        }
    }

    #[test]
    fn test_normalize_ket() {
        let ket: Vector<Ket, 2> = Vector::from_arr([C64::new(1.0, 0.0), C64::new(0.0, 1.0)]);
        let mut ket_norm = ket.clone();
        ket_norm.normalize();
        let one_over_sqrt2 = 1.0 / f64::sqrt(2.0);
        let elem0_real = ket_norm[0].real();
        let elem0_imag = ket_norm[0].imag();
        assert!((one_over_sqrt2 - elem0_real).abs() < 0.0001);
        assert!(elem0_imag.abs() < 0.0001);

        let elem1_real = ket_norm[1].real();
        let elem1_imag = ket_norm[1].imag();
        assert!(elem1_real.abs() < 0.0001);
        assert!((one_over_sqrt2 - elem1_imag).abs() < 0.0001);
    }

    #[test]
    fn test_normalize_bra() {
        let bra: Vector<Bra, 2> = Vector::from_arr([C64::new(1.0, 0.0), C64::new(0.0, 1.0)]);
        let mut bra_norm = bra.clone();
        bra_norm.normalize();
        let one_over_sqrt2 = 1.0 / f64::sqrt(2.0);
        let elem0_real = bra_norm[0].real();
        let elem0_imag = bra_norm[0].imag();
        assert!((one_over_sqrt2 - elem0_real).abs() < 0.0001);
        assert!(elem0_imag.abs() < 0.0001);

        let elem1_real = bra_norm[1].real();
        let elem1_imag = bra_norm[1].imag();
        assert!(elem1_real.abs() < 0.0001);
        assert!((one_over_sqrt2 - elem1_imag).abs() < 0.0001);
    }
}
