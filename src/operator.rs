use core::fmt;
use core::ops::Mul;

use crate::complex::C64;
use crate::vector::{Ket, Vector};

#[derive(Debug)]
pub enum OperatorError {
    HermitianPropertiesNotSatisfied,
}

impl fmt::Display for OperatorError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let err_msg = match self {
            OperatorError::HermitianPropertiesNotSatisfied => {
                "Attempt to contruct non-Hermitian operator"
            }
        };
        write!(f, "{}", err_msg)
    }
}

/// DxD Hermitian operator.
#[derive(Debug, Copy, Clone)]
pub struct HermitianMatrix<const D: usize> {
    pub(crate) inner: [[C64; D]; D],
}

impl<const D: usize> HermitianMatrix<D> {
    pub fn from_arr(arr: [[C64; D]; D]) -> Result<Self, OperatorError> {
        for ridx in 0..D {
            for cidx in 0..D {
                if arr[ridx][cidx] != arr[cidx][ridx].conj() {
                    return Err(OperatorError::HermitianPropertiesNotSatisfied);
                }
            }
        }
        Ok(Self { inner: arr })
    }
}

impl<const D: usize> fmt::Display for HermitianMatrix<D> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "\n[")?;
        for ridx in 0..D {
            if ridx > 0 {
                write!(f, " [")?;
            } else {
                write!(f, "[")?;
            }
            for cidx in 0..D - 1 {
                write!(f, "{}, ", self.inner[ridx][cidx])?;
            }
            if ridx < D - 1 {
                writeln!(f, "{}]", self.inner[ridx][D - 1])?;
            } else {
                write!(f, "{}]]", self.inner[ridx][D - 1])?;
            }
        }
        Ok(())
    }
}

impl<const D: usize> Mul<Vector<Ket, D>> for HermitianMatrix<D> {
    type Output = Vector<Ket, D>;

    fn mul(self, rhs: Vector<Ket, D>) -> Vector<Ket, D> {
        let mut out_ket: Vector<Ket, D> = Vector::default();
        for ridx in 0..D {
            let mut out = C64::zero();
            for (m, v) in self.inner[ridx].into_iter().zip(rhs.into_iter()) {
                out += m * v;
            }
            out_ket[ridx] = out;
        }
        out_ket
    }
}

#[cfg(test)]
mod tests {
    use crate::complex::C64;
    use crate::operator::HermitianMatrix;
    use crate::vector::{InnerProductDualSpace, Ket, Vector};

    #[test]
    fn test_hm_identity_on_ket() {
        let ket: Vector<Ket, 3> = Vector::from_arr([
            C64::new(1.3, 0.01),
            C64::new(0.25, 0.75),
            C64::new(0.3, 0.7),
        ]);
        let zero = C64::zero();
        let one = C64::new(1.0, 0.0);
        let identity_op: HermitianMatrix<3> =
            HermitianMatrix::from_arr([[one, zero, zero], [zero, one, zero], [zero, zero, one]])
                .unwrap();
        let out = &identity_op * &ket;
        for (vi, vo) in ket.into_iter().zip(out.into_iter()) {
            let diff = vo - vi;
            assert!(diff.real().abs() < 0.0001 && diff.imag().abs() < 0.0001);
        }
    }

    #[test]
    fn test_cannot_construct_non_hermitian_operator() {
        let op_result = HermitianMatrix::<2>::from_arr([
            [C64::new(1.0, 0.0), C64::new(0.25, 0.75)],
            [C64::new(0.75, 0.25), C64::new(1.0, 0.0)],
        ]);
        assert!(op_result.is_err());
    }
}
