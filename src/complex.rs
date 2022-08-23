use core::fmt;
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub, SubAssign};

/// Complex number with each component represented with `f64`s.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct C64 {
    re: f64,
    im: f64,
}

impl C64 {
    pub fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }
    pub fn conj(self) -> Self {
        Self {
            re: self.re,
            im: -self.im,
        }
    }
    pub fn to_polar(self) -> (f64, f64) {
        let r = (self.re * self.re + self.im * self.im).sqrt();
        let theta = f64::atan2(self.im, self.re);
        (r, theta)
    }
    pub fn from_polar(r: f64, theta: f64) -> Self {
        let (re, im) = (r * f64::cos(theta), r * f64::sin(theta));
        Self { re, im }
    }
    pub fn real(&self) -> f64 {
        self.re
    }
    pub fn imag(&self) -> f64 {
        self.im
    }
    pub const fn zero() -> Self {
        Self { re: 0.0, im: 0.0 }
    }
    pub const fn one() -> Self {
        Self { re: 1.0, im: 0.0 }
    }
    pub const fn i() -> Self {
        Self { re: 0.0, im: 1.0 }
    }
}

impl Default for C64 {
    fn default() -> Self {
        Self::zero()
    }
}

impl fmt::Display for C64 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let join_op = if self.imag() >= 0.0 { "+" } else { "-" };
        write!(f, "{} {} {}i", self.real(), join_op, self.imag().abs())
    }
}

impl Add for C64 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self {
            re: self.re + rhs.re,
            im: self.im + rhs.im,
        }
    }
}

impl AddAssign for C64 {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl Sub for C64 {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self {
            re: self.re - rhs.re,
            im: self.im - rhs.im,
        }
    }
}

impl SubAssign for C64 {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl Mul<C64> for C64 {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self {
            re: self.re * rhs.re - self.im * rhs.im,
            im: self.re * rhs.im + self.im * rhs.re,
        }
    }
}

impl MulAssign<C64> for C64 {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl Mul<f64> for C64 {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self {
        Self {
            re: self.re * rhs,
            im: self.im * rhs,
        }
    }
}

impl MulAssign<f64> for C64 {
    fn mul_assign(&mut self, rhs: f64) {
        *self = *self * rhs;
    }
}

impl Div<C64> for C64 {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        let denom = self.re * self.re + self.im * self.im;
        let re = (rhs.re * self.re + rhs.im * self.im) / denom;
        let im = (rhs.im * self.re - rhs.re * self.im) / denom;
        Self { re, im }
    }
}

impl DivAssign<C64> for C64 {
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs;
    }
}

impl Div<f64> for C64 {
    type Output = Self;

    fn div(self, rhs: f64) -> Self {
        Self {
            re: self.re / rhs,
            im: self.im / rhs,
        }
    }
}

impl DivAssign<f64> for C64 {
    fn div_assign(&mut self, rhs: f64) {
        *self = *self / rhs;
    }
}

#[cfg(test)]
mod tests {
    use crate::C64;

    #[test]
    fn test_complex_cartesian_polar_round_trip() {
        let c = C64::new(0.3, -0.65);
        let (r, theta) = c.to_polar();
        let c_reconstituted = C64::from_polar(r, theta);
        let diff = c_reconstituted - c;
        assert!(diff.real().abs() < 0.0001 && diff.imag().abs() < 0.0001);
    }
}
