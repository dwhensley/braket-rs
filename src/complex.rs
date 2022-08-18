use core::ops::{Add, Div, Mul, Sub};

/// Complex number with each component represented with `f64`s.
#[derive(Debug, Copy, Clone)]
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
}

impl Default for C64 {
    fn default() -> Self {
        Self { re: 0.0, im: 0.0 }
    }
}

impl Add for C64 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            re: self.re + other.re,
            im: self.im + other.im,
        }
    }
}

impl Sub for C64 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            re: self.re - other.re,
            im: self.im - other.im,
        }
    }
}

impl Mul<C64> for C64 {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        Self {
            re: self.re * other.re - self.im * other.im,
            im: self.re * other.im + self.im * other.re,
        }
    }
}

impl Mul<f64> for C64 {
    type Output = Self;

    fn mul(self, other: f64) -> Self {
        Self {
            re: self.re * other,
            im: self.im * other,
        }
    }
}

impl Div<C64> for C64 {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        let denom = self.re * self.re + self.im * self.im;
        let re = (other.re * self.re + other.im * self.im) / denom;
        let im = (other.im * self.re - other.re * self.im) / denom;
        Self { re, im }
    }
}

impl Div<f64> for C64 {
    type Output = Self;

    fn div(self, other: f64) -> Self {
        Self {
            re: self.re / other,
            im: self.im / other,
        }
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
