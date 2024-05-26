use std::fmt::{Display, Formatter, Write};

use num::{bigint::Sign, BigInt, BigUint, Integer, One, Signed, Zero};

use crate::discriminant;

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct QI {
    pub a: BigInt,
    pub b: BigInt,
}

impl QI {
    pub const ZERO: Self = Self {
        a: BigInt::ZERO,
        b: BigInt::ZERO,
    };

    #[inline]
    pub fn is_zero(&self) -> bool {
        self.a.is_zero() && self.b.is_zero()
    }

    #[inline]
    pub fn imaginary() -> Self {
        let e = discriminant::is4kp1();
        Self {
            a: BigInt::ZERO,
            b: if e { BigInt::from(2) } else { BigInt::one() },
        }
    }

    pub fn reduce2(a: &mut Self, b: &mut Self) -> bool {
        if a.is_zero() {
            core::mem::swap(a, b);
            return false;
        }
        if b.is_zero() {
            return false;
        }
        let mut na = a.gaussian_norm();
        let mut nb = b.gaussian_norm();
        if a.collinear_with(b) {
            let ca = a.a.gcd(&b.a);
            let mut cb = a.b.gcd(&b.b);
            if a.a.sign() * b.a.sign() == Sign::Minus {
                cb = -cb;
            }
            a.a = ca;
            a.b = cb;
            *b = Self::ZERO;
            return true;
        }
        let mut dot = a.gaussian_dot(b);
        let mut changed = false;
        while !dot.is_zero() {
            let mut progress = false;
            let qb = BigInt::from_biguint(dot.sign(), (dot.magnitude() + (&nb - 1u32) / 2u32) / &nb);
            if !qb.is_zero() {
                progress = true;
                changed = true;
                a.a -= &b.a * &qb;
                a.b -= &b.b * &qb;
                na = a.gaussian_norm();
                dot -= BigInt::from_biguint(qb.sign(), &nb * qb.magnitude());
            }
            let qa = BigInt::from_biguint(dot.sign(), (dot.magnitude() + (&na - 1u32) / 2u32) / &na);
            if !qa.is_zero() {
                progress = true;
                changed = true;
                b.a -= &a.a * &qa;
                b.b -= &a.b * &qa;
                nb = b.gaussian_norm();
                dot -= BigInt::from_biguint(qa.sign(), &na * qa.magnitude());
            }
            if !progress {
                break;
            }
        }
        if b.is_zero() {
            core::mem::swap(a, b);
        }
        changed
    }

    pub fn reduce3(a: &mut Self, b: &mut Self, c: &mut Self) {
        loop {
            let ab = Self::reduce2(a, b);
            let bc = Self::reduce2(b, c);
            let ca = Self::reduce2(c, a);
            if !(ab || bc || ca) {
                if !c.is_zero() {
                    if b.is_zero() {
                        core::mem::swap(b, c);
                    } else if a.is_zero() {
                        core::mem::swap(a, c);
                    }
                }
                assert!(c.is_zero());
                if a.is_zero() {
                    core::mem::swap(a, b);
                }
                return;
            }
        }
    }

    #[inline]
    fn collinear_with(&self, other: &Self) -> bool {
        &self.a * &other.b == &self.b * &other.a
    }

    #[inline]
    fn gaussian_norm(&self) -> BigUint {
        let ma = self.a.magnitude();
        let mb = self.b.magnitude();
        ma * ma + mb * mb
    }

    #[inline]
    fn gaussian_dot(&self, other: &Self) -> BigInt {
        &self.a * &other.a + &self.b * &other.b
    }

    #[inline]
    pub fn norm(&self) -> BigInt {
        let d = discriminant::get();
        let pre_norm = &self.a * &self.a - &self.b * &self.b * d.get();
        if d.get() & 3 == 1 {
            pre_norm / 4
        } else {
            pre_norm
        }
    }

    pub fn common(&self) -> BigUint {
        let e = discriminant::is4kp1();
        let common = self.a.magnitude().gcd(self.b.magnitude());
        if e && common.is_even() && self.a.trailing_zeros() != self.b.trailing_zeros() {
            common / 2u32
        } else {
            common
        }
    }

    #[inline]
    pub fn dot(&self, other: &Self) -> BigInt {
        let d = discriminant::get();
        let pre_dot = &self.a * &other.a - &self.b * &other.b * d.get();
        if d.get() & 3 == 1 {
            pre_dot / 2
        } else {
            pre_dot * 2
        }
    }

    pub fn latex(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let e = discriminant::is4kp1();
        let s = discriminant::get_latex();
        if self.b.is_zero() {
            if e {
                let a: BigInt = &self.a / 2;
                a.fmt(f)
            } else {
                self.a.fmt(f)
            }
        } else if self.a.is_zero() {
            let b_: BigInt;
            let b = if e {
                b_ = &self.b / 2;
                &b_
            } else {
                &self.b
            };
            if !b.magnitude().is_one() {
                b.fmt(f)?;
            } else if b.is_negative() {
                f.write_char('-')?;
            }
            f.write_str(s)
        } else {
            let div2 = e && self.a.is_odd();
            if div2 {
                f.write_str("\\frac{")?;
            }
            let (a_, b_): (BigInt, BigInt);
            let (a, b) = if e && !self.a.is_odd() {
                a_ = &self.a / 2;
                b_ = &self.b / 2;
                (&a_, &b_)
            } else {
                (&self.a, &self.b)
            };
            a.fmt(f)?;
            if b.magnitude().is_one() {
                f.write_char(if b.is_negative() { '-' } else { '+' })?;
            } else {
                write!(f, "{b:+}")?;
            }
            f.write_str(s)?;
            if div2 {
                f.write_str("}2")?;
            }
            Ok(())
        }
    }
}

impl Display for QI {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let e = discriminant::is4kp1();
        let s = discriminant::get_str();
        if !e {
            write!(f, "{}{:+}{s}", self.a, self.b)?;
        } else if self.a.is_odd() {
            write!(f, "{}.5{:+}.5{s}", &self.a / 2, &self.b / 2)?;
        } else {
            write!(f, "{}{:+}{s}", &self.a / 2, &self.b / 2)?;
        }
        Ok(())
    }
}

impl From<BigInt> for QI {
    fn from(n: BigInt) -> Self {
        Self {
            a: n,
            b: BigInt::ZERO,
        }
    }
}

#[cfg(test)]
mod tests {
    use std::num::NonZeroI64;

    use num::BigInt;

    use super::{discriminant, QI};

    fn check(qi: QI, s: &str) {
        let mut t = String::new();
        let mut fmt = core::fmt::Formatter::new(&mut t);
        qi.latex(&mut fmt);
        assert_eq!(s, t);
    }

    #[test]
    #[rustfmt::skip]
    fn latex_test_1() {
        let _guard = discriminant::DISC_TEST_LOCK.lock().unwrap();
        unsafe { discriminant::set(NonZeroI64::new(-6).unwrap()).unwrap() };

        check(QI::from(BigInt::from(-15)), "-15");
        check(QI::from(BigInt::from(-1)), "-1");
        check(QI::from(BigInt::from(0)), "0");
        check(QI::from(BigInt::from(1)), "1");
        check(QI::from(BigInt::from(42)), "42");

        check(QI { a: 0.into(), b: (-2).into() }, r"-2\sqrt{-6}");
        check(QI { a: 0.into(), b: (-1).into() }, r"-\sqrt{-6}");
        check(QI { a: 0.into(), b: 1.into() }, r"\sqrt{-6}");
        check(QI { a: 0.into(), b: 2.into() }, r"2\sqrt{-6}");

        check(QI { a: (-1).into(), b: (-2).into() }, r"-1-2\sqrt{-6}");
        check(QI { a: 1.into(), b: (-1).into() }, r"1-\sqrt{-6}");
        check(QI { a: 2.into(), b: 1.into() }, r"2+\sqrt{-6}");
        check(QI { a: (-2).into(), b: 2.into() }, r"-2+2\sqrt{-6}");
    }

    #[test]
    #[rustfmt::skip]
    fn latex_test_2() {
        let _guard = discriminant::DISC_TEST_LOCK.lock().unwrap();
        unsafe { discriminant::set(NonZeroI64::new(-7).unwrap()).unwrap() };

        check(QI::from(BigInt::from(-30)), "-15");
        check(QI::from(BigInt::from(-2)), "-1");
        check(QI::from(BigInt::from(0)), "0");
        check(QI::from(BigInt::from(2)), "1");
        check(QI::from(BigInt::from(42)), "21");

        check(QI { a: 0.into(), b: (-4).into() }, r"-2\sqrt{-7}");
        check(QI { a: 0.into(), b: (-2).into() }, r"-\sqrt{-7}" );
        check(QI { a: 0.into(), b: 2.into() }, r"\sqrt{-7}" );
        check(QI { a: 0.into(), b: 4.into() }, r"2\sqrt{-7}" );

        check(QI { a: (-2).into(), b: (-4).into() }, r"-1-2\sqrt{-7}" );
        check(QI { a: 2.into(), b: (-2).into() }, r"1-\sqrt{-7}" );
        check(QI { a: 4.into(), b: 2.into() }, r"2+\sqrt{-7}" );
        check(QI { a: (-4).into(), b: 4.into() }, r"-2+2\sqrt{-7}" );

        check(QI { a: 3.into(), b: (-3).into() }, r"\frac{3-3\sqrt{-7}}2" );
        check(QI { a: (-3).into(), b: (-1).into() }, r"\frac{-3-\sqrt{-7}}2" );
        check(QI { a: (-1).into(), b: 1.into() }, r"\frac{-1+\sqrt{-7}}2" );
        check(QI { a: 1.into(), b: 3.into() }, r"\frac{1+3\sqrt{-7}}2" );
    }
}
