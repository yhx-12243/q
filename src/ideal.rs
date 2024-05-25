use std::{
    fmt::{Display, Formatter},
    str::FromStr,
};

use num::{
    bigint::{IntDigits, Sign},
    BigInt, BigUint, Integer, One, Signed, Zero,
};
use smallvec::{smallvec_inline, SmallVec};

use crate::{discriminant, factor::factor, pell, qi::QI, qr::quadratic_residue};

#[derive(Clone, Debug, Default, PartialEq, Eq)]
#[repr(transparent)]
pub struct Ideal(SmallVec<[QI; 2]>);

impl Ideal {
    #[inline]
    pub const fn principal(x: QI) -> Self {
        Self(unsafe { SmallVec::from_const_with_len_unchecked([x, QI::ZERO], 1) })
    }

    #[inline]
    pub fn is_principal(&self) -> bool {
        self.0.len() == 1
    }

    pub fn reduce(&mut self) {
        if self.0.len() < 2 {
            return;
        }
        let mut i = self.0.len();
        while i >= 3 {
            let [a, b, c] = unsafe { self.0.get_many_unchecked_mut([i - 3, i - 2, i - 1]) };
            QI::reduce3(a, b, c);
            i -= 1;
        }
        let [a, b] = unsafe { self.0.get_many_unchecked_mut([0, 1]) };
        QI::reduce2(a, b);
        let a = core::mem::take(a);
        let b = core::mem::take(b);
        *self = if b.is_zero() {
            Self::principal(a)
        } else {
            Self(smallvec_inline![a, b])
        }
    }

    pub fn norm(&self) -> BigUint {
        match &*self.0 {
            [q] => q.norm().into_parts().1,
            [q, r] => {
                let a = q.norm();
                let b = r.norm();
                let c = q.dot(r);
                a.magnitude().gcd(b.magnitude()).gcd(c.magnitude())
            }
            _ => panic!("please reduce first before calling norm"),
        }
    }

    pub fn is_multiple_of(&self, other: &Self) -> bool {
        let a = self.norm();
        let b = other.norm();
        if !a.is_multiple_of(&b) {
            return false;
        }
        for qa in &self.0 {
            for qb in &other.0 {
                if !qa.dot(qb).magnitude().is_multiple_of(&b) {
                    return false;
                }
            }
        }
        true
    }

    /// p should be a prime, otherwise UB.
    pub fn factor_prime(p: &BigUint) -> SmallVec<[Self; 2]> {
        let d = discriminant::get();
        let e = d.get() & 3 == 1;
        let q = if e { p * 2u32 } else { p.clone() };

        if e && p.digits() == [2] {
            // special case
            return if d.get() & 7 == 1 {
                // splits
                let ideal = Self(smallvec_inline![
                    BigInt::from(q.clone()).into(),
                    QI { a: BigInt::one(), b: BigInt::one() },
                ]);
                let ideal_ = Self(smallvec_inline![
                    BigInt::from(q).into(),
                    QI { a: -BigInt::one(), b: BigInt::one() },
                ]);
                smallvec_inline![ideal, ideal_]
            } else {
                // inert
                let ideal = Self::principal(BigInt::from(q).into());
                unsafe { SmallVec::from_const_with_len_unchecked([ideal, Self::default()], 1) }
            };
        }

        let mut a = BigUint::from(d.get().unsigned_abs()) % p;
        // ramified
        if a.is_zero() {
            let ideal = Self(smallvec_inline![BigInt::from(q).into(), QI::imaginary()]);
            let ideal_ = ideal.clone();
            return smallvec_inline![ideal, ideal_];
        }

        if d.get() < 0 {
            a = p - a;
        }

        if let Some(mut x) = quadratic_residue(&a, p) {
            // splits
            let one = if e {
                x *= 2u32;
                BigInt::from(2)
            } else {
                BigInt::one()
            };
            let px = BigInt::from(x.clone());
            let nx = BigInt::from_biguint(
                if p.digits() == [2] {
                    Sign::Plus
                } else {
                    Sign::Minus
                },
                x,
            );
            let ideal = Self(smallvec_inline![
                BigInt::from(q.clone()).into(),
                QI { a: px, b: one.clone() }
            ]);
            let ideal_ = Self(smallvec_inline![
                BigInt::from(q).into(),
                QI { a: nx, b: one }
            ]);
            smallvec_inline![ideal, ideal_]
        } else {
            // inert
            let ideal = Self::principal(BigInt::from(q).into());
            unsafe { SmallVec::from_const_with_len_unchecked([ideal, Self::default()], 1) }
        }
    }

    pub fn factor(&mut self) -> anyhow::Result<Vec<(Self, u32)>> {
        self.reduce();

        let d = discriminant::get();
        let e = d.get() & 3 == 1;
        if let [x] = &*self.0 && x.b.is_zero() && x.a.is_positive() {
            let n = x.a.magnitude();
            let factors = if e { factor(&(n / 2u32)) } else { factor(n) }?;
            let mut result = Vec::with_capacity(factors.len() * 2);
            for (p, a) in factors {
                let mut is = Self::factor_prime(&p);
                let num = is.len();
                unsafe { is.set_len(2) };
                let [i1, i2] = is
                    .into_inner()
                    .map_err(|_| anyhow::anyhow!("wrong implementation of Ideal::factor_prime"))?;

                match num {
                    1 => result.push((i1, a)), // principal, skipped
                    2 if i1 == i2 => result.push((
                        if let Some(qi) = pell::work(d.get(), &p) {
                            Self::principal(qi)
                        } else {
                            i1
                        },
                        2 * a,
                    )),
                    2 => {
                        if let Some(qi) = pell::work(d.get(), &p) {
                            let mut qj = qi.clone();
                            qj.a = -qj.a;
                            result.push((Self::principal(qi), a));
                            result.push((Self::principal(qj), a));
                        } else {
                            result.push((i1, a));
                            result.push((i2, a));
                        }
                    }
                    _ => anyhow::bail!("unknown factor_prime error"),
                }
            }
            Ok(result)
        } else {
            unimplemented!()
        }
    }

    pub fn latex(&self, f: &mut Formatter<'_>) {
        let _ = f.write_str("\\left(");
        let mut maybe_comma = "";
        for qi in &self.0 {
            let _ = f.write_str(maybe_comma);
            qi.latex(f);
            maybe_comma = ",";
        }
        let _ = f.write_str("\\right)");
    }
}

impl Display for Ideal {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let mut builder = f.debug_tuple("");
        for qi in &self.0 {
            builder.field_with(|f| qi.fmt(f));
        }
        builder.finish()
    }
}

impl FromStr for Ideal {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> anyhow::Result<Self> {
        s.parse().map(Self::principal)
    }
}

#[cfg(test)]
mod tests {
    use std::num::NonZeroI64;

    use num::BigInt;
    use smallvec::smallvec;

    use super::{Ideal, QI};
    use crate::discriminant;

    fn check(ideal: Ideal, s: &str) {
        let mut t = String::new();
        let mut fmt = core::fmt::Formatter::new(&mut t);
        ideal.latex(&mut fmt);
        assert_eq!(s, t);
    }

    #[test]
    #[rustfmt::skip]
    fn latex_test_1() {
        let _guard = discriminant::DISC_TEST_LOCK.lock().unwrap();
        unsafe { discriminant::set(NonZeroI64::new(-6).unwrap()).unwrap() };

        let ideal = Ideal(smallvec![
            QI::from(BigInt::from(-15)),
            QI::from(BigInt::from(-1)),
            QI::from(BigInt::from(0)),
            QI::from(BigInt::from(1)),
            QI::from(BigInt::from(42)),

            QI { a: 0.into(), b: (-2).into() },
            QI { a: 0.into(), b: (-1).into() },
            QI { a: 0.into(), b: 1.into() },
            QI { a: 0.into(), b: 2.into() },

            QI { a: (-1).into(), b: (-2).into() },
            QI { a: 1.into(), b: (-1).into() },
            QI { a: 2.into(), b: 1.into() },
            QI { a: (-2).into(), b: 2.into() },
        ]);

        check(ideal, r"\left(-15,-1,0,1,42,-2\sqrt{-6},-\sqrt{-6},\sqrt{-6},2\sqrt{-6},-1-2\sqrt{-6},1-\sqrt{-6},2+\sqrt{-6},-2+2\sqrt{-6}\right)");
    }

    #[test]
    #[rustfmt::skip]
    fn latex_test_2() {
        let _guard = discriminant::DISC_TEST_LOCK.lock().unwrap();
        unsafe { discriminant::set(NonZeroI64::new(-7).unwrap()).unwrap() };

        let ideal = Ideal(smallvec![
            QI::from(BigInt::from(-30)),
            QI::from(BigInt::from(-2)),
            QI::from(BigInt::from(0)),
            QI::from(BigInt::from(2)),
            QI::from(BigInt::from(42)),

            QI { a: 0.into(), b: (-4).into() },
            QI { a: 0.into(), b: (-2).into() },
            QI { a: 0.into(), b: 2.into() },
            QI { a: 0.into(), b: 4.into() },

            QI { a: (-2).into(), b: (-4).into() },
            QI { a: 2.into(), b: (-2).into() },
            QI { a: 4.into(), b: 2.into() },
            QI { a: (-4).into(), b: 4.into() },

            QI { a: 3.into(), b: (-3).into() },
            QI { a: (-3).into(), b: (-1).into() },
            QI { a: (-1).into(), b: 1.into() },
            QI { a: 1.into(), b: 3.into() },
        ]);

        check(ideal, r"\left(-15,-1,0,1,21,-2\sqrt{-7},-\sqrt{-7},\sqrt{-7},2\sqrt{-7},-1-2\sqrt{-7},1-\sqrt{-7},2+\sqrt{-7},-2+2\sqrt{-7},\frac{3-3\sqrt{-7}}2,\frac{-3-\sqrt{-7}}2,\frac{-1+\sqrt{-7}}2,\frac{1+3\sqrt{-7}}2\right)");
    }
}
