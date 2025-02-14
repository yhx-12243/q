use core::{
    cmp::Ordering,
    fmt::{self, Formatter},
    mem::MaybeUninit,
    str::FromStr,
};
use std::io;

use num::{
    BigInt, BigUint, Integer, One, Zero,
    bigint::{IntDigits, Sign},
};
use smallvec::{SmallVec, smallvec_inline};

use crate::{discriminant, factor::factor, pell, qi::QI, qr::quadratic_residue};

#[derive(Clone, Debug, Default, PartialEq, Eq)]
#[repr(transparent)]
pub struct Ideal(SmallVec<QI, 2>);

fn extract<R>(input: &mut R) -> io::Result<Option<BigInt>>
where
    R: io::BufRead,
{
    #[inline]
    const fn check(x: u8) -> bool {
        x.is_ascii_digit() || x == b'-'
    }

    #[inline]
    fn finalize(data: &[u8]) -> io::Result<Option<BigInt>> {
        let data = unsafe { core::str::from_utf8_unchecked(data) };
        match BigInt::from_str(data) {
            Ok(i) => Ok(Some(i)),
            Err(e) => Err(io::Error::other(e)),
        }
    }

    let mut buf = Vec::new();

    loop {
        let available = match input.fill_buf() {
            Ok(n) => n,
            Err(ref e) if e.is_interrupted() => continue,
            Err(e) => return Err(e),
        };
        if buf.is_empty() {
            if let Some(pos) = available.iter().position(|x| check(*x)) {
                let next = &available[pos..];
                if let Some(right) = next.iter().position(|x| !check(*x)) {
                    let ret = finalize(&next[..right]);
                    input.consume(pos + right);
                    return ret;
                }
                buf.extend_from_slice(next);
            }
        } else {
            if !available.first().is_some_and(|x| check(*x)) {
                return finalize(&buf);
            }
            if let Some(right) = available.iter().position(|x| !check(*x)) {
                buf.extend_from_slice(&available[..right]);
                input.consume(right);
                return finalize(&buf);
            }
            buf.extend_from_slice(available);
        }
        if available.is_empty() {
            return Ok(None);
        }
        let len = available.len();
        input.consume(len);
    }
}

impl Ideal {
    #[inline]
    pub const fn zero() -> Self {
        Self(SmallVec::new())
    }

    #[inline]
    pub const fn principal(x: QI) -> Self {
        Self(unsafe { SmallVec::from_buf_and_len_unchecked(MaybeUninit::new([x, QI::ZERO]), 1) })
    }

    #[inline]
    pub fn is_zero(&self) -> bool {
        self.0.is_empty()
    }

    pub fn reduce(&mut self) {
        if self.0.len() < 2 {
            return;
        }
        let mut i = self.0.len();
        while i >= 3 {
            let [a, b, c] = unsafe { self.0.get_disjoint_unchecked_mut([i - 3, i - 2, i - 1]) };
            QI::reduce3(a, b, c);
            i -= 1;
        }
        let [a, b] = unsafe { self.0.get_disjoint_unchecked_mut([0, 1]) };
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
    pub fn factor_prime(p: &BigUint) -> SmallVec<Self, 2> {
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
                SmallVec::from_buf_and_len([ideal, const { Self::zero() }], 1)
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
                QI { a: px, b: one.clone() },
            ]);
            let ideal_ = Self(smallvec_inline![
                BigInt::from(q).into(),
                QI { a: nx, b: one },
            ]);
            smallvec_inline![ideal, ideal_]
        } else {
            // inert
            let ideal = Self::principal(BigInt::from(q).into());
            SmallVec::from_buf_and_len([ideal, const { Self::zero() }], 1)
        }
    }

    #[allow(clippy::too_many_lines)]
    pub fn factor(&mut self) -> anyhow::Result<Vec<(Self, u32)>> {
        self.reduce();

        let mut common = match &mut *self.0 {
            [] => return Ok(vec![(Self::default(), 1)]), // zero ideal
            [q] => q.common(),
            [q, r] => q.common().gcd(&r.common()),
            // SAFETY: we already reduced before
            _ => unsafe { core::hint::unreachable_unchecked() },
        };

        {
            let common_ = BigInt::from(common);
            for each in &mut self.0 {
                each.a /= &common_;
                each.b /= &common_;
            }
            common = common_.into_parts().1;
        }

        let d = discriminant::get();

        let norm = self.norm();

        // eprintln!("common = {common}, remain = {self}, norm = {norm}");

        let [common_factors, norm_factors] = factor([&common, &norm])?;
        let mut common_factors = common_factors.iter().peekable();
        let mut norm_factors = norm_factors.iter().peekable();

        let mut result = Vec::with_capacity((common_factors.len() + norm_factors.len()) * 2);
        loop {
            let l = common_factors.peek();
            let r = norm_factors.peek();
            let (p, a1, a2) = if let Some(l) = l && let Some(r) = r {
                match l.0.cmp(&r.0) {
                    Ordering::Equal => {
                        let (p, a1) = unsafe { common_factors.next().unwrap_unchecked() };
                        let (_, a2) = unsafe { norm_factors.next().unwrap_unchecked() };
                        (p, *a1, *a2)
                    }
                    Ordering::Less => {
                        let (p, a) = unsafe { common_factors.next().unwrap_unchecked() };
                        (p, *a, 0)
                    }
                    Ordering::Greater => {
                        let (p, a) = unsafe { norm_factors.next().unwrap_unchecked() };
                        (p, 0, *a)
                    }
                }
            } else if let Some((p, a)) = common_factors.next() {
                (p, *a, 0)
            } else if let Some((p, a)) = norm_factors.next() {
                (p, 0, *a)
            } else {
                break;
            };

            let mut is = Self::factor_prime(p);
            let num = is.len();
            unsafe { is.set_len(2) };
            let [mut i1, mut i2] = is
                .into_inner()
                .map_err(|_| anyhow::anyhow!("wrong implementation of Ideal::factor_prime"))?;

            // eprintln!("processing {p} {a1} {a2}");

            match num {
                1 => {
                    if a2 != 0 {
                        anyhow::bail!("contradiction on prime {p} (kind 1)");
                    }
                    // principal, skipped
                    result.push((i1, a1));
                }
                2 if i1 == i2 => {
                    if a2 > 1 {
                        anyhow::bail!("contradiction on prime {p} (kind 2)");
                    }
                    if let Some(qi) = pell::work(d.get(), p) {
                        i1 = Self::principal(qi);
                    }
                    result.push((i1, 2 * a1 + a2));
                }
                2 => {
                    if let Some(qi) = pell::work(d.get(), p) {
                        let mut qj = qi.clone();
                        qj.a = -qj.a;
                        i1 = Self::principal(qi);
                        i2 = Self::principal(qj);
                    }
                    if a2 == 0 {
                        result.push((i1, a1));
                        result.push((i2, a1));
                    } else {
                        match (self.is_multiple_of(&i1), self.is_multiple_of(&i2)) {
                            (true, false) => {
                                result.push((i1, a1 + a2));
                                if a1 != 0 { result.push((i2, a1)); }
                            }
                            (false, true) => {
                                if a1 != 0 { result.push((i1, a1)); }
                                result.push((i2, a1 + a2));
                            }
                            (x, y) => anyhow::bail!("contradiction on prime {p} (kind 3) : ({x}, {y})"),
                        }
                    }
                }
                _ => anyhow::bail!("unknown factor_prime error"),
            }
        }
        Ok(result)
    }

    fn tex_common(&self, f: &mut Formatter<'_>, inner: fn(&QI, &mut Formatter) -> fmt::Result) -> fmt::Result {
        f.write_str("\\left(")?;
        if self.is_zero() {
            f.write_str("0")?;
        } else {
            let mut maybe_comma = "";
            for qi in &self.0 {
                f.write_str(maybe_comma)?;
                inner(qi, f)?;
                maybe_comma = ",";
            }
        }
        f.write_str("\\right)")
    }

    pub fn latex(&self, f: &mut Formatter<'_>) -> fmt::Result {
        self.tex_common(f, QI::latex)
    }

    pub fn tex(&self, f: &mut Formatter<'_>) -> fmt::Result {
        self.tex_common(f, QI::tex)
    }

    pub fn read<R>(mut input: R) -> io::Result<Self>
    where
        R: io::BufRead,
    {
        let e = discriminant::is4kp1();
        let mut qis = const { SmallVec::new() };
        let mut last: Option<BigInt> = None;

        while let Some(cur) = extract(&mut input)? {
            if let Some(last) = last.take() {
                if e && last.is_even() ^ cur.is_even() {
                    return Err(io::Error::other("not a integer".to_owned()));
                }
                if !(last.is_zero() && cur.is_zero()) {
                    qis.push(QI { a: last, b: cur });
                }
            } else {
                last = Some(cur);
            }
        }

        if last.is_some() {
            Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "odd number of integers were inputted".to_owned(),
            ))
        } else {
            Ok(Self(qis))
        }
    }
}

impl fmt::Display for Ideal {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let mut builder = f.debug_tuple("");
        for qi in &self.0 {
            builder.field_with(|f| qi.fmt(f));
        }
        builder.finish()
    }
}

#[cfg(test)]
mod tests {
    use core::num::NonZeroI64;

    use num::BigInt;
    use smallvec::smallvec;

    use super::{Ideal, QI};
    use crate::discriminant;

    fn check(ideal: Ideal, s: &str) {
        let mut t = String::new();
        let mut fmt = core::fmt::Formatter::new(&mut t, core::fmt::FormattingOptions::new());
        ideal.latex(&mut fmt).unwrap();
        assert_eq!(s, t);
    }

    #[test]
    #[rustfmt::skip]
    fn latex_test_1() {
        let _guard = discriminant::DISC_TEST_LOCK.lock().unwrap();
        unsafe { discriminant::set(NonZeroI64::new(-6).unwrap(), false).unwrap() };

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
        unsafe { discriminant::set(NonZeroI64::new(-7).unwrap(), false).unwrap() };

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
