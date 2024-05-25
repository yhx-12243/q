use std::num::Wrapping;

use hashbrown::HashSet;
use num::{
    bigint::{IntDigits, Sign},
    BigInt, BigUint, Integer, One, Signed, Zero,
};

use crate::{qi::QI, qr::quadratic_residue};

#[inline]
fn mod_2_64_signed(x: &BigInt) -> Wrapping<u64> {
    let y = mod_2_64(x.magnitude());
    if x.is_negative() {
        -y
    } else {
        y
    }
}

#[inline]
fn mod_2_64(x: &BigUint) -> Wrapping<u64> {
    Wrapping(if let Some(r) = x.digits().first() {
        *r
    } else {
        0
    })
}

struct Mat2By2 {
    a: BigUint,
    b: BigUint,
    c: BigUint,
    d: BigUint,
}

/// inner method of divide & conquer.
fn convergent_inner(qs: &[BigUint]) -> Mat2By2 {
    if qs.len() == 1 {
        return Mat2By2 {
            a: qs[0].clone(),
            b: BigUint::one(),
            c: BigUint::one(),
            d: BigUint::ZERO,
        };
    }

    let (left, right) = unsafe { qs.split_at_unchecked(qs.len() / 2) };
    let left = convergent_inner(left);
    let right = convergent_inner(right);

    // clippy is so stupid to discover that it is a matrix multiplication.
    #[allow(clippy::suspicious_operation_groupings)]
    Mat2By2 {
        a: &left.a * &right.a + &left.b * &right.c,
        b: &left.a * &right.b + &left.b * &right.d,
        c: &left.c * &right.a + &left.d * &right.c,
        d: &left.c * &right.b + &left.d * &right.d,
    }
}

/// generate convergents of given continued fraction.
fn convergent(qs: &[BigUint]) -> (BigUint, BigUint) {
    let mat = convergent_inner(qs);
    (mat.a, mat.c)
}

/// solve P x^2 - Q x y + R y^2 where P, R = 0, 1 (R can't be zero).
#[allow(clippy::option_option)]
fn solve_simple(P: &BigUint, Q: &BigUint, R: &BigUint) -> Option<Option<(BigUint, BigUint)>> {
    if R.is_one() {
        Some(Some((BigUint::ZERO, BigUint::one())))
    } else if P.is_one() {
        Some(Some((BigUint::one(), BigUint::ZERO)))
    } else if P.is_zero() {
        let (q, r) = (R - 1u32).div_rem(Q);
        Some(r.is_zero().then(|| (q, BigUint::one())))
    } else {
        None
    }
}

/// Solve P x^2 - Q x y + R y^2 = 1 by convergents of x/y.
fn solve_by_convergents(
    P: &BigUint,
    Q: &BigUint,
    R: &BigUint,
    mut x: BigUint,
    mut y: BigUint,
) -> Option<(BigUint, BigUint)> {
    if let Some(res) = solve_simple(P, Q, R) {
        return res;
    }

    if Q.is_zero() {
        return None;
    }

    let (mut cx, mut cy, mut ccx, mut ccy) = (
        Wrapping(1u64),
        Wrapping(0u64),
        Wrapping(0u64),
        Wrapping(1u64),
    );

    let P_ = mod_2_64(P);
    let Q_ = mod_2_64(Q);
    let R_ = mod_2_64(R);

    let empirical = (x.len() + y.len() + 5) * 2;
    let mut qs = Vec::with_capacity(empirical);
    let mut indices = Vec::with_capacity(empirical);
    loop {
        let (q, r) = x.div_rem(&y);

        let q_ = mod_2_64(&q);
        ccx += cx * q_;
        ccy += cy * q_;

        qs.push(q);

        if (ccx * (ccx * P_ - ccy * Q_) + ccy * ccy * R_).0 == 1 {
            indices.push(qs.len());
        }

        if r.is_zero() {
            for index in indices {
                let (x, y) = convergent(unsafe { qs.get_unchecked(..index) });
                if &x * &x * P + &y * &y * R == &x * &y * Q + 1u32 {
                    return Some((x, y));
                }
            }
            return None;
        }

        core::mem::swap(&mut cx, &mut ccx);
        core::mem::swap(&mut cy, &mut ccy);
        x = y;
        y = r;
    }
}

/// A quadratic irrational (a + √d) / b, where b | a^2 - d.
#[derive(Clone, Hash, PartialEq, Eq)]
struct Qin {
    a: BigInt,
    b: BigInt,
}

// (try -d 1263816031 --input 2 !)
/// Solve P x^2 - Q x y + R y^2 = ±1 by convergents of x.
fn solve_by_convergents_QIN(
    P: &BigInt,
    Q: &BigUint,
    R: &BigUint,
    D: u64,
    mut x: Qin,
) -> Option<(/* bool, */ BigUint, BigUint)> {
    if R.is_one() {
        return Some((/* false, */ BigUint::ZERO, BigUint::one()));
    } else if P.magnitude().is_one() {
        return Some((/* P.is_negative(), */ BigUint::one(), BigUint::ZERO));
    } else if P.is_zero() {
        let (q, r) = (R - 1u32).div_rem(Q);
        if r.is_zero() {
            return Some((/* false, */ q, BigUint::one()));
        }
        let (q, r) = (R + 1u32).div_rem(Q);
        return r.is_zero().then(|| (/* true, */ q, BigUint::one()));
    }

    let (mut cx, mut cy, mut ccx, mut ccy) = (
        Wrapping(1u64),
        Wrapping(0u64),
        Wrapping(0u64),
        Wrapping(1u64),
    );

    let D_sqrt = D.isqrt();
    let P_ = mod_2_64_signed(P);
    let Q_ = mod_2_64(Q);
    let R_ = mod_2_64(R);

    let empirical = ((D.ilog2() + 5) * 2) as usize;
    let mut hash = HashSet::with_capacity(empirical);
    let mut qs = Vec::with_capacity(empirical);
    let mut indices = Vec::with_capacity(empirical);
    loop {
        let q = (&x.a + (D_sqrt + x.b.is_negative() as u64)).magnitude() / x.b.magnitude();

        let q_ = mod_2_64(&q);
        ccx += cx * q_;
        ccy += cy * q_;

        let trial = ccx * (ccx * P_ - ccy * Q_) + ccy * ccy * R_;
        if trial.0 == 1 || trial.0 == u64::MAX {
            indices.push(qs.len() + 1);
        }

        if !hash.insert(x.clone()) {
            qs.push(q);

            // println!("duplicate at {} tries, indices = {indices:?}", hash.len());
            for index in indices {
                let (x, y) = convergent(unsafe { qs.get_unchecked(..index) });
                let x = BigInt::from(x);
                let trial = &x * (&x * P - BigInt::from(&y * Q)) + BigInt::from(&y * &y * R);
                if trial.magnitude().is_one() {
                    return Some((/* trial.is_negative(), */ x.into_parts().1, y));
                }
            }
            return None;
        }

        core::mem::swap(&mut cx, &mut ccx);
        core::mem::swap(&mut cy, &mut ccy);

        let q_ = BigInt::from(q);
        x.a -= &x.b * &q_;
        x.b = (D - &x.a * &x.a) / &x.b;
        x.a = -x.a;

        qs.push(q_.into_parts().1);
    }
}

/// Solve x^2 + D y^2 = p.
pub fn work_neg(D: u64, p: &BigUint) -> Option<QI> {
    #[allow(clippy::verbose_bit_mask)]
    let e = !D & 3 == 0;

    let mut qr = {
        let mut d_p = BigUint::from(D) % p;
        if !d_p.is_zero() {
            d_p = p - d_p;
        }
        quadratic_residue(&d_p, p)?
    };

    if e {
        #[rustfmt::skip]
        if p.digits() == [2] {
            return (D == 7).then(|| QI { a: BigInt::one(), b: BigInt::one() });
        }

        if qr.is_even() {
            qr = p - qr;
        }

        let Q = qr;
        let R = p;
        let P = (&Q * &Q + D) / R / 4u32;

        // solve P y^2 - Q y z + R z^2 = 1
        let p_2 = &P * 2u32;
        let r_2 = R * 2u32;
        let (y, z) = if (&Q).max(&p_2) <= (&Q).max(&r_2) {
            solve_by_convergents(&P, &Q, R, Q.clone(), p_2)?
        } else {
            solve_by_convergents(&P, &Q, R, r_2, Q.clone())?
        };

        let x = {
            let a = &Q * &y;
            let b = p * &z * 2u32;
            if a < b {
                b - a
            } else {
                a - b
            }
        };

        Some(QI {
            a: x.into(),
            b: y.into(),
        })
    } else {
        let Q = qr;
        let R = p;
        let P = (&Q * &Q + D) / R;
        // P = (Q * Q - D) / R

        // solve P y^2 - 2 Q y z + R z^2 = 1
        let (y, z) = if (&Q).max(&P) <= (&Q).max(R) {
            solve_by_convergents(&P, &(&Q * 2u32), R, Q.clone(), P.clone())?
        } else {
            solve_by_convergents(&P, &(&Q * 2u32), R, R.clone(), Q.clone())?
        };

        let x = {
            let a = &Q * &y;
            let b = p * &z;
            if a < b {
                b - a
            } else {
                a - b
            }
        };

        Some(QI {
            a: x.into(),
            b: y.into(),
        })
    }
}

/// Solve x^2 - D y^2 = ±p.
pub fn work_pos(D: u64, p: &BigUint) -> Option<QI> {
    let e = D & 3 == 1;

    let qr = quadratic_residue(&(BigUint::from(D) % p), p)?;

    if e {
        if p.digits() == [2] {
            return if D > 16 && D & 7 == 1 {
                // solving x^2 + x y - (D-1)/4 y^2 = 2,
                // (let x = -2 z)
                // -(D-1)/8 y^2 - y z + 2 z^2 = 1
                let qin = Qin {
                    a: BigInt::from(-1),
                    b: BigInt::from(D / 4),
                };
                let (/* _is_neg, */ y, mut z) = solve_by_convergents_QIN(
                    &BigInt::from(-((D / 8) as i64)),
                    &BigUint::one(),
                    &BigUint::from(2u32),
                    D,
                    qin,
                )?;

                z *= 4u32;
                let x = if y < z { &z - &y } else { &y - &z };

                Some(QI {
                    a: x.into(),
                    b: y.into(),
                })
            } else {
                None
            };
        }

        let Q = if qr.is_even() { p - qr } else { qr };
        let R = p;
        let P = {
            let (sign, mut part) = (BigInt::from(&Q * &Q) - D).into_parts();
            part /= R;
            part /= 4u32;
            BigInt::from_biguint(sign, part)
        };

        // solve P y^2 - Q y z + R z^2 = 1
        let qin = Qin {
            a: if P.is_negative() {
                BigInt::from_biguint(Sign::Minus, Q.clone())
            } else {
                Q.clone().into()
            },
            b: P.abs() * 2,
        };

        let (/* _is_neg, */ y, z) = solve_by_convergents_QIN(&P, &Q, R, D, qin)?;

        let x = {
            let a = &Q * &y;
            let b = p * &z * 2u32;
            if a < b {
                b - a
            } else {
                a - b
            }
        };

        Some(QI {
            a: x.into(),
            b: y.into(),
        })
    } else {
        let Q = qr;
        let R = p;
        let P = {
            let (sign, mut part) = (BigInt::from(&Q * &Q) - D).into_parts();
            part /= R;
            BigInt::from_biguint(sign, part)
        };

        // solve P y^2 - 2 Q y z + R z^2 = ±1
        let qin = Qin {
            a: if P.is_negative() {
                BigInt::from_biguint(Sign::Minus, Q.clone())
            } else {
                Q.clone().into()
            },
            b: P.abs(),
        };

        let (/* _is_neg, */ y, z) = solve_by_convergents_QIN(&P, &(&Q * 2u32), R, D, qin)?;

        let x = {
            let a = &Q * &y;
            let b = p * &z;
            if a < b {
                b - a
            } else {
                a - b
            }
        };

        Some(QI {
            a: x.into(),
            b: y.into(),
        })
    }
}

/// Solve x^2 - D y^2 = ±p.
pub fn work(D: i64, p: &BigUint) -> Option<QI> {
    if D < 0 {
        work_neg(D.unsigned_abs(), p)
    } else {
        work_pos(D as u64, p)
    }
}
