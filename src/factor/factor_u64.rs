use core::hint::unlikely;
use std::collections::{
    BTreeMap,
    btree_map::Entry::{Occupied, Vacant},
};

use num_integer::Integer;
use rand::{
    Rng, SeedableRng, TryRng,
    rngs::{SmallRng, SysRng},
};

#[allow(clippy::large_const_arrays, clippy::many_single_char_names)]
mod small_primes {
    const __GENERATE_SMALL_PRIMES: ([u16; 6542], [u16; 65536]) = {
        let mut p = [0; 6542];
        let mut c = [u16::MAX; 0x10000];
        let mut n = 0;
        let mut i = 2u16;
        loop {
            if c[i as usize] == u16::MAX {
                p[n] = i;
                c[i as usize] = n as u16;
                n += 1;
            }
            let mut j = 0usize;
            while j <= c[i as usize] as usize {
                let Some(v) = i.checked_mul(p[j]) else {
                    break;
                };
                c[v as usize] = j as u16;
                j += 1;
            }
            if let Some(nxt) = i.checked_add(1) {
                i = nxt;
            } else {
                break;
            }
        }
        assert!(n == 6542, "sieve prime error");
        (p, c)
    };
    pub const SMALL_PRIMES: [u16; 6542] = __GENERATE_SMALL_PRIMES.0;
    pub const LEAST_FACTORS: [u16; 65536] = __GENERATE_SMALL_PRIMES.1;
}

use small_primes::{LEAST_FACTORS, SMALL_PRIMES};

const SIEVE_MAX: u64 = 1 << 16;

#[cfg(target_arch = "x86_64")]
#[inline(always)]
fn rem(h: u64, l: u64, b: u64) -> u64 {
    let r: u64;
    unsafe {
        core::arch::asm!(
            "div {b}",
            inlateout("rdx") h => r,
            in("rax") l,
            b = in(reg) b,
            options(pure, nomem, nostack)
        );
        r
    }
}

#[inline]
const fn sub_mod(a: u64, b: u64, m: u64) -> u64 {
    if a < b {
        a.wrapping_sub(b).wrapping_add(m)
    } else {
        a - b
    }
}

#[inline]
fn mul_mod(a: u64, b: u64, m: u64) -> u64 {
    let p = u128::from(a) * u128::from(b);
    #[cfg(target_arch = "x86_64")]
    {
        rem((p >> 64) as u64, p as u64, m)
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        (p % u128::from(m)) as u64
    }
}

#[inline]
fn sqr_add_mod(a: u64, b: u64, m: u64) -> u64 {
    let p = u128::from(a) * u128::from(a) + u128::from(b);
    #[cfg(target_arch = "x86_64")]
    {
        rem((p >> 64) as u64, p as u64, m)
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        (p % u128::from(m)) as u64
    }
}

pub(super) fn power_mod(mut a: u64, mut n: u64, m: u64, mut c: u64) -> u64 {
    while n != 0 {
        if (n & 1) != 0 {
            c = mul_mod(c, a, m);
        }
        n >>= 1;
        a = mul_mod(a, a, m);
    }
    c
}

#[allow(clippy::manual_is_multiple_of, clippy::many_single_char_names)]
pub(super) fn miller_rabin(n: u64) -> bool {
    const TEST: [u64; 8] = [2, 3, 5, 7, 11, 13, 359, 911];

    assert!(n >= SIEVE_MAX);

    let c = (n - 1).trailing_zeros();
    let s = (n - 1) >> c;
    for t in TEST {
        if n % t == 0 {
            return false;
        }
        let mut u = power_mod(t, s, n, 1);
        for _ in 0..c {
            let v = mul_mod(u, u, n);
            if u != 1 && u != n - 1 && v == 1 {
                return false;
            }
            u = v;
        }
        if u != 1 {
            return false;
        }
    }
    true
}

#[allow(clippy::many_single_char_names, clippy::verbose_bit_mask)]
fn pollard_rho(n: u64, c: u64, mut x: u64) -> u64 {
    let mut y = x;
    let mut p = 1;
    let mut i = 1;
    let mut k = 2;
    while k <= 0x20000 {
        x = sqr_add_mod(x, c, n);
        let d = sub_mod(y, x, n);
        let q = mul_mod(p, d, n);
        if unlikely(q == 0) {
            let x = p.gcd(&n);
            if x != 1 {
                return x;
            }
        } else {
            p = q;
        }
        if (i & 127) == 0 {
            let x = p.gcd(&n);
            if x != 1 {
                return x;
            }
        }
        i += 1;
        if i == k {
            y = x;
            k <<= 1;
        }
    }
    0
}

#[inline]
fn push(p: u64, a: u32, result: &mut BTreeMap<u64, u32>) {
    match result.entry(p) {
        Occupied(mut occupied) => *occupied.get_mut() += a,
        Vacant(vacant) => { vacant.insert(a); },
    }
}

#[allow(clippy::manual_is_multiple_of)]
pub fn factor_inner<T: Rng + ?Sized>(mut n: u64, result: &mut BTreeMap<u64, u32>, rng: &mut T) {
    while n != 1 {
        if n >= SIEVE_MAX {
            if miller_rabin(n) {
                return push(n, 1, result);
            }
            let mut c = 97;
            let mut p;
            // let r = unsafe { UniformInt::<u64>::new(0, n).unwrap_unchecked() };
            loop {
                // let x = r.sample(rng);
                let x = ((u128::from(rng.next_u64()) * u128::from(n)) >> 64) as u64;
                p = pollard_rho(n, c, x);
                if p != 0 {
                    break;
                }
                c += 1;
            }
            factor_inner(p, result, rng);
            n /= p;
        } else {
            let pid = LEAST_FACTORS[n as usize];
            let p = SMALL_PRIMES[pid as usize].into();
            let mut a = 0;
            while n % p == 0 {
                n /= p;
                a += 1;
            }
            push(p, a, result);
        }
    }
}

pub fn factor(n: u64) -> BTreeMap<u64, u32> {
    let mut result = BTreeMap::new();
    let mut rng = SmallRng::seed_from_u64(SysRng.try_next_u64().unwrap_or_default());
    factor_inner(n, &mut result, &mut rng);
    result
}

#[cfg(test)]
mod test {
    use super::{miller_rabin, mul_mod, pollard_rho};

    #[test]
    fn mul_mod_test() {
        let a = 0x0123_4567_dead_beef;
        let b = 0xdead_beef_0123_4567;
        let m = 0xeaeb_eced_eeef_e0e1;
        let r = mul_mod(a, b, m);
        assert_eq!(r, 0x2cf4_0dfe_13ac_0fac);
    }

    #[test]
    fn miller_rabin_test() {
        let p = 1_099_999_999_999_999_999;
        assert!(miller_rabin(p));
    }

    #[test]
    fn pollard_rho_test() {
        let n = 31_372_013_730_798_557;
        let d = pollard_rho(n, 97, 12_744_588_943_418_861);
        assert!(d != 0 && n % d == 0);
    }
}
