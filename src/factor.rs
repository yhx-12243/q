use core::{error::Error, time::Duration};

use ciborium::Value as Cbor;
use num_bigint::{BigUint, IntDigits};
use num_traits::One;
use reqwest::blocking::Client;
use serde::Deserialize;

use crate::CONFIG;

mod factor_u64;

pub const DEFAULT_SERVER: &str = "https://43.138.56.99/factordb";

#[derive(Deserialize)]
struct FactorDBMirrorResp {
    #[serde(rename = "e")]
    exponent: u32,
    // #[serde(rename = "s")]
    // status: String,
    #[serde(rename = "E")]
    error: Option<String>,
    #[serde(rename = "f")]
    factors: Option<Cbor>,
}

/// factor the product of ns, take advantage of the known product.
#[allow(clippy::type_complexity)]
pub fn factor<const N: usize>(ns: [&BigUint; N]) -> Result<[Vec<(BigUint, u32)>; N], Box<dyn Error>> {
    #[inline]
    fn from_cbor(cbor: Cbor) -> Option<BigUint> {
        match cbor {
            Cbor::Integer(n) => Some(BigUint::from(u128::try_from(n).ok()?)),
            Cbor::Tag(2, box Cbor::Bytes(mut bytes)) => {
                bytes.reverse();
                Some(BigUint::from_bytes_le(&bytes))
            }
            _ => None,
        }
    }

    let mut result = [const { Vec::new() }; N];

    if ns.iter().all(|n| n.is_one()) {
        return Ok(result);
    }

    let url = CONFIG
        .get()
        .map_or_else(|| DEFAULT_SERVER, |config| &*config.factordb_mirror_server);

    let client = Client::builder().connect_timeout(const { Duration::from_secs(5) }).build()?;

    'outer: for (i, n) in ns.into_iter().enumerate() {
        if n.is_one() { continue; } // just leave empty.
        for (j, &m) in unsafe { ns.get_unchecked(..i) }.iter().enumerate() {
            if *m == *n {
                let [ri, rj] = unsafe { result.get_disjoint_unchecked_mut([i, j]) };
                ri.clone_from(rj);
                continue 'outer;
            }
        }

        if let [n] = *n.digits() {
            let factors = factor_u64::factor(n);
            *unsafe { result.get_unchecked_mut(i) } = factors
                .into_iter()
                .map(|(p, a)| (p.into(), a))
                .collect();
            continue;
        }

        let res = client.post(url).query(&[("factor", "")]).body(n.to_bytes_be()).send()?.bytes()?;
        let FactorDBMirrorResp {
            exponent,
            error,
            factors,
        } = ciborium::de::from_reader(&*res)?;

        if let Some(Cbor::Map(factors)) = factors
            && let Some(factors) = factors.into_iter().map(
                |(cbor, exp)| Some((from_cbor(cbor)?, u32::try_from(exp.as_integer()?).ok()? * exponent))
            ).collect() {
            // SAFETY: result and ns are both array of size N.
            *unsafe { result.get_unchecked_mut(i) } = factors;
            continue;
        }

        return Err(error.unwrap_or_else(|| format!("factorization of {n} failed")).into());
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use core::str::FromStr;

    use num_bigint::BigUint;

    use super::factor;

    #[test]
    #[rustfmt::skip]
    fn test_factor() {
        let result = factor([
            &BigUint::from(2_021_010_889u32),
            &BigUint::from(2_021_011_832u32),
            &BigUint::from(1u32),
            &BigUint::from_str(&"4".repeat(67)).unwrap(),
            &BigUint::from_str(&"4".repeat(67)).unwrap(),
            &BigUint::from_str(&"4".repeat(67)).unwrap(),
            &BigUint::from_str(&"4".repeat(67)).unwrap(),
        ]).unwrap();
        assert_eq!(
            result,
            [
                [
                    (BigUint::from(2_021_010_889u32), 1),
                ].as_slice(),
                [
                    (BigUint::from(2u32), 3),
                    (BigUint::from(7u32), 1),
                    (BigUint::from(1789u32), 1),
                    (BigUint::from(20173u32), 1),
                ].as_slice(),
                [
                ].as_slice(),
                [
                    (BigUint::from(2u32), 2),
                    (BigUint::from_str("493121").unwrap(), 1),
                    (BigUint::from_str("79863595778924342083").unwrap(), 1),
                    (BigUint::from_str("28213380943176667001263153660999177245677").unwrap(), 1),
                ].as_slice(),
                [
                    (BigUint::from(2u32), 2),
                    (BigUint::from_str("493121").unwrap(), 1),
                    (BigUint::from_str("79863595778924342083").unwrap(), 1),
                    (BigUint::from_str("28213380943176667001263153660999177245677").unwrap(), 1),
                ].as_slice(),
                [
                    (BigUint::from(2u32), 2),
                    (BigUint::from_str("493121").unwrap(), 1),
                    (BigUint::from_str("79863595778924342083").unwrap(), 1),
                    (BigUint::from_str("28213380943176667001263153660999177245677").unwrap(), 1),
                ].as_slice(),
                [
                    (BigUint::from(2u32), 2),
                    (BigUint::from_str("493121").unwrap(), 1),
                    (BigUint::from_str("79863595778924342083").unwrap(), 1),
                    (BigUint::from_str("28213380943176667001263153660999177245677").unwrap(), 1),
                ].as_slice(),
            ],
        );
    }
}
