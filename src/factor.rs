use core::time::Duration;

use bytes::Bytes;
use ciborium::Value as Cbor;
use futures_util::{
    FutureExt,
    future::{err, join_all},
};
use num_bigint::BigUint;
use num_traits::One;
use reqwest::{Client, Response};
use serde::Deserialize;

use crate::CONFIG;

#[derive(Deserialize)]
struct FactorDBMirrorResp {
    #[serde(rename = "e")]
    exponent: u32,
    // #[serde(rename = "s")]
    // status: String,
    #[serde(rename = "E")]
    error: Option<String>,
    #[serde(rename = "f")]
    factors: Option<Vec<(Cbor, u32)>>,
}

/// factor the product of ns, take advantage of the known product.
pub async fn factor<const N: usize>(ns: [&BigUint; N]) -> anyhow::Result<[Vec<(BigUint, u32)>; N]> {
    #[inline]
    fn extract_body(res: reqwest::Result<Response>) -> impl Future<Output = reqwest::Result<Bytes>> {
        match res {
            Ok(res) => res.bytes().left_future(),
            Err(e) => err(e).right_future(),
        }
    }

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
    let mut dup = [usize::MAX; N];

    if ns.iter().all(|n| n.is_one()) {
        return Ok(result);
    }

    let url = CONFIG.get().map_or_else(
        || "https://3.145.200.48/factordb",
        |config| &*config.factordb_mirror_server,
    );

    let client = Client::builder().connect_timeout(const { Duration::from_secs(5) }).build()?;

    let futs = ns.into_iter().enumerate().filter_map(|(i, n)| {
        if n.is_one() { return None; }
        for (j, &m) in unsafe { ns.get_unchecked(..i).iter().enumerate() } {
            if m == n {
                dup[i] = j;
                return None;
            }
        }
        Some(client.post(url).query(&[("factor", "")]).body(n.to_bytes_be()).send().then(extract_body))
    });
    let mut ress = join_all(futs).await.into_iter();

    for (i, n) in ns.into_iter().enumerate() {
        if n.is_one() { continue; } // just leave empty.
        unsafe {
            let j = *dup.get_unchecked(i);
            if j != usize::MAX {
                let [ri, rj] = result.get_disjoint_unchecked_mut([i, j]);
                ri.clone_from(rj);
                continue;
            }
        }
        let res = unsafe { ress.next().unwrap_unchecked() }?;
        let FactorDBMirrorResp {
            exponent,
            error,
            factors,
        } = ciborium::de::from_reader(&*res)?;

        let factors = if let Some(factors) = factors {
            factors.into_iter().map(|(cbor, exp)| Some((from_cbor(cbor)?, exp * exponent))).collect::<Option<_>>()
        } else {
            return Err(anyhow::Error::msg(error.unwrap_or_else(|| "factorization failed".to_string())));
        };
        // SAFETY: result and ns are both array of size N.
        unsafe {
            *result.get_unchecked_mut(i) = factors.ok_or_else(||anyhow::anyhow!("failed to construct factorization expression"))?;
        }
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use core::str::FromStr;

    use futures_executor::block_on;
    use num_bigint::BigUint;

    use super::factor;

    #[test]
    #[rustfmt::skip]
    fn test_factor() {
        let result = block_on(factor([
            &BigUint::from(2_021_010_889u32),
            &BigUint::from(2_021_011_832u32),
            &BigUint::from(1u32),
            &BigUint::from_str(&"4".repeat(67)).unwrap(),
            &BigUint::from_str(&"4".repeat(67)).unwrap(),
            &BigUint::from_str(&"4".repeat(67)).unwrap(),
            &BigUint::from_str(&"4".repeat(67)).unwrap(),
        ])).unwrap();
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
