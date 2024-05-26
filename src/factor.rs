use num::{BigUint, One};
use serde::Deserialize;
use std::{
    collections::{
        btree_map::Entry::{Occupied, Vacant},
        BTreeMap,
    },
    fs,
    io::{BufRead, BufReader, Write},
    path::PathBuf,
    process::{Command, Stdio},
    time::{Duration, Instant},
};

use crate::YAFU_DIR;

struct DropGuard(PathBuf);

impl Drop for DropGuard {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.0);
    }
}

#[derive(Deserialize)]
struct YafuOutput {
    #[serde(rename = "factors-prime")]
    factors: Vec<String>,
}

/// factor the product of ns, take advantage of the known product.
pub fn factor<const N: usize>(ns: [&BigUint; N]) -> anyhow::Result<[Vec<(BigUint, u32)>; N]> {
    let mut result = [const { Vec::new() }; N];

    if ns.iter().all(|n| n.is_one()) {
        return Ok(result);
    }

    let t1 = std::time::Instant::now();
    #[allow(clippy::transmute_undefined_repr)]
    let td = unsafe { std::mem::transmute::<Instant, Duration>(t1) };
    let mut path = YAFU_DIR
        .get()
        .map_or_else(|| PathBuf::from("../data/factor"), PathBuf::clone);
    path.push(format!("q{}.json", td.as_nanos()));
    let mut child = Command::new("yafu")
        .arg("factor(@)")
        .arg("-factorjson")
        .arg(&path)
        .stdin(Stdio::piped())
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()?;
    let stdin = child
        .stdin
        .as_mut()
        .ok_or_else(|| anyhow::anyhow!("no stdin"))?;
    for n in ns {
        if !n.is_one() {
            writeln!(stdin, "{n}")?;
        }
    }
    child.wait()?;

    let file = fs::File::open(&path)?;

    let _guard = DropGuard(path);
    let mut reader = BufReader::new(file);

    let mut buf = String::new();
    for (i, n) in ns.into_iter().enumerate() {
        if n.is_one() { continue; } // just leave empty.
        buf.clear();
        reader.read_line(&mut buf)?;
        let YafuOutput { factors } = serde_json::from_str(&buf)?;
        let mut map = BTreeMap::new();
        for factor in factors {
            let factor = factor.parse()?;
            match map.entry(factor) {
                Occupied(mut occupied) => *occupied.get_mut() += 1,
                Vacant(vacant) => { vacant.insert(1); },
            }
        }
        // SAFETY: result and ns are both array of size N.
        unsafe {
            *result.get_unchecked_mut(i) = map.into_iter().collect();
        }
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use num::BigUint;

    use super::factor;

    #[test]
    #[rustfmt::skip]
    fn test_factor() {
        let result = factor([
            &BigUint::from(2021010889u32),
            &BigUint::from(2021011832u32),
            &BigUint::from(1u32),
            &BigUint::from_str(&"4".repeat(67)).unwrap(),
        ]).unwrap();
        assert_eq!(
            result,
            [
                [
                    (BigUint::from(2021010889u32), 1),
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
            ],
        );
    }
}
