use num::BigUint;
use serde::Deserialize;
use std::{
    collections::{
        btree_map::Entry::{Occupied, Vacant},
        BTreeMap,
    },
    fs,
    io::{BufReader, Write},
    process::{Command, Stdio},
    time::{Duration, Instant},
};

struct DropGuard(String);

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

pub fn factor(n: &BigUint) -> anyhow::Result<Vec<(BigUint, u32)>> {
    let t1 = std::time::Instant::now();
    #[allow(clippy::transmute_undefined_repr)]
    let td = unsafe { std::mem::transmute::<Instant, Duration>(t1) };
    let result = format!("../data/factor/q{}.json", td.as_nanos());
    let mut child = Command::new("yafu")
        .arg("factor(@)")
        .arg("-factorjson")
        .arg(&result)
        .stdin(Stdio::piped())
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()?;
    let stdin = child
        .stdin
        .as_mut()
        .ok_or_else(|| anyhow::anyhow!("no stdin"))?;
    write!(stdin, "{n}")?;
    child.wait()?;

    let file = fs::File::open(&result)?;

    let _guard = DropGuard(result);
    let reader = BufReader::new(file);
    let output = serde_json::from_reader::<_, YafuOutput>(reader)?;

    let mut factors = BTreeMap::new();
    for factor in output.factors {
        let factor = factor.parse()?;
        match factors.entry(factor) {
            Occupied(mut occupied) => *occupied.get_mut() += 1,
            Vacant(vacant) => { vacant.insert(1); },
        }
    }

    Ok(factors.into_iter().collect())
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use num::BigUint;

    use super::factor;

    #[test]
    #[rustfmt::skip]
    fn test_factor() {
        let n = BigUint::from_str(&"4".repeat(67)).unwrap();
        let result = factor(&n).unwrap();
        assert_eq!(
            result,
            &[
                (BigUint::from(2u32), 2),
                (BigUint::from_str("493121").unwrap(), 1),
                (BigUint::from_str("79863595778924342083").unwrap(), 1),
                (BigUint::from_str("28213380943176667001263153660999177245677").unwrap(), 1),
            ]
        );
    }
}
