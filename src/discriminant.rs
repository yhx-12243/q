use core::num::NonZeroI64;

static mut DISCRIMINANT: NonZeroI64 = unsafe { NonZeroI64::new_unchecked(-1) };
static mut DISC_STR: String = String::new();
static mut DISC_LATEX: String = String::new();
#[cfg(test)]
pub static DISC_TEST_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());

pub fn get() -> NonZeroI64 {
    unsafe { DISCRIMINANT }
}

pub fn get_str() -> &'static str {
    #[allow(clippy::explicit_auto_deref)] // false positive
    unsafe { &*DISC_STR }
}

pub fn get_latex() -> &'static str {
    #[allow(clippy::explicit_auto_deref)] // false positive
    unsafe { &*DISC_LATEX }
}

pub fn is4kp1() -> bool {
    unsafe { DISCRIMINANT.get() & 3 == 1 }
}

pub unsafe fn set(d: NonZeroI64, plain: bool) -> anyhow::Result<()> {
    // this is a critical section, but we ensure only called once or using lock.
    unsafe {
        if let Some(e) = d.get().checked_isqrt() && e * e == d.get() {
            anyhow::bail!("discriminant can't be a square");
        }
        DISCRIMINANT = d;
        #[allow(clippy::deref_addrof)]
        if d.get() == -1 {
            "i".clone_into(&mut *&raw mut DISC_STR);
            if plain { "{\\rm i}" } else { "\\mathrm i" }.clone_into(&mut *&raw mut DISC_LATEX);
        } else {
            DISC_STR = format!("√{d}");
            DISC_LATEX = format!("\\sqrt{{{d}}}");
        };
        Ok(())
    }
}
