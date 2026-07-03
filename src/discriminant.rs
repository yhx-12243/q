use core::{error::Error, num::NonZeroI64};

static mut DISCRIMINANT: NonZeroI64 = NonZeroI64::new(-1).unwrap();
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

pub unsafe fn set(d: NonZeroI64, plain: bool) -> Result<(), Box<dyn Error>> {
    // this is a critical section, but we ensure only called once or using lock.
    unsafe {
        if let Some(e) = d.get().checked_isqrt() && e * e == d.get() {
            return Err("discriminant can't be a square".into());
        }
        DISCRIMINANT = d;
        #[allow(clippy::deref_addrof)]
        if d.get() == -1 {
            "i".clone_into(&mut *&raw mut DISC_STR);
            if plain { "{\\rm i}" } else { "\\mathrm i" }.clone_into(&mut *&raw mut DISC_LATEX);
        } else {
            DISC_STR = format!("√{d}");
            DISC_LATEX = if matches!(d.get(), 0..10) {
                format!("\\sqrt{d}")
            } else {
                format!("\\sqrt{{{d}}}")
            };
        }
        Ok(())
    }
}
