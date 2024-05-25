#![feature(
    debug_closure_helpers,
    fmt_helpers_for_derive,
    fmt_internals,
    get_many_mut,
    isqrt,
    let_chains,
    raw_ref_op,
    slice_ptr_get,
    stmt_expr_attributes,
)]

mod discriminant;
mod factor;
mod ideal;
mod pell;
mod qi;
mod qr;

use core::num::NonZeroI64;

use ideal::Ideal;

#[derive(clap::Parser)]
struct Args {
    #[arg(short = 'D', value_name = "discriminant")]
    D: NonZeroI64,
    #[arg(short, long, value_name = "input")]
    input: String,
}

fn main() -> anyhow::Result<()> {
    use clap::Parser;

    let args = Args::parse();
    unsafe { discriminant::set(args.D)? };

    let mut ideal = args.input.parse::<Ideal>()?;

    let ideals = ideal.factor()?;

    {
        use core::fmt::{rt::Argument, Arguments, Formatter};
        use std::io::Write;

        #[allow(clippy::unit_arg, clippy::unnecessary_wraps)]
        fn latex_wrapper(ideal: &Ideal, fmt: &mut Formatter) -> core::fmt::Result {
            Ok(ideal.latex(fmt))
        }

        let mut stdout = std::io::stdout();
        stdout.write_fmt(Arguments::new_v1(
            &[""],
            &[Argument::new(&ideal, latex_wrapper)],
        ))?;
        stdout.write_all(b"=")?;
        for (ideal, exp) in ideals {
            stdout.write_fmt(Arguments::new_v1(
                &[""],
                &[Argument::new(&ideal, latex_wrapper)],
            ))?;
            if exp > 1 {
                write!(stdout, "^{exp}")?;
            }
        }
    }

    Ok(())
}
