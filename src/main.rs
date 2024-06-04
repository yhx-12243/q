#![feature(
    array_try_map,
    debug_closure_helpers,
    exact_size_is_empty,
    exit_status_error,
    fmt_helpers_for_derive,
    fmt_internals,
    get_many_mut,
    hint_assert_unchecked,
    io_error_more,
    isqrt,
    iter_next_chunk,
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

#[derive(clap::Parser)]
#[command(
    version,
    about = "ideal unique factorization of quadratic integer domains"
)]
struct Args {
    #[arg(
        short = 'D',
        visible_alias = "discriminant",
        value_name = "discriminant",
        help = "The discriminant of the quadratic integer domain"
    )]
    D: core::num::NonZeroI64,
    #[arg(
        long,
        default_value = "./",
        value_name = "dir",
        help = "The directory of YAFU output"
    )]
    dir: std::path::PathBuf,
    #[arg(
        short = 'M',
        long,
        default_value_t = 10240,
        value_name = "number",
        help = "Max number of continuous fraction iteration to solve Pell equation/check for principal ideals"
    )]
    max_iter: usize,
}

static CONFIG: std::sync::OnceLock<Args> = std::sync::OnceLock::new();

fn main() -> anyhow::Result<()> {
    use clap::Parser;
    use ideal::Ideal;

    let args = Args::parse();
    unsafe { discriminant::set(args.D)? };

    std::fs::create_dir_all(&args.dir)?;
    CONFIG.set(args).map_err(|_| anyhow::anyhow!("unable to set config"))?;

    let mut ideal = Ideal::read(std::io::stdin().lock())?;

    ideal.reduce();

    {
        use core::fmt::{rt::Argument, Arguments};
        use std::io::Write;

        let mut stdout = std::io::stdout().lock();
        stdout.write_fmt(Arguments::new_v1(
            &["", "="],
            &[Argument::new(&ideal, Ideal::latex)],
        ))?;

        let ideals = ideal.factor()?;

        if ideals.is_empty() {
            stdout.write_all(b"\\left(1\\right)")?;
        }

        for (ideal, exp) in ideals {
            stdout.write_fmt(Arguments::new_v1(
                &[""],
                &[Argument::new(&ideal, Ideal::latex)],
            ))?;
            if exp > 1 {
                write!(stdout, "^{{{exp}}}")?;
            }
        }
    }

    Ok(())
}
