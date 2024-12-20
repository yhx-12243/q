#![feature(
    debug_closure_helpers,
    exit_status_error,
    fmt_internals,
    formatting_options,
    get_many_mut,
    integer_sign_cast,
    let_chains,
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
    #[arg(long, help = "Whether to output Plain TeX code instead of LaTeX")]
    plain_tex: bool,
}

static CONFIG: std::sync::OnceLock<Args> = std::sync::OnceLock::new();

fn main() -> anyhow::Result<()> {
    use clap::Parser;
    use ideal::Ideal;

    let args @ Args { plain_tex, .. } = Args::parse();
    unsafe { discriminant::set(args.D, args.plain_tex)? };

    std::fs::create_dir_all(&args.dir)?;
    CONFIG
        .set(args)
        .map_err(|_| anyhow::anyhow!("unable to set config"))?;

    let mut ideal = Ideal::read(std::io::stdin().lock())?;

    ideal.reduce();

    {
        use core::fmt::{Arguments, rt::Argument};
        use std::io::Write;

        let fmt = if plain_tex { Ideal::tex } else { Ideal::latex };

        let mut stdout = std::io::stdout().lock();
        stdout.write_fmt(Arguments::new_v1(&["", "="], &[Argument::new(&ideal, fmt)]))?;

        let ideals = ideal.factor()?;

        if ideals.is_empty() {
            stdout.write_all(b"\\left(1\\right)")?;
        }

        for (ideal, exp) in ideals {
            stdout.write_fmt(Arguments::new_v1(&[""], &[Argument::new(&ideal, fmt)]))?;
            if exp > 1 {
                write!(stdout, "^{{{exp}}}")?;
            }
        }
    }

    Ok(())
}
