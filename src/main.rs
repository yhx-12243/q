#![feature(
    box_patterns,
    const_trait_impl,
    core_io_internals,
    debug_closure_helpers,
    fmt_internals,
    integer_casts,
    likely_unlikely,
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
        default_value = factor::DEFAULT_SERVER,
        value_name = "url",
        help = "The URL of the factordb mirror server"
    )]
    factordb_mirror_server: String,
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

fn main() -> Result<(), Box<dyn core::error::Error>> {
    use clap::Parser;
    use ideal::Ideal;

    let args @ Args { plain_tex, .. } = Args::parse();
    unsafe { discriminant::set(args.D, args.plain_tex)? };

    CONFIG
        .set(args)
        .map_err(|_| Box::<dyn core::error::Error>::from("unable to set config"))?;

    let mut ideal = Ideal::read(std::io::stdin().lock())?;

    ideal.reduce();

    {
        use core::fmt::from_fn;
        use std::io::Write;

        let fmt = if plain_tex { Ideal::tex } else { Ideal::latex };

        let mut stdout = std::io::stdout().lock();
        write!(stdout, "{}=", from_fn(|f| fmt(&ideal, f)))?;

        let ideals = ideal.factor()?;

        if ideals.is_empty() {
            stdout.write_all(b"\\left(1\\right)")?;
        }

        for (ideal, exp) in ideals {
            write!(stdout, "{}", from_fn(|f| fmt(&ideal, f)))?;
            if exp > 1 {
                write!(stdout, "^{{{exp}}}")?;
            }
        }
    }

    Ok(())
}
