fn main() {
    println!("cargo::rustc-link-search=native=/opt/homebrew/Cellar/gmp/6.3.0/lib");
    println!("cargo::rustc-link-lib=gmp");
}
