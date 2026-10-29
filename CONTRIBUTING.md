# Contributing

## Toolchain

Use Rust 1.89 or newer.

## Building

```bash
cargo build --release
./target/release/vartovcf --help
```

## Mandatory Pre-submission Checks

All must pass before a change is submitted, and CI runs the same three:

```bash
cargo ci-fmt
cargo ci-lint
cargo ci-test
```
