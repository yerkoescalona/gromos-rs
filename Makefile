# gromos-rs — local CI equivalent
# Usage: make check   (fast, before commit)
#        make ci      (full suite, matches GitHub CI)

.PHONY: check ci fmt fmt-check lint-rust lint test

# Run before every commit (matches .claude/overview.md rule)
check:
	cargo fmt
	cargo check --tests

# Full CI equivalent — run before pushing
ci: fmt-check lint-rust test

fmt:
	cargo fmt

fmt-check:
	cargo fmt --check

lint-rust:
	cargo clippy --workspace --all-targets -- \
		-W warnings \
		-A clippy::upper_case_acronyms \
		-A non_snake_case

test:
	cargo test --workspace
