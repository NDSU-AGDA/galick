# Migration to Galick

The R package and GitHub repository are now `galick`. Use
`DidierMurilloF/galick` for installation and repository links.

## Code and saved results

- Load `library(galick)` and use `galick::` for qualified calls.
- Public function names and argument names are unchanged.
- Result classes now use `galick_`: `galick_result`, `galick_run`,
  `galick_simulation`, `galick_evaluation` and `galick_panel_design`.
- Raised conditions now use `galick_error_`. Collected issue codes are unchanged.
- Update `system.file()` and `vignette()` calls to `package = "galick"`.
- Recreate saved pre-rename S3 objects rather than editing their classes manually.

The rename is a pre-release namespace break; no compatibility aliases are installed.

## Tools and records

The Beagle cache is now `tools::R_user_dir("galick", "cache")`. To reuse a JAR
stored elsewhere, set `BEAGLE_PATH` or pass `beagle_jar` explicitly. The
`BEAGLE_PATH` variable describes the external tool and has not changed.

The internal run-directory option is `galick.run_dir`. The opt-in network test
switch is `GALICK_TEST_NETWORK=true`.

New Galick run records identify the package under `versions$galick` and use the
version 2.0.0 run-record schema. The schema bump records the changed required
package key. Simulation, evaluation and panel-design schemas remain at 1.0.0.
Historical records remain historical artifacts; do not
rewrite them. Old chromosome resume records are not reusable as Galick jobs:
package versions and fingerprints are part of the computation identity.

## Repository administration

The GitHub repository has been renamed to `DidierMurilloF/galick`. GitHub Pages
configuration now targets `https://didiermurillof.github.io/galick/`; deployment
is not automatic. CI builds a downloadable website preview; hosting is configured separately.
