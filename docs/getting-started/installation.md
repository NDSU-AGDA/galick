# Installation

Galick requires **R 4.1 or later**. The R package has no compiled code of its own;
some dependencies and external tools contain compiled components.

## R package

```r
install.packages("pak")
pak::pak("DidierMurilloF/galick")
library(galick)
```

The installed namespace and repository are both named `galick`.
Until the package rename reaches the default branch, install
`DidierMurilloF/galick@feat/evaluate-imputation` instead. See
[Migration](migration.md) for changed classes, caches and records.

Dependencies are `cli`, `data.table`, `digest`, `jsonlite`, `processx`, `tibble`,
`vroom` and `withr`. Validators, plain-text conversion, simulation and genotype
extraction do not require Java. Do not install external tools unless your workflow
needs them.

## Java

Beagle execution requires Java 8 or later. Check the executable visible to R:

```r
Sys.which("java")
```

Install a JRE/JDK through your system package manager or
[Adoptium](https://adoptium.net/). Set `JAVA_HOME` if Java is not on `PATH`, or
pass an executable path through `run_beagle(java = ...)`. Galick checks the Java
version before launching a job.

## Beagle

The Beagle JAR is not bundled. Galick resolves it in this order:

1. The explicit `beagle_jar` argument.
2. `BEAGLE_PATH` in the R process environment.
3. `file.path(tools::R_user_dir("galick", "cache"), "beagle.jar")`.

With Java installed:

```r
jar <- install_beagle()
beagle_version(beagle_jar = jar)
```

`install_beagle()` supports the pinned `beagle.29Oct24.c8e.jar` build of version
5.4. It downloads over HTTPS and checks the version by running the JAR.
The installer does not verify a published cryptographic checksum.
It does not select the latest upstream release. For another build, use
`BEAGLE_PATH` or an explicit JAR path. Consult the
[official Beagle site](https://faculty.washington.edu/browning/beagle/beagle.html)
for its documentation and license.

To configure an existing JAR for the current R session:

```r
Sys.setenv(BEAGLE_PATH = "/absolute/path/to/beagle.jar")
beagle_version()
```

## Optional tools

| Tool | Required for |
|---|---|
| `bgzip` | Writing compressed `.vcf.gz` exports and merged outputs |
| `tabix` or `bcftools` | Optional post-Beagle indexing; indexing is otherwise recorded as skipped |
| PLINK 1.9-compatible executable | `vcf_to_plink()` and `plink_to_vcf()` BED/BIM/FAM bridge commands |
| `future.apply` | Chromosome execution with `parallel = "future"` |

Reading `.vcf.gz` does not require bgzip. The PLINK bridge uses legacy
`--make-bed` / `--recode vcf` commands; use a compatible executable rather than
assuming every PLINK 2 build accepts the same options.

Check which executables R can find:

```r
Sys.which(c("java", "bgzip", "tabix", "bcftools", "plink"))
```

An empty value means the executable is not on R's `PATH`. Once the required
tools are available, continue with the [Quick tour](quick-tour.md).
