---
title: Install Galick and configure Beagle in R
description: Install the public Galick R package, configure Java and Beagle, and locate optional genetics PLINK, bgzip, tabix and bcftools tools.
---
# Install Galick

Galick requires **R 4.1 or later**. The R package has no compiled code of its own;
some dependencies and external tools contain compiled components.

## Install from CRAN

**Galick is not yet published on CRAN.** The following command is for use
after the first CRAN release; it does not install Galick today. Use GitHub below
for the current version.

```r
# Available only after Galick is published on CRAN.
install.packages("galick")
```

## Install from GitHub

This method is available now and installs the current default branch from the
public `NDSU-AGDA/galick` repository, together with its R dependencies.

```r
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}
pak::pak("NDSU-AGDA/galick")
```

## Load and check the installation

After installation, start a new R session if an older Galick version was loaded.

```r
library(galick)
packageVersion("galick")
```

Dependencies are `cli`, `data.table`, `digest`, `jsonlite`, `processx`, `tibble`,
`vroom` and `withr`. Validators, plain-text conversion, simulation and genotype
extraction do not require Java. Do not install external tools unless your workflow
needs them.

## Choose external tools for your workflow

| Workflow | Additional software |
|---|---|
| Validate, repair, convert plain VCFs, simulate or extract genotypes | None |
| Run Beagle, benchmark or evaluate imputation, design a marker panel | Java 8+ and a Beagle 5 JAR |
| Write compressed `.vcf.gz` exports or merged VCFs | `bgzip` |
| Index a Beagle output | `tabix` or `bcftools`; optional for a successful run |
| Convert binary BED/BIM/FAM files | A genetics PLINK 1.9-compatible executable |
| Run chromosomes through a future plan | The R package `future.apply` |

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

## Other tools

Reading `.vcf.gz` does not require bgzip. The PLINK bridge uses legacy
`--make-bed` / `--recode vcf` commands; use a compatible executable rather than
assuming every PLINK 2 build accepts the same options.

On Windows, PuTTY also installs an unrelated SSH client named `plink.exe`.
Galick verifies the genetics PLINK version banner and rejects that namesake.
It checks `plink1.9` before `plink` on `PATH`; alternatively pass the full path
to the genetics executable through `plink =` in the bridge functions.

Check which executables R can find:

```r
Sys.which(c("java", "bgzip", "tabix", "bcftools", "plink"))
```

An empty value means that name is not on R's `PATH`; it does not rule out a tool
at an explicit path. A nonempty `plink` result alone does not establish that it
is the genetics tool. Once your workflow's requirements are available, continue
with the [Quick tour](quick-tour.md).
