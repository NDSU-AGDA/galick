# Installation

beagler is an R package that shells out to **Java**, **Beagle**, and
optionally **PLINK** / **bcftools** / **tabix**. The R package itself
is pure R; the external tools handle the heavy lifting.

## R package

=== "pak (recommended)"

    ```r
    install.packages("pak")
    pak::pak("DidierMurilloF/beagleR")
    ```

=== "remotes"

    ```r
    install.packages("remotes")
    remotes::install_github("DidierMurilloF/beagleR")
    ```

=== "From source"

    ```bash
    git clone https://github.com/DidierMurilloF/beagleR.git
    R CMD INSTALL beagleR
    ```

The package targets R ≥ 4.0. It pulls in `cli`, `jsonlite`, `processx`,
`tibble`, `vroom`, and `withr` — all on CRAN.

## Java

Required for `run_beagle()`. v0.1.0 needs **Java 8 or later**.

=== "macOS"

    ```bash
    brew install openjdk
    ```

=== "Ubuntu / Debian"

    ```bash
    sudo apt install default-jre
    ```

=== "Windows"

    Install [Adoptium Temurin](https://adoptium.net/) and ensure
    `java.exe` is on `PATH` (or set `JAVA_HOME`).

Verify from R:

```r
beagler:::find_java()
beagler:::detect_java_version(beagler:::find_java())
```

## Beagle

The Beagle JAR cannot ship with the R package (CRAN policy). beagler
discovers it in this order:

1. The `beagle_jar = ` argument to `run_beagle()`.
2. The `BEAGLE_PATH` environment variable.
3. `tools::R_user_dir("beagler", "cache")/beagle.jar`.

The easiest setup:

```r
library(beagler)
install_beagle()              # downloads pinned tested 5.4 to (3)
beagle_version()              # "5.4"
```

!!! note "What `install_beagle()` actually downloads"

    The pinned URL points at `beagle.29Oct24.c8e.jar` from
    [Brian Browning's official Beagle page](https://faculty.washington.edu/browning/beagle/beagle.html).
    We treat this as a *pinned tested* 5.4 jar — upstream may list a
    different build by the time you read this. If you need a specific
    build, set `BEAGLE_PATH` or pass `beagle_jar = ` directly.

If you'd rather manage Beagle outside R:

```bash
# Download once
curl -L https://faculty.washington.edu/browning/beagle/beagle.29Oct24.c8e.jar \
  -o /opt/beagle/beagle.jar

# Tell R where it lives
echo 'BEAGLE_PATH=/opt/beagle/beagle.jar' >> ~/.Renviron
```

## Optional tools

These are not required, but `run_beagle()` and the chromosome workflows
take advantage of them when present:

- **bcftools** — used to index Beagle's output (`run_beagle()` falls
  back to `tabix` or marks the index step as `skipped`).
- **tabix** — alternative to bcftools for indexing.
- **plink (≥ 1.9)** — required only if you use `vcf_to_plink()` /
  `plink_to_vcf()`.

=== "macOS"

    ```bash
    brew install bcftools htslib plink
    ```

=== "Ubuntu / Debian"

    ```bash
    sudo apt install bcftools tabix plink1.9
    ```

## Verify everything

```r
library(beagler)

# Java + Beagle
beagle_version()

# Optional tools (each returns a path or NULL)
beagler:::find_bcftools()
beagler:::find_tabix()
beagler:::find_plink()
```

Once `beagle_version()` returns a string, you're ready for the
[Quick tour](quick-tour.md).
