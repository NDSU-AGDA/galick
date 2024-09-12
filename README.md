# Install

```
devtools::install()
```


# Run beagleR

Make sure you have the environment variable "BEAGLE_PATH"

```
# .Renviron
BEAGLE_PATH=/path/to/beagle.jar
```

```R
beagler::run_beagle(gt = "my.vcf", out = "my_out.vcf.gz", n_threads=1)
```
