---
title: Publish and maintain the Galick documentation
description: Build the MkDocs site, deploy to GitHub Pages and verify canonical URLs, public API coverage, sitemaps and search-engine indexing.
---
# Website publishing

The canonical documentation address is
<https://ndsu-agda.github.io/galick/>. The public source is
[NDSU-AGDA/galick](https://github.com/NDSU-AGDA/galick).

## GitHub Pages setup

The Python-only MkDocs workflow builds the site once, validates it and reuses
that artifact for deployment. Pull requests produce previews but cannot publish.
Changes to documentation on `main`, or a manual workflow run, deploy to Pages.
There is no R installation, platform matrix or second site build in this workflow.

The workflow attempts to enable Pages automatically. If organization permissions
prevent that, a repository administrator must open **Settings > Pages** and select
**GitHub Actions** as the build and deployment source, then rerun **MkDocs**.
If the `github-pages` environment requires approval, an administrator must approve
the deployment or adjust its protection rules. No personal access token or custom
domain is needed for the default project URL.

## Local verification

```sh
python -m pip install -r docs/requirements.txt
python -m mkdocs build --strict
python tools/check-site.py
Rscript --vanilla tools/build-reference.R
Rscript --vanilla tools/check-documentation.R
```

The site checker verifies unique page titles and descriptions, canonical and
social-sharing URLs, software metadata, sitemap coverage and all exported
function names in the search index. Update `extra.version` in `mkdocs.yml` when
the package version changes; the checker compares it with `DESCRIPTION`.

## Search discovery

Useful, accurate pages and public links support discovery; no configuration
guarantees a particular ranking or immediate indexing. The site uses distinct
descriptions and readable workflow URLs rather than keyword stuffing, fake
ratings or fabricated performance claims.

After publication, a maintainer can add a URL-prefix property for
`https://ndsu-agda.github.io/galick/` in Google Search Console. Ownership
verification requires the maintainer's account; add only the verification tag or
file supplied there. Submit `https://ndsu-agda.github.io/galick/sitemap.xml` and
inspect the homepage and key workflow URLs. Do not add a project-level robots
file expecting it to control the `ndsu-agda.github.io` host's root robots policy.

Keep the repository's About website URL pointed at the public documentation.
Relevant repository topics include R, VCF, genotype imputation, Beagle,
plant breeding, animal breeding and SNP panels. Only repository administrators
can edit those settings.
