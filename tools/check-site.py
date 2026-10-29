"""Check the built documentation's canonical URLs, search metadata and API coverage."""

import json
import sys
from email.parser import Parser
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urljoin
from xml.etree import ElementTree

BASE = "https://ndsu-agda.github.io/galick/"


class Head(HTMLParser):
    def __init__(self):
        super().__init__()
        self.title = ""
        self.in_title = False
        self.in_schema = False
        self.schema = ""
        self.meta = {}
        self.canonical = []
        self.h1 = 0

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "title":
            self.in_title = True
        if tag == "h1":
            self.h1 += 1
        if tag == "meta":
            key = attrs.get("name", attrs.get("property"))
            self.meta.setdefault(key, []).append(attrs.get("content", ""))
        if tag == "link" and attrs.get("rel") == "canonical":
            self.canonical.append(attrs.get("href"))
        if tag == "script" and attrs.get("type") == "application/ld+json":
            self.in_schema = True

    def handle_endtag(self, tag):
        if tag == "title":
            self.in_title = False
        if tag == "script":
            self.in_schema = False

    def handle_data(self, data):
        if self.in_title:
            self.title += data
        if self.in_schema:
            self.schema += data


def check(directory):
    titles, descriptions, urls = set(), set(), set()
    for path in sorted(directory.rglob("*.html")):
        relative = path.relative_to(directory).as_posix()
        if relative == "404.html":
            continue
        if relative.startswith(("overrides/", "about/")):
            raise AssertionError("Theme source and internal maintainer pages must not be published")
        head = Head()
        head.feed(path.read_text(encoding="utf-8"))
        url = urljoin(BASE, relative.removesuffix("index.html"))
        assert head.canonical == [url], f"{relative}: incorrect canonical URL"
        assert head.title and head.title not in titles, f"{relative}: missing/duplicate title"
        assert head.h1 == 1, f"{relative}: expected one main heading"
        description = head.meta.get("description", [])
        assert len(description) == 1 and 40 <= len(description[0]) <= 200, relative
        assert description[0] not in descriptions, f"{relative}: duplicate description"
        assert head.meta.get("og:url") == [url], f"{relative}: incorrect sharing URL"
        assert head.meta.get("og:description") == description, relative
        assert head.meta.get("twitter:card") == ["summary"], relative
        assert "noindex" not in ",".join(head.meta.get("robots", [])), relative
        titles.add(head.title)
        descriptions.add(description[0])
        urls.add(url)
        if relative == "index.html":
            schema = json.loads(head.schema)
            version = Parser().parsestr(Path("DESCRIPTION").read_text())["Version"]
            assert schema["@type"] == "SoftwareSourceCode"
            assert schema["version"] == version, "Schema version must match DESCRIPTION"
            assert schema["codeRepository"] == "https://github.com/NDSU-AGDA/galick"
    sitemap = ElementTree.parse(directory / "sitemap.xml")
    locations = {node.text for node in sitemap.findall(".//{*}loc")}
    assert locations == urls, "Sitemap must cover every indexable page at the canonical origin"
    search = json.loads((directory / "search/search_index.json").read_text())
    indexed = " ".join(row["text"] for row in search["docs"])
    exports = [line[7:-1] for line in Path("NAMESPACE").read_text().splitlines()
               if line.startswith("export(")]
    assert all(name in indexed for name in exports), "Search index is missing public functions"
    assert (directory / "assets/logo.png").is_file(), "Missing social preview image"
    assert urls, "No built documentation found"
    print(f"Checked {len(urls)} pages: unique metadata, canonical URLs, sitemap, sharing and {len(exports)} exports.")


if __name__ == "__main__":
    check(Path(sys.argv[1]) if len(sys.argv) > 1 else Path("site"))
