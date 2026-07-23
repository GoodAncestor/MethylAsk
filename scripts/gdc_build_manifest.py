"""Generate a GDC download manifest for the full methylation-beta corpus.

Produces a TSV of every "Methylation Beta Value" file (id, name, size, project,
md5) so the 294 GB pull becomes one resumable job on the server rather than
something discovered mid-run (docs/DESIGN.md §8.2).

Run on the server that will hold the mirror:
    python scripts/gdc_build_manifest.py --out data/mirror/gdc/manifest.tsv
Then download with the GDC client or curl loop against the manifest.

Output is also a valid GDC Data Transfer Tool manifest (id/filename/md5/size/state).
"""
from __future__ import annotations
import argparse, json, sys, urllib.parse, urllib.request

_API = "https://api.gdc.cancer.gov/files"
_FILTER = {"op": "in", "content": {"field": "data_type", "value": ["Methylation Beta Value"]}}
_FIELDS = "file_id,file_name,file_size,md5sum,state,cases.project.project_id"


def page(from_: int, size: int) -> dict:
    params = {"filters": json.dumps(_FILTER), "fields": _FIELDS,
              "format": "json", "size": str(size), "from": str(from_)}
    url = _API + "?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "methylask"}), timeout=120) as r:
        return json.loads(r.read())


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="gdc_manifest.tsv")
    ap.add_argument("--page-size", type=int, default=2000)
    args = ap.parse_args(argv)

    from_, total, rows, bytes_total = 0, None, [], 0
    while True:
        j = page(from_, args.page_size)
        pg = j["data"]["pagination"]
        total = pg["total"]
        for h in j["data"]["hits"]:
            proj = (h.get("cases") or [{}])[0].get("project", {}).get("project_id", "NA")
            rows.append((h["file_id"], h["file_name"], h.get("md5sum", "NA"),
                         h["file_size"], h.get("state", "NA"), proj))
            bytes_total += h["file_size"]
        from_ += args.page_size
        sys.stderr.write(f"\r{len(rows)}/{total} files, {bytes_total/1e9:.1f} GB")
        if from_ >= total:
            break
    sys.stderr.write("\n")

    with open(args.out, "w") as fh:
        fh.write("id\tfilename\tmd5\tsize\tstate\tproject_id\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print(f"wrote {len(rows)} files ({bytes_total/1e9:.1f} GB) -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
