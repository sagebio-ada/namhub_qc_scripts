"""Small parsing helpers."""


def to_num(x):
    """'1,234' / '87.81%' / '49055009.0' -> int or float; '' stays ''."""
    s = str(x).replace(",", "").strip().rstrip("%").strip()
    if s == "":
        return ""
    try:
        return int(s)
    except ValueError:
        return float(s)


def pct(num, den):
    return num / den * 100 if den else ""


def read_table(path, comment="#"):
    """Tab-separated file with a header row -> list of dicts; lines starting with `comment` skipped."""
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip() and not ln.startswith(comment)]
    header = lines[0].split("\t")
    return [dict(zip(header, ln.split("\t"))) for ln in lines[1:]]
