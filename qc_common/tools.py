"""Run external QC tools so their output lands where fastq_sections.py reads it.

Needs on PATH: `fq` (stjude-rust-labs), `kraken2`, `bracken` (e.g. conda env namhub_atac_qc, or
`module load kraken/... bracken/...`). Every run_* function skips work whose output already
exists, so a rerun resumes, and writes through a temp file so a failed run leaves no output.
"""
import gzip
import re
import shutil
import subprocess
import sys
import tempfile
from itertools import islice
from pathlib import Path


def need(*tools):
    """Exit with a clear message if any tool is not on PATH."""
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        sys.exit(f"Not on PATH: {', '.join(missing)} (see qc_common/tools.py for how to get them)")


def _run(cmd, log=None):
    res = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    if log is not None:
        Path(log).write_text(res.stdout + res.stderr)
    if res.returncode != 0:
        sys.exit(f"Command failed ({res.returncode}): {' '.join(str(c) for c in cmd)}\n{res.stderr[-2000:]}")


def run_fq_lint(fastqs, out, fq="fq"):
    """`fq lint` on one FASTQ (single-end) or a pair, log -> `out` (<root>/fq_lint/<stage>/<sample>.fq_lint.txt).

    Same options as nf-core (paired-name validator P001 off); fq logs to stderr.
    """
    out = Path(out)
    if out.exists():
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.name + ".tmp")
    # fq exits non-zero when it finds invalid records; the log (with the ERROR lines) is the result.
    res = subprocess.run([fq, "lint", "--disable-validator", "P001", *map(str, fastqs)],
                         capture_output=True, text=True)
    tmp.write_text(res.stdout + res.stderr)
    tmp.rename(out)
    return out


def _head_fastq(src, dst, n_reads):
    """First n_reads records of a (gzipped) FASTQ -> plain FASTQ."""
    opener = gzip.open if str(src).endswith(".gz") else open
    with opener(src, "rt") as fh, open(dst, "w") as out:
        out.writelines(islice(fh, 4 * n_reads))


def bracken_length(db, read_length):
    """Bracken needs database<N>mers.kmer_distrib; pick the N closest to the read length."""
    lens = sorted(int(m[1]) for f in Path(db).glob("database*mers.kmer_distrib")
                  if (m := re.fullmatch(r"database(\d+)mers\.kmer_distrib", f.name)))
    if not lens:
        sys.exit(f"No database*mers.kmer_distrib in {db}: this Kraken2 database has no Bracken files")
    return min(lens, key=lambda n: abs(n - read_length))


def run_kraken_bracken(fastqs, sample, contam_dir, db, read_length, threads=8, max_reads=0,
                       kraken2="kraken2", bracken="bracken"):
    """kraken2 then bracken (species level) for one sample.

    Writes the layout fastq_sections.contaminants reads:
      <contam_dir>/kraken2/kraken_reports/<sample>.kraken2.report.txt
      <contam_dir>/bracken/<sample>.tsv
    fastqs: one file (single-end) or two (paired). max_reads > 0 screens only the first
    max_reads reads (pairs) of each file; 0 = all.
    """
    contam_dir = Path(contam_dir)
    report = contam_dir / "kraken2" / "kraken_reports" / f"{sample}.kraken2.report.txt"
    tsv = contam_dir / "bracken" / f"{sample}.tsv"
    report.parent.mkdir(parents=True, exist_ok=True)
    tsv.parent.mkdir(parents=True, exist_ok=True)

    if not report.exists():
        with tempfile.TemporaryDirectory(prefix=f"kraken_{sample}_", dir=contam_dir) as tmp:
            inputs = [str(f) for f in fastqs]
            if max_reads > 0:
                inputs = []
                for i, f in enumerate(fastqs, 1):
                    dst = Path(tmp) / f"r{i}.fastq"
                    _head_fastq(f, dst, max_reads)
                    inputs.append(str(dst))
            partial = report.with_name(report.name + ".tmp")
            cmd = [kraken2, "--db", db, "--threads", threads, "--report", partial,
                   "--output", "/dev/null"]
            if len(inputs) == 2:
                cmd.append("--paired")
            _run(cmd + inputs, log=contam_dir / "kraken2" / f"{sample}.kraken2.log")
            partial.rename(report)

    if not tsv.exists():
        partial = tsv.with_name(tsv.name + ".tmp")
        _run([bracken, "-d", db, "-i", report, "-o", partial, "-r", bracken_length(db, read_length),
              "-l", "S", "-t", "10"], log=contam_dir / "bracken" / f"{sample}.bracken.log")
        partial.rename(tsv)
    return report, tsv
