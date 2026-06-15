import argparse
import re
from collections.abc import Sequence
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from .core import PmlstOptions, run
from .utils import PmlstError


def package_version() -> str:
    try:
        return version("pmlst_ssi")
    except PackageNotFoundError:
        pyproject = Path(__file__).resolve().parents[2] / "pyproject.toml"
        version_match = re.search(
            r'^version = "([^"]+)"',
            pyproject.read_text(encoding="utf-8"),
            re.MULTILINE,
        )
        if version_match is None:
            raise
        return version_match.group(1)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "-i",
        "--infile",
        help="FASTA or FASTQ files to do pMLST on.",
        nargs="+",
        required=True,
    )
    parser.add_argument("-o", "--outdir", help="Output directory.", default=".")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-s", "--scheme", help="scheme database used for pMLST prediction"
    )
    group.add_argument(
        "-pf", "--pf_results", help="Path to precomputed PlasmidFinder results"
    )
    parser.add_argument("-p", "--database", help="Directory containing the databases.")
    parser.add_argument(
        "-t",
        "--tmp_dir",
        help="Temporary directory for storage of the results\
                              from the external software.",
        default="tmp_pMLST",
    )
    parser.add_argument(
        "-c",
        "--coverage",
        help="Minimum template coverage threshold",
        type=float,
        default=0.6,
    )
    parser.add_argument(
        "-id",
        "--identity",
        help="Minimum template identity threshold",
        type=float,
        default=0.95,
    )
    parser.add_argument(
        "-mp",
        "--method_path",
        help=(
            "Path to the method to use (kma or blastn) if assembled contigs "
            "are inputted the path to executable blastn should be given, if "
            "fastq files are given path to executable kma should be given"
        ),
    )
    parser.add_argument(
        "-x",
        "--extented_output",
        help=(
            "Give extented output with allignment files, template and query "
            "hits in fasta and a tab seperated file with allele profile results"
        ),
        action="store_true",
    )
    parser.add_argument(
        "-xm",
        "--extented_output_with_scheme",
        help=(
            "Similar to -x (--extended_output) but it adds specific scheme "
            "name to the start of the file name when multiple schemes are "
            "selected"
        ),
        action="store_true",
    )
    parser.add_argument(
        "-so",
        "--simple_output",
        help=(
            "Give simple output with only the allele profile results in tab "
            "seperated file"
        ),
        action="store_true",
    )
    parser.add_argument("-q", "--quiet", action="store_true")
    parser.add_argument(
        "--version", action="version", version=f"pMLST {package_version()}"
    )
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    options = PmlstOptions(
        infile=args.infile,
        outdir=args.outdir,
        scheme=args.scheme,
        pf_results=args.pf_results,
        database=args.database,
        tmp_dir=args.tmp_dir,
        coverage=args.coverage,
        identity=args.identity,
        method_path=args.method_path,
        extented_output=args.extented_output,
        extented_output_with_scheme=args.extented_output_with_scheme,
        simple_output=args.simple_output,
        quiet=args.quiet,
    )
    try:
        run(options)
    except PmlstError as error:
        raise SystemExit(str(error)) from None


if __name__ == "__main__":
    main()
