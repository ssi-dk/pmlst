#!/usr/bin/env python3

import os
import re
from collections.abc import Sequence
from contextlib import redirect_stdout
from dataclasses import dataclass
from pathlib import Path
from typing import Any, TypeAlias

from .database import (
    resolve_database_path,
    validate_database_root,
    validate_database_schemes,
)
from .methods import run_methods
from .output import (
    build_service_data,
    write_empty_simple_output,
    write_extended_outputs,
    write_json_output,
    write_simple_output,
)
from .utils import (
    PmlstError,
    extended_cigar,
    get_file_format,
    import_profile,
    load_scheme_list_config,
    plasmidfinder_parsing,
    st_typing,
)

SchemeList: TypeAlias = dict[str, str]
AlleleInfo: TypeAlias = dict[str, Any]
AlleleMatches: TypeAlias = dict[str, AlleleInfo]
MethodResults: TypeAlias = dict[str, Any]
Alignments: TypeAlias = dict[str, dict[str, str]]


@dataclass(frozen=True)
class PmlstOptions:
    infile: Sequence[str]
    outdir: str = "."
    scheme: str | None = None
    pf_results: str | None = None
    database: str | None = None
    tmp_dir: str = "tmp_pMLST"
    coverage: float = 0.6
    identity: float = 0.95
    method_path: str | None = None
    extented_output: bool = False
    extented_output_with_scheme: bool = False
    simple_output: bool = False
    quiet: bool = False


def run(options: PmlstOptions) -> None:
    if options.quiet:
        with Path(os.devnull).open("w") as quiet_output, redirect_stdout(quiet_output):
            _run(options)
        return

    _run(options)


def _run(options: PmlstOptions) -> None:
    infile = list(options.infile)
    for input_file in infile:
        if not Path(input_file).exists():
            raise PmlstError(f"Input file does not exist: {input_file}")

    outdir_path = Path(options.outdir).resolve()
    outdir = str(outdir_path)
    if not outdir_path.exists():
        raise PmlstError(f"Output folder does not exist: {outdir}")

    database_path = resolve_database_path(options.database)
    validate_database_root(database_path)
    database = str(database_path)

    scheme_list = load_scheme_list_config({}, database)

    schemes_to_run, profile_names_to_run, list_of_plasmids = resolve_schemes(
        options, scheme_list, outdir
    )
    if schemes_to_run == []:
        return

    tmp_dir = str(Path(options.tmp_dir).resolve())

    if 0 <= options.coverage <= 1:
        min_cov = options.coverage
    else:
        raise PmlstError("Coverage threshold must be between 0 and 1.")
    if 0 <= options.identity <= 1:
        threshold = options.identity
    else:
        raise PmlstError("Identity threshold must be between 0 and 1.")

    # Check file format (fasta, fastq or other format)
    file_format = get_file_format(infile)
    validate_database_schemes(
        database_path, schemes_to_run, require_kma=file_format == "fastq"
    )

    (
        list_of_loci_list,
        results_list,
        query_aligns_list,
        homol_aligns_list,
        sbjct_aligns_list,
    ) = run_methods(
        infile,
        outdir,
        schemes_to_run,
        database,
        tmp_dir,
        file_format,
        min_cov,
        threshold,
        options.method_path,
    )

    data: dict[str, Any] = {}
    for i, scheme in enumerate(schemes_to_run):
        loci_list = list_of_loci_list[i]
        profile_name = profile_names_to_run[i]
        results = results_list[i]
        query_aligns = query_aligns_list[i]
        homol_aligns = homol_aligns_list[i]
        sbjct_aligns = sbjct_aligns_list[i]

        allele_matches, warning = type_alleles(
            scheme,
            loci_list,
            results,
            query_aligns,
            homol_aligns,
            sbjct_aligns,
        )

        # Import all possible st profiles into dict
        st_profiles = import_profile(database, scheme, loci_list)

        # Find st or neatest sts
        st, note, nearest_sts = st_typing(st_profiles, allele_matches, loci_list)

        # Give warning of mlst schene if no loci were found
        if note == "" and warning != "":
            note = warning

        # Set ST for incF
        if scheme.lower() == "incf":
            st = format_incf_st(allele_matches)

        clpx = find_clonal_complex(database, scheme, st, nearest_sts)

        # Get run info for JSON file
        service = "pmlst"
        if i > 0:
            service += f"_{i}"

        # Make JSON output file
        data[service] = build_service_data(
            infile,
            scheme,
            profile_name,
            file_format,
            st,
            allele_matches,
            nearest_sts,
            clpx,
            note,
        )

        if options.extented_output or options.extented_output_with_scheme:
            write_extended_outputs(
                outdir,
                schemes_to_run,
                options.extented_output_with_scheme,
                service,
                scheme,
                profile_name,
                st,
                nearest_sts,
                clpx,
                allele_matches,
                query_aligns,
                homol_aligns,
                sbjct_aligns,
                note,
            )

    # Save json output
    if not options.simple_output:
        write_json_output(outdir, data)
    else:
        write_simple_output(outdir, data, list_of_plasmids, options.pf_results)


def resolve_schemes(
    options: PmlstOptions, scheme_list: SchemeList, outdir: str
) -> tuple[list[str], list[str], list[str]]:
    schemes_to_run: list[str] = []
    profile_names_to_run: list[str] = []
    list_of_plasmids: list[str] = []
    if options.pf_results:
        schemes_to_run, profile_names_to_run, list_of_plasmids = plasmidfinder_parsing(
            options.pf_results, scheme_list
        )
        if schemes_to_run == []:
            write_empty_simple_output(outdir)
            print("PlasmidFinder results file does not contain any data, only header.")
            return [], [], []
    elif options.scheme == "all":
        schemes_to_run = list(scheme_list.keys())
        profile_names_to_run = list(scheme_list.values())
    elif options.scheme and options.scheme != "all":
        for each_scheme in options.scheme.split(","):
            if each_scheme in scheme_list:
                schemes_to_run.append(each_scheme)
                profile_names_to_run.append(scheme_list[each_scheme])
            else:
                raise PmlstError(
                    f"{options.scheme}, is not a valid scheme. \n\n"
                    "Please choose a scheme available in the database:\n"
                    f"{', '.join(scheme_list)}"
                )
    else:
        raise PmlstError(
            "No schemes provided, please provide a scheme to run pMLST on with -s flag."
        )

    return schemes_to_run, profile_names_to_run, list_of_plasmids


def type_alleles(
    scheme: str,
    loci_list: list[str],
    results: MethodResults,
    query_aligns: Alignments,
    homol_aligns: Alignments,
    sbjct_aligns: Alignments,
) -> tuple[AlleleMatches, str]:
    # Check that the results dict is not empty
    warning = ""
    if results[scheme] == "No hit found":
        results[scheme] = {}
        warning = (
            "No MLST loci was found in the input data, "
            "make sure that the correct pMLST scheme was chosen."
        )

    allele_matches: AlleleMatches = {}

    # Get the found allele profile contained in the results dict
    for hit, locus_hit in results[scheme].items():
        # Get allele number for locus
        allele_name = locus_hit["sbjct_header"]
        allele_obj = re.search(r"(\w+)[_|-](\w+$)", allele_name)
        if allele_obj is None:
            raise PmlstError(f"Unable to parse allele name: {allele_name}")

        # Get variable to later storage in the results dict
        locus = allele_obj.group(1)
        allele = allele_obj.group(2)
        coverage = float(locus_hit["perc_coverage"])
        identity = float(locus_hit["perc_ident"])
        score = float(locus_hit["cal_score"])
        gaps = int(locus_hit["gaps"])
        align_len = locus_hit["HSP_length"]
        sbj_len = int(locus_hit["sbjct_length"])
        sbjct_seq = locus_hit["sbjct_string"]
        query_seq = locus_hit["query_string"]
        homol_seq = locus_hit["homo_string"]
        cigar = extended_cigar(sbjct_aligns[scheme][hit], query_aligns[scheme][hit])

        if coverage == 100 and identity == 100:
            try:
                allele_matches[locus]["alternative_hit"][allele_name] = {
                    "allele": allele + "!",
                    "align_len": align_len,
                    "sbj_len": sbj_len,
                    "coverage": coverage,
                    "identity": identity,
                    "hit_name": hit,
                }
                if allele_matches[locus]["allele"][-1] != "!":
                    allele_matches[locus]["allele"] += "!"
            except KeyError:
                allele_matches[locus] = {
                    "score": score,
                    "allele": allele,
                    "coverage": coverage,
                    "identity": identity,
                    "match_priority": 1,
                    "align_len": align_len,
                    "gaps": gaps,
                    "sbj_len": sbj_len,
                    "allele_name": allele_name,
                    "sbjct_seq": sbjct_seq,
                    "query_seq": query_seq,
                    "homol_seq": homol_seq,
                    "hit_name": hit,
                    "cigar": cigar,
                    "alternative_hit": {},
                }
        else:
            if locus not in allele_matches:
                allele_matches[locus] = {"score": 0, "match_priority": 4}

            if coverage == 100 and identity != 100:
                if allele_matches[locus]["match_priority"] > 2 or (
                    allele_matches[locus]["match_priority"] == 2
                    and score > allele_matches[locus]["score"]
                ):
                    allele_matches[locus] = {
                        "score": score,
                        "allele": allele + "*",
                        "coverage": coverage,
                        "identity": identity,
                        "match_priority": 2,
                        "align_len": align_len,
                        "gaps": gaps,
                        "sbj_len": sbj_len,
                        "allele_name": allele_name,
                        "sbjct_seq": sbjct_seq,
                        "query_seq": query_seq,
                        "homol_seq": homol_seq,
                        "hit_name": hit,
                        "cigar": cigar,
                    }
            elif coverage != 100 and identity == 100:
                # Check that higher prioritized hit was not already stored
                if allele_matches[locus]["match_priority"] > 3 or (
                    allele_matches[locus]["match_priority"] == 3
                    and score > allele_matches[locus]["score"]
                ):
                    allele_matches[locus] = {
                        "score": score,
                        "allele": allele + "?",
                        "coverage": coverage,
                        "identity": identity,
                        "match_priority": 3,
                        "align_len": align_len,
                        "gaps": gaps,
                        "sbj_len": sbj_len,
                        "allele_name": allele_name,
                        "sbjct_seq": sbjct_seq,
                        "query_seq": query_seq,
                        "homol_seq": homol_seq,
                        "hit_name": hit,
                        "cigar": cigar,
                    }
            else:  # coverage != 100 and identity != 100:
                if (
                    allele_matches[locus]["match_priority"] == 4
                    and score > allele_matches[locus]["score"]
                ):
                    allele_matches[locus] = {
                        "score": score,
                        "allele": allele + "?*",
                        "coverage": coverage,
                        "identity": identity,
                        "match_priority": 4,
                        "align_len": align_len,
                        "gaps": gaps,
                        "sbj_len": sbj_len,
                        "allele_name": allele_name,
                        "sbjct_seq": sbjct_seq,
                        "query_seq": query_seq,
                        "homol_seq": homol_seq,
                        "hit_name": hit,
                        "cigar": cigar,
                    }
    for locus in loci_list:
        if locus not in allele_matches:
            allele_matches[locus] = {
                "identity": "",
                "coverage": "",
                "allele": "",
                "allele_name": "No hit found",
                "align_len": "",
                "gaps": "",
                "sbj_len": "",
            }

    return allele_matches, warning


def format_incf_st(allele_matches: AlleleMatches) -> str:
    st = ["F", "A", "B"]
    if "FII" in allele_matches and allele_matches["FII"]["identity"] == 100.0:
        st[0] += allele_matches["FII"]["allele_name"].split("_")[-1]
    elif "FIC" in allele_matches and allele_matches["FIC"]["identity"] == 100.0:
        st[0] = "C" + allele_matches["FIC"]["allele_name"].split("_")[-1]
    elif "FIIK" in allele_matches and allele_matches["FIIK"]["identity"] == 100.0:
        st[0] = "K" + allele_matches["FIIK"]["allele_name"].split("_")[-1]
    elif "FIIS" in allele_matches and allele_matches["FIIS"]["identity"] == 100.0:
        st[0] = "S" + allele_matches["FIIS"]["allele_name"].split("_")[-1]
    elif "FIIY" in allele_matches and allele_matches["FIIY"]["identity"] == 100.0:
        st[0] = "Y" + allele_matches["FIIY"]["allele_name"].split("_")[-1]
    else:
        st[0] += "-"

    if "FIA" in allele_matches and allele_matches["FIA"]["identity"] == 100.0:
        st[1] += allele_matches["FIA"]["allele_name"].split("_")[-1]
    else:
        st[1] += "-"

    if "FIB" in allele_matches and allele_matches["FIB"]["identity"] == 100.0:
        st[2] += allele_matches["FIB"]["allele_name"].split("_")[-1]
    else:
        st[2] += "-"

    return "[" + ":".join(st) + "]"


def find_clonal_complex(database: str, scheme: str, st: str, nearest_sts: str) -> str:
    # Check if ST is associated with a clonal complex.
    clpx = ""
    if st != "Unknown" or nearest_sts != "":
        clpx_path = Path(database) / f"{scheme}.clpx"
        with clpx_path.open() as clpx_file:
            for raw_line in clpx_file:
                line = raw_line.split("\t")
                if st[0] == line[0] or nearest_sts == line[0]:
                    clpx = line[1].strip()
    return clpx
