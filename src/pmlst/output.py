import json
import pprint
import time
from pathlib import Path
from typing import Any, TextIO, TypeAlias

AlleleInfo: TypeAlias = dict[str, Any]
AlleleMatches: TypeAlias = dict[str, AlleleInfo]
Alignments: TypeAlias = dict[str, dict[str, str]]
JsonDict: TypeAlias = dict[str, Any]
Rows: TypeAlias = list[list[str]]


def make_allele_results(allele_matches: AlleleMatches) -> JsonDict:
    allele_results: JsonDict = {}
    for locus, locus_info in allele_matches.items():
        allele_results[locus] = {
            "identity": 0,
            "coverage": 0,
            "allele": [],
            "allele_name": [],
            "align_len": [],
            "gaps": 0,
            "sbj_len": [],
        }
        for key, value in locus_info.items():
            if key in allele_results[locus] or (
                key == "alternative_hit" and value != {}
            ):
                allele_results[locus][key] = value
    return allele_results


def build_service_data(
    infile: list[str],
    scheme: str,
    profile_name: str,
    file_format: str,
    st: str,
    allele_matches: AlleleMatches,
    nearest_sts: str,
    clpx: str,
    note: str,
) -> JsonDict:
    # Get run info for JSON file
    date = time.strftime("%d.%m.%Y")
    time_in_hours = time.strftime("%H:%M:%S")

    userinput = {
        "filename": infile,
        "scheme": scheme,
        "profile": profile_name,
        "file_format": file_format,
    }
    run_info = {
        "date": date,
        "time": time_in_hours,
    }  # , "database":{"remote_db":remote_db, "last_commit_hash":head_hash}}
    server_results = {
        "sequence_type": st,
        "allele_profile": make_allele_results(allele_matches),
        "nearest_sts": nearest_sts,
        "clonal_complex": clpx,
        "notes": note,
    }

    return {
        "user_input": userinput,
        "run_info": run_info,
        "results": server_results,
    }


def write_extended_outputs(
    outdir: str,
    schemes_to_run: list[str],
    extented_output_with_scheme: bool,
    service: str,
    scheme: str,
    profile_name: str,
    st: str,
    nearest_sts: str,
    clpx: str,
    allele_matches: AlleleMatches,
    query_aligns: Alignments,
    homol_aligns: Alignments,
    sbjct_aligns: Alignments,
    note: str,
) -> None:
    outdir_path = Path(outdir)
    # Define extented output
    if len(schemes_to_run) > 1 or extented_output_with_scheme:
        table_filename = outdir_path / f"{scheme}_results_tab.tsv"
        query_filename = outdir_path / f"{scheme}_Hit_in_genome_seq.fsa"
        sbjct_filename = outdir_path / f"{scheme}_pMLST_allele_seq.fsa"
        result_filename = outdir_path / f"{scheme}_results.txt"
    else:
        table_filename = outdir_path / "results_tab.tsv"
        query_filename = outdir_path / "Hit_in_genome_seq.fsa"
        sbjct_filename = outdir_path / "pMLST_allele_seq.fsa"
        result_filename = outdir_path / "results.txt"

    with (
        table_filename.open("w") as table_file,
        query_filename.open("w") as query_file,
        sbjct_filename.open("w") as sbjct_file,
        result_filename.open("w") as result_file,
    ):
        # Make results file
        result_file.write(f"{service} Results\n\n")
        result_file.write(f"pMLST profile: {profile_name}\n\nSequence Type: {st}\n")
        # If ST is unknown report nearest ST
        if st == "Unknown" and nearest_sts != "":
            if len(nearest_sts.split(",")) == 1:
                result_file.write(f"Nearest ST: {nearest_sts}\n")
            else:
                result_file.write(f"Nearest STs: {nearest_sts}\n")

        # Report clonal complex if one was associated with ST:
        if clpx != "":
            result_file.write(f"Clonal complex: {clpx}\n")

        # Write tsv table header
        table_header = [
            "Locus",
            "Identity",
            "Coverage",
            "Alignment Length",
            "Allele Length",
            "Gaps",
            "Allele",
        ]
        table_file.write("\t".join(table_header) + "\n")
        rows = []
        for locus, allele_info in allele_matches.items():
            identity = str(allele_info["identity"])
            coverage = str(allele_info["coverage"])
            allele = allele_info["allele"]
            allele_name = allele_info["allele_name"]
            align_len = str(allele_info["align_len"])
            sbj_len = str(allele_info["sbj_len"])
            gaps = str(allele_info["gaps"])

            # Write alleles names with indications of imperfect hits
            if allele_name != "No hit found":
                allele_name_w_mark = locus + "_" + allele
            else:
                allele_name_w_mark = allele_name

            # Write allele results to tsv table
            row = [
                locus,
                identity,
                coverage,
                align_len,
                sbj_len,
                gaps,
                allele_name_w_mark,
            ]
            rows.append(row)
            if "alternative_hit" in allele_info:
                for allele_name, dic in allele_info["alternative_hit"].items():
                    row = [
                        locus,
                        identity,
                        coverage,
                        str(dic["align_len"]),
                        str(dic["sbj_len"]),
                        "0",
                        allele_name + "!",
                    ]
                    rows.append(row)
            #

            if allele_name == "No hit found":
                continue

            # Write query fasta output
            hit_name = allele_info["hit_name"]
            query_seq = query_aligns[scheme][hit_name]
            sbjct_seq = sbjct_aligns[scheme][hit_name]

            match = "PERFECT MATCH" if allele_info["match_priority"] == 1 else "WARNING"
            header = (
                f">{locus}:{match} ID:{allele_info['identity']}% "
                f"COV:{allele_info['coverage']}% "
                f"Best_match:{allele_info['allele_name']}\n"
            )
            query_file.write(header)
            for i in range(0, len(query_seq), 60):
                query_file.write(query_seq[i : i + 60] + "\n")

            # Write template fasta output
            header = f">{allele_info['allele_name']}\n"
            sbjct_file.write(header)
            for i in range(0, len(sbjct_seq), 60):
                sbjct_file.write(sbjct_seq[i : i + 60] + "\n")

            if "alternative_hit" in allele_info:
                for allele_name in allele_info["alternative_hit"]:
                    header = (
                        f">{locus}:PERFECT MATCH ID:100% COV:100% "
                        f"Best_match:{allele_name}\n"
                    )
                    hit_name = allele_info["alternative_hit"][allele_name]["hit_name"]
                    query_seq = query_aligns[scheme][hit_name]
                    sbjct_seq = sbjct_aligns[scheme][hit_name]
                    query_file.write(header)
                    for i in range(0, len(query_seq), 60):
                        query_file.write(query_seq[i : i + 60] + "\n")

                    # Write template fasta output
                    header = f">{allele_name}\n"
                    sbjct_file.write(header)
                    for i in range(0, len(sbjct_seq), 60):
                        sbjct_file.write(sbjct_seq[i : i + 60] + "\n")

        # Write Allele profile results tables in results file and table file
        rows.sort(key=lambda x: x[0])
        result_file.write(text_table(table_header, rows))
        for row in rows:
            table_file.write("\t".join(row) + "\n")
        # Write any notes
        if note != "":
            result_file.write(f"\nNotes: {note}\n\n")

        # Write allignment output
        result_file.write("\n\nExtended Output:\n\n")
        make_aln(
            scheme,
            result_file,
            allele_matches,
            query_aligns,
            homol_aligns,
            sbjct_aligns,
        )


def write_json_output(outdir: str, data: JsonDict) -> None:
    data_result_file = Path(outdir) / "data.json"
    with data_result_file.open("w") as outfile:
        json.dump(data, outfile)

    # Legacy CLI behavior prints generated results to stdout.
    pprint.pprint(data)


def build_simple_output_rows(
    data: JsonDict, list_of_plasmids: list[str], pf_results: str | None
) -> Rows:
    list_of_plasmids_str = ",".join(list_of_plasmids)
    simple_output_list = [
        [
            "plasmids",
            "IncF",
            "IncI1",
            "IncA/C",
            "IncHI1",
            "IncHI2",
            "IncN",
            "pMLST summary",
        ],
        [list_of_plasmids_str, "", "", "", "", "", "", ""],
    ]

    for service in data:
        scheme = data[service]["user_input"]["scheme"].lower()

        sequence_type = data[service]["results"]["sequence_type"]
        if sequence_type == "Unknown":
            sequence_type = ""
        if sequence_type == "":
            continue
        if scheme == "incf":
            simple_output_list[1][1] = sequence_type.strip("[]")
        elif scheme == "inci1":
            simple_output_list[1][2] = sequence_type.strip("[]")
        elif scheme == "incac":
            simple_output_list[1][3] = sequence_type.strip("[]")
        elif scheme == "inchi1":
            simple_output_list[1][4] = sequence_type.strip("[]")
        elif scheme == "inchi2":
            simple_output_list[1][5] = sequence_type.strip("[]")
        elif scheme == "incn":
            simple_output_list[1][6] = sequence_type.strip("[]")
        if simple_output_list[1][7] == "":
            simple_output_list[1][7] = scheme + sequence_type
            if not pf_results:
                simple_output_list[1][0] = scheme
        else:
            if sequence_type != "":
                simple_output_list[1][7] += "," + scheme + sequence_type
                if not pf_results:
                    simple_output_list[1][0] += "," + scheme

    return simple_output_list


def write_simple_output(
    outdir: str, data: JsonDict, list_of_plasmids: list[str], pf_results: str | None
) -> None:
    simple_output_list = build_simple_output_rows(data, list_of_plasmids, pf_results)
    simple_output_file = Path(outdir) / "simple_output.tsv"
    write_rows(simple_output_file, simple_output_list)

    # Legacy CLI behavior prints generated results to stdout.
    pprint.pprint(simple_output_list)


def write_empty_simple_output(outdir: str) -> None:
    simple_output_list = [
        [
            "plasmids",
            "IncF",
            "IncI1",
            "IncA/C",
            "IncHI1",
            "IncHI2",
            "IncN",
            "pMLST summary",
        ],
        ["", "", "", "", "", "", "", ""],
    ]
    simple_output_file = Path(outdir) / "simple_output.tsv"
    write_rows(simple_output_file, simple_output_list)


def write_rows(filename: Path, rows: Rows) -> None:
    with filename.open("w") as file_handle:
        for row in rows:
            file_handle.write("\t".join(row) + "\n")


def make_aln(
    scheme: str,
    file_handle: TextIO,
    allele_matches: AlleleMatches,
    query_aligns: Alignments,
    homol_aligns: Alignments,
    sbjct_aligns: Alignments,
) -> None:
    for locus_info in allele_matches.values():
        allele_name = locus_info["allele_name"]
        if allele_name == "No hit found":
            continue
        hit_name = locus_info["hit_name"]

        seqs = ["", "", ""]
        seqs[0] = sbjct_aligns[scheme][hit_name]
        seqs[1] = homol_aligns[scheme][hit_name]
        seqs[2] = query_aligns[scheme][hit_name]

        write_align(seqs, allele_name, file_handle)

        # write alternative seq
        if "alternative_hit" in locus_info:
            for allele_name in locus_info["alternative_hit"]:
                hit_name = locus_info["alternative_hit"][allele_name]["hit_name"]
                seqs = ["", "", ""]
                seqs[0] = sbjct_aligns[scheme][hit_name]
                seqs[1] = homol_aligns[scheme][hit_name]
                seqs[2] = query_aligns[scheme][hit_name]

                write_align(seqs, allele_name, file_handle)


def write_align(seq: list[str], seq_name: str, file_handle: TextIO) -> None:
    file_handle.write(f"# {seq_name}\n")
    sbjct_seq = seq[0]
    homol_seq = seq[1]
    query_seq = seq[2]
    for i in range(0, len(sbjct_seq), 60):
        file_handle.write(f"{'template:':<10}\t{sbjct_seq[i : i + 60]}\n")
        file_handle.write(f"{'':<10}\t{homol_seq[i : i + 60]}\n")
        file_handle.write(f"{'query:':<10}\t{query_seq[i : i + 60]}\n\n")


def text_table(headers: list[str], rows: Rows, empty_replace: str = "-") -> str:
    """Create text table

    USAGE:
       >>> from tabulate import tabulate
       >>> headers = ['A','B']
       >>> rows = [[1,2],[3,4]]
       >>> print(text_table(headers, rows))
       **********
         A    B
       **********
         1    2
         3    4
       ==========
    """
    # Lazy import keeps pmlst --help/--version available without runtime deps.
    from tabulate import tabulate  # noqa: PLC0415

    # Replace empty cells with placeholder
    normalized_rows = (
        (value if value else empty_replace for value in row) for row in rows
    )
    # Create table
    table = tabulate(normalized_rows, headers, tablefmt="simple").split("\n")
    # Prepare title injection
    width = len(table[0])
    # Switch horisontal line
    table[1] = "*" * (width + 2)
    # Update table with title
    return ("%s\n" * 3) % ("*" * (width + 2), "\n".join(table), "=" * (width + 2))
