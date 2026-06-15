import shutil
from pathlib import Path
from typing import Any, TypeAlias

from .utils import PmlstError, get_read_filename

MethodResults: TypeAlias = dict[str, Any]
Alignments: TypeAlias = dict[str, dict[str, str]]


def resolve_method_path(file_format: str, method_path: str | None) -> str:
    if method_path:
        resolved_method_path = shutil.which(method_path)
        if resolved_method_path is None:
            raise PmlstError(
                "Invalid -mp/--method_path executable path: "
                f"{method_path}. Provide a valid path to blastn or kma."
            )
        return resolved_method_path

    if file_format == "fasta":
        executable = "blastn"
        input_mode = "FASTA/assembled-contig input"
    elif file_format == "fastq":
        executable = "kma"
        input_mode = "FASTQ/read input"
    else:
        raise PmlstError("Input file must be fastq or fasta format, not " + file_format)

    resolved_method_path = shutil.which(executable)
    if resolved_method_path is None:
        raise PmlstError(
            f"Missing required runtime executable '{executable}' for {input_mode}. "
            f"Install {executable} or provide its path with -mp/--method_path."
        )
    return resolved_method_path


def run_methods(
    infile: list[str],
    outdir: str,
    schemes_to_run: list[str],
    database: str,
    tmp_dir: str,
    file_format: str,
    min_cov: float,
    threshold: float,
    method_path: str | None,
) -> tuple[
    list[list[str]],
    list[MethodResults],
    list[Alignments],
    list[Alignments],
    list[Alignments],
]:
    db_path = f"{database}/"
    resolved_method_path = resolve_method_path(file_format, method_path)

    list_of_method_obj: list[Any] = []
    list_of_loci_list: list[list[str]] = []
    database_path = Path(database)
    for scheme in schemes_to_run:
        # Get loci list from allele profile file
        with (database_path / f"{scheme}.txt.clean").open() as st_file:
            file_header = st_file.readline().strip().split("\t")
            list_of_loci_list.append(file_header[1:])

        # Call appropriate method (kma or blastn) based on file format
        if file_format == "fastq":
            from cgecore.cgefinder import CGEFinder

            # Check the number of files
            if len(infile) == 1:
                infile_1 = infile[0]
                infile_2 = None
            elif len(infile) == 2:
                infile_1 = infile[0]
                infile_2 = infile[1]
            else:
                raise PmlstError("Only 2 input file accepted for raw read data,\
                        if data from more runs is avaliable for the same\
                        sample, please concatinate the reads into two files")

            sample_name = get_read_filename(infile)

            # Call KMA
            method_obj = CGEFinder.kma(
                infile_1,
                outdir,
                [scheme],
                db_path,
                min_cov=min_cov,
                threshold=threshold,
                kma_path=resolved_method_path,
                sample_name=sample_name,
                inputfile_2=infile_2,
                kma_mrs=0.75,
                kma_gapopen=-5,
                kma_gapextend=-1,
                kma_penalty=-3,
                kma_reward=1,
            )

            list_of_method_obj.append(method_obj)

        elif file_format == "fasta":
            from cgecore.blaster.blaster import Blaster

            # Check that only one fasta file is inputted
            if len(infile) != 1:
                raise PmlstError("Only one input file accepted for assembled data.")
            infile_1 = infile[0]
            # Call BLASTn
            method_obj = Blaster(
                infile_1,
                [scheme],
                db_path,
                tmp_dir,
                min_cov,
                threshold,
                resolved_method_path,
                cut_off=False,
            )
            # allewed_overlap=50)

            list_of_method_obj.append(method_obj)
    return collect_method_results(list_of_method_obj, list_of_loci_list)


def collect_method_results(
    list_of_method_obj: list[Any], list_of_loci_list: list[list[str]]
) -> tuple[
    list[list[str]],
    list[MethodResults],
    list[Alignments],
    list[Alignments],
    list[Alignments],
]:
    # Get the results from the method objects
    results_list: list[MethodResults] = []
    query_aligns_list: list[Alignments] = []
    homol_aligns_list: list[Alignments] = []
    sbjct_aligns_list: list[Alignments] = []
    for method_obj in list_of_method_obj:
        results_list.append(method_obj.results)
        query_aligns_list.append(method_obj.gene_align_query)
        homol_aligns_list.append(method_obj.gene_align_homo)
        sbjct_aligns_list.append(method_obj.gene_align_sbjct)

    return (
        list_of_loci_list,
        results_list,
        query_aligns_list,
        homol_aligns_list,
        sbjct_aligns_list,
    )
