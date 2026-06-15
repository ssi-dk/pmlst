import gzip
from collections.abc import Sequence
from pathlib import Path
from typing import Any, TypeAlias

AlleleInfo: TypeAlias = dict[str, Any]
AlleleMatches: TypeAlias = dict[str, AlleleInfo]
StProfiles: TypeAlias = dict[str, dict[str, list[str]]]
SchemeList: TypeAlias = dict[str, str]


class PmlstError(Exception):
    def __init__(self, message: str, exit_code: int = 1) -> None:
        super().__init__(message)
        self.exit_code = exit_code


def get_read_filename(infiles: Sequence[str]) -> str:
    """Infiles must be a list with 1 or 2 input files.
    Removes path from given string and removes extensions:
    .fq .fastq .gz and .trim
    extract the common sample name i 2 files are given.
    """
    # Remove common fastq extensions
    seq_path = infiles[0]
    seq_file = Path(seq_path).name
    seq_file = seq_file.replace(".fq", "")
    seq_file = seq_file.replace(".fastq", "")
    seq_file = seq_file.replace(".gz", "")
    seq_file = seq_file.replace(".trim", "")
    if len(infiles) == 1:
        return seq_file.rstrip()

    # If two files are given get the common sample name
    sample_name = ""
    seq_file_2 = Path(infiles[1]).name
    for i in range(len(seq_file)):
        if seq_file_2[i] == seq_file[i]:
            sample_name += seq_file[i]
        else:
            break
    if sample_name == "":
        raise PmlstError(
            f"Input error: sample names of input files, {infiles[0]} and {infiles[1]}, \
                   does not share a common sample name. If these files \
                   are paired end reads from the same sample, please rename \
                   them with a common sample name (e.g. 's22_R1.fq', 's22_R2.fq') \
                   or input them seperately."
        )

    return sample_name.rstrip("-").rstrip("_")


def is_gzipped(file_path: str) -> bool:
    """Returns True if file is gzipped and False otherwise.
    The result is inferred from the first two bits in the file read
    from the input path.
    On unix systems this should be: 1f 8b
    Theoretically there could be exceptions to this test but it is
    unlikely and impossible if the input files are otherwise expected
    to be encoded in utf-8.
    """
    with Path(file_path).open("rb") as fh:
        bit_start = fh.read(2)
    return bit_start == b"\x1f\x8b"


def get_file_format(input_files: Sequence[str]) -> str:
    """
    Takes all input files and checks their first character to assess
    the file format. Returns one of the following strings; fasta, fastq,
    other or mixed. fasta and fastq indicates that all input files are
    of the same format, either fasta or fastq. other indiates that all
    files are not fasta nor fastq files. mixed indicates that the inputfiles
    are a mix of different file formats.
    """

    # Open all input files and get the first character
    file_format: list[str] = []
    for infile in input_files:
        if is_gzipped(infile):  # [-3:] == ".gz":
            with gzip.open(infile, "rb") as f:
                fst_char = f.read(1)
        else:
            with Path(infile).open("rb") as f:
                fst_char = f.read(1)
        # Assess the first character
        if fst_char == b"@":
            file_format.append("fastq")
        elif fst_char == b">":
            file_format.append("fasta")
    if len(set(file_format)) != 1:
        return "mixed"
    return ",".join(set(file_format))


def import_profile(database: str, scheme: str, loci_list: list[str]) -> StProfiles:
    """Import all possible allele profiles with corresponding st's
    for the scheme into a dict. The profiles are stored in a dict
    of dicts, to easily look up what st types are accosiated with
    a specific allele number of each loci.
    """
    # Open allele profile file from databaseloci
    with (Path(database) / f"{scheme}.txt.clean").open() as profile_file:
        profile_header = (
            profile_file.readline().strip().split("\t")[1 : len(loci_list) + 1]
        )

        # Create dict for looking up st-types with locus/allele combinations
        st_profiles: StProfiles = {}
        # For each locus initate make an inner dict to store allele and st's
        for locus in loci_list:
            st_profiles[locus] = {}

        # Fill inner dict with allele no as key and seen st-types as value.
        for line in profile_file:
            profile = line.strip().split("\t")
            st_name = profile[0]
            allele_list = profile[1 : len(loci_list) + 1]

            # Save each locus-allele combination with the st-type.
            for i in range(len(allele_list)):
                allele = allele_list[i]
                locus = profile_header[i]
                if allele in st_profiles[locus]:
                    st_profiles[locus][allele] += [st_name]
                else:
                    st_profiles[locus][allele] = [st_name]

    return st_profiles


def st_typing(
    st_profiles: StProfiles, allele_matches: AlleleMatches, loci_list: list[str]
) -> tuple[str, str, str]:
    """
    Takes the path to a dictionary, the inp list of the allele
    number that each loci has been assigned, and an output file string
    where the found st type and similaity is written into it.
    """

    # Find best ST type for all allele profiles
    note = ""

    # First line contains matrix column headers, which are the specific loci
    st_hits: list[str] = []
    st_marks: list[str] = []
    note = ""

    # Check the quality of the alle hits
    for locus in allele_matches:
        allele = allele_matches[locus]["allele"]

        # Check if allele is marked as a non-perfect match. Save mark and write note.
        if "?*" in allele:
            note += f"?* {locus}: Imperfect hit, ST can not be trusted!\n"
            st_marks = ["?", "*"]
        elif "?" in allele:
            note += f"? {locus}: Uncertain hit, ST can not be trusted.\n"
            st_marks.append("?")
        elif "*" in allele:
            note += f"* {locus}: Novel allele, ST may indicate nearest ST.\n"
            st_marks.append("*")

        # Remove mark from allele so it can be used to look up nearest st types
        allele = allele.rstrip("*?!")

        # Get all st's that have the alleles in it's allele profile
        st_hits += st_profiles[locus].get(allele, ["None"])
        if (
            "alternative_hit" in allele_matches[locus]
            and allele_matches[locus]["alternative_hit"] != {}
        ):
            note += f"! {locus}: Multiple perfect hits found\n"
            st_marks.append("!")
            for _allele_name, hit_info in allele_matches[locus][
                "alternative_hit"
            ].items():
                allele = hit_info["allele"].rstrip("!")
                st_hits += st_profiles[locus].get(allele, ["None"])

    # Save allele marks to be transfered to the ST
    st_mark = "".join(set(st_marks))
    notes = st_mark
    # Add marks information to notes
    if "!" in st_mark:
        notes += (
            " alleles with multiple perfect hits found, multiple STs might be found\n"
        )
    if "*" in st_mark and "?" in st_mark:
        notes += " alleles with less than 100% identity and 100% coverages found\n"
    elif st_mark == "*":
        notes = st_mark + " alleles with less than 100% identity found\n"
    elif st_mark == "?":
        notes = st_mark + " alleles with less than 100% coverage found\n"
    notes += note

    # Find most frequent st in st_hits
    st_hits_counter: dict[str, int] = {}
    max_count = 0
    for hit in st_hits:
        if hit != "None":
            if hit in st_hits_counter:
                st_hits_counter[hit] += 1
            else:
                st_hits_counter[hit] = 1
            if max_count < st_hits_counter[hit]:
                max_count = st_hits_counter[hit]

    # Check if allele profile match found st 100 %
    similarity = round(float(max_count) / len(loci_list) * 100, 2)

    if similarity != 100:
        st = "Unknown"
        nearest_st_list: list[str] = []
        # If st is not perfect find nearest st's
        for st_hit, allele_score in sorted(
            st_hits_counter.items(), key=lambda x: x[1], reverse=True
        ):
            if allele_score < max_count:
                break
            nearest_st_list.append(st_hit)
        nearest_sts = ",".join(nearest_st_list)  # + st_mark
    else:
        # The allele profile has a perfect ST hit, but allele marks may indicate
        # imperfect hits.
        sts = [st for st, no in st_hits_counter.items() if no == max_count]
        st = f"{st_mark},".join(sts) + st_mark
        nearest_sts = ""

    return st, notes, nearest_sts


def load_scheme_list_config(scheme_list: SchemeList, database: str) -> SchemeList:
    with (Path(database) / "config").open() as config_file:
        for raw_line in config_file:
            if raw_line.startswith("#"):
                continue
            line = raw_line.split("\t")
            scheme_list[line[0]] = line[1]
    return scheme_list


def plasmidfinder_parsing(
    pf_results: str, scheme_list: SchemeList
) -> tuple[list[str], list[str], list[str]]:
    pf_results_path = Path(pf_results)
    if not pf_results_path.exists():
        raise PmlstError(
            "PlasmidFinder path provided does not exist, "
            "please provide a valid database path."
        )
    with pf_results_path.open() as file:
        data = file.read()

    # Split the header and check if it contains the expected columns
    header = data.split("\n")[0].split("\t")
    expected_columns = ["Database", "Plasmid", "Identity", "Contig"]
    if not all(column in header for column in expected_columns):
        raise PmlstError(
            "PlasmidFinder results file does not contain the expected columns."
        )
    # Check if there are any rows besides header in the file
    if len(data.split("\n")) < 3:
        return ([], [], [])
    list_of_plasmids: list[str] = []
    for raw_line in data.split("\n")[1:]:
        if raw_line:
            line = raw_line.split("\t")
            list_of_plasmids.append(line[1])

    schemes_to_run: list[str] = []
    profile_names_to_run: list[str] = []
    for each_scheme in list(scheme_list.keys()):
        for plasmid in list_of_plasmids:
            if plasmid.lower().startswith(each_scheme):
                if each_scheme not in schemes_to_run:
                    schemes_to_run.append(each_scheme)
                    profile_names_to_run.append(plasmid)
                else:
                    index = schemes_to_run.index(each_scheme)
                    profile_names_to_run[index] += "," + plasmid

    return schemes_to_run, profile_names_to_run, list_of_plasmids


# Keep a local copy because cgecore 2.0.1 contains cgecore/alignment.py, but
# normal imports resolve to cgecore/alignment/__init__.py and extended_cigar is
# not importable from that package path.
def extended_cigar(aligned_template: str, aligned_query: str) -> str:
    """Convert mutation annotations to extended cigar format

    https://github.com/lh3/minimap2#the-cs-optional-tag

    USAGE:
       >>> template = 'CGATCGATAAATAGAGTAG---GAATAGCA'
       >>> query = 'CGATCG---AATAGAGTAGGTCGAATtGCA'
       >>> extended_cigar(template, query) == ':6-ata:10+gtc:4*at:3'
       True
    """
    #   - Go through each position in the alignment
    insertion: list[str] = []
    deletion: list[str] = []
    matches: list[str] = []
    cigar: list[str] = []
    for r_aa, q_aa in zip(
        aligned_template.lower(), aligned_query.lower(), strict=False
    ):
        gap_ref = r_aa == "-"
        gap_que = q_aa == "-"
        match = r_aa == q_aa
        if matches and not match:
            # End match block
            cigar.append(f":{len(matches)}")
            matches = []
        if insertion and not gap_ref:
            # End insertion
            cigar.append(f"+{''.join(insertion)}")
            insertion = []
        elif deletion and not gap_que:
            # End deletion
            cigar.append(f"-{''.join(deletion)}")
            deletion = []
        if gap_ref:
            if insertion:
                # Extend insertion
                insertion.append(q_aa)
            else:
                # Start insertion
                insertion = [q_aa]
        elif gap_que:
            if deletion:
                # Extend deletion
                deletion.append(r_aa)
            else:
                # Start deletion
                deletion = [r_aa]
        elif match:
            if matches:
                # Extend match block
                matches.append(r_aa)
            else:
                # Start match block
                matches = [r_aa]
        else:
            # Add SNP annotation
            cigar.append(f"*{r_aa}{q_aa}")

    if matches:
        cigar.append(f":{len(matches)}")
        del matches
    if insertion:
        # End insertion
        cigar.append(f"+{''.join(insertion)}")
        del insertion
    elif deletion:
        # End deletion
        cigar.append(f"-{''.join(deletion)}")
        del deletion

    return "".join(cigar)


#################################################
