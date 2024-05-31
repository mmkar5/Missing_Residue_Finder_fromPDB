from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1
from Bio.PDB.PDBList import PDBList
import os, glob
import argparse
import csv
import sys


def main():
    args = get_args()

    if args.path:
        path = args.path
    else:
        path = r"temp_files"

    pdb_ids = get_input(args.i, args.all_in_path)

    if args.o:
        if os.path.exists(args.o):
            sys.exit("Given output file already exists!")
        else:
            output_file = open(args.o, "w", newline="")
            if args.save_as_csv:
                header = 0
    else:
        output_file = None

    for filename in get_PDBx(path, pdb_ids, all_in_path=args.all_in_path):
        results = []
        pdb_info = MMCIF2Dict(filename)
        chains = pdb_info["_pdbx_poly_seq_scheme.pdb_strand_id"]
        res_names = pdb_info["_pdbx_poly_seq_scheme.pdb_mon_id"]
        all_res_names = pdb_info["_pdbx_poly_seq_scheme.mon_id"]
        res_num = pdb_info["_pdbx_poly_seq_scheme.pdb_seq_num"]
        seqid_res_num = pdb_info["_pdbx_poly_seq_scheme.seq_id"]

        results.append(f"PDB id:{os.path.basename(filename).split('.cif')[0].upper()}")

        if args.seq:
            results.append(
                f"Sequence:{pdbseq_with_missing_res(chains, res_names, all_res_names)}"
            )
        if args.res:
            results.append(
                f"Missing residues:{missing_res_name(chains, res_names, all_res_names)}"
            )
        if args.num:
            results.append(
                f"Missing residue position:{missing_res_num(chains, res_names, res_num, all_res_names)}"
            )
        if args.num_range:
            results.append(
                f"Missing residue position range:{convert_to_ranges(missing_res_num(chains, res_names, res_num, all_res_names))[0]}"
            )
        if args.s_num:
            results.append(
                f"Missing residue serial number:{missing_res_num(chains, res_names, seqid_res_num, all_res_names)}"
            )
        if args.len:
            results.append(
                f"Missing residue range length:{convert_to_ranges(missing_res_num(chains, res_names, res_num, all_res_names))[1]}"
            )

        for result in results:
            if output_file:
                if args.save_as_csv:
                    writer = csv.writer(output_file)
                    if header == 0:
                        filednames = ["ID", "Chain", "residue", "range", "length"]
                        writer.writerow(filednames)
                        header = 1
                    id = [os.path.basename(filename).split(".cif")[0].upper()]
                    res = missing_res_name(chains, res_names, all_res_names)
                    res_position_range = convert_to_ranges(
                        missing_res_num(chains, res_names, res_num, all_res_names)
                    )[0]
                    res_length = convert_to_ranges(
                        missing_res_num(chains, res_names, res_num, all_res_names)
                    )[1]
                    writer.writerows(
                        zip(
                            id * len(res.keys()),
                            res.keys(),
                            res.values(),
                            res_position_range.values(),
                            res_length.values(),
                        )
                    )

                else:
                    output_file.write(result + "\n")
            else:
                print(result)

    if output_file:
        output_file.close()


def get_args():
    """Parse command line arguments given by the user for obtaining missing residue information from protein mmCIF files (save as corresponding PDB file mssing residue information)"""

    parser = argparse.ArgumentParser(
        description="Arguments to determine missing residues in mmCIF file"
    )
    parser.add_argument(
        "-i",
        required=False,
        type=str,
        help="Enter input filepath containg pdb_ids in each line or comma seperated",
    )
    parser.add_argument(
        "-all_in_path",
        action="store_true",
        help="uses all the mmCIF files in the given path as input, ignores -pdb_id_list",
    )
    parser.add_argument(
        "-path",
        required=False,
        type=str,
        help="Enter the path where the mmCIF input files are present, or will be dowloaded in",
    )

    parser.add_argument(
        "-o",
        required=False,
        type=str,
        help="Enter output filepath, if output to be saved in a file",
    )
    parser.add_argument(
        "-seq",
        action="store_true",
        help="obtain sequence of pdb coordinates in single line with missing residues as small letters",
    )
    parser.add_argument("-res", action="store_true", help="obtain the missing residues")
    parser.add_argument(
        "-num",
        action="store_true",
        help="obtain the positions of the missing residues in the pdb(mmCIF) file",
    )
    parser.add_argument(
        "-num_range",
        action="store_true",
        help="obtain the range of missing residue position in the pdb(mmCIF) file",
    )
    parser.add_argument(
        "-s_num",
        action="store_true",
        help="obtain the missing residue position assuming the first residue position is one",
    )
    parser.add_argument(
        "-len",
        action="store_true",
        help="obtain the length of ranges of missing residue position in the pdb(mmCIF) file",
    )
    parser.add_argument(
        "-save_as_csv",
        action="store_true",
        help="saves the output in the output filepath provided as csv [ID,Chain,residue,range,length] ,only works if output filepath provided;ignores aruments other than -i or -all_in_path, -o, -path",
    )

    args = parser.parse_args()
    return args


def get_input(filepath=None, all_in_path=False):
    """Returns the pdb_ids as a list. Returns the pdb_ids as a list. Takes input from user as a input file containing pdb_ids in seperate lines or comma-seperated,
    or prompts user to provide the pdb id in the terminal.
    filepath= name of the input file(with filepath, if not in current working directory) (optional)
    all_in_path= True if all the mmCIF files in the path needs to be processed (Default=False)
    """
    pdb_ids = []
    if filepath:
        with open(filepath) as file:
            for line in file:
                if "," in line:
                    pdb_ids = line.split(",")
                else:
                    if not line.isspace():
                        pdb_ids.append(line.strip().lower())
    else:
        print("Press Ctrl+C, when no input left to enter")
        while True:
            try:
                pdb_ids.append(input("Enter pdb_id:").lower())
            except KeyboardInterrupt:
                print()
                break
    return pdb_ids


def get_PDBx(path, pdb_ids=None, all_in_path=False):
    """Locates or downloads the required mmCIF files in the given path and returns the filepath
    path = The location of the mmCIF files
    pdb_ids = List of the pdb_ids (optional)
    all_in_path = True if all the mmCIF files in the path needs to be processed.(Default=False)
    """

    if pdb_ids:
        pdbl = PDBList()
        pdbl.download_pdb_files(pdb_ids, file_format="mmCif", pdir=path)
    if all_in_path:
        for filename in glob.glob(os.path.join(path, "*.cif")):
            yield filename
    else:
        for filename in glob.glob(os.path.join(path, "*.cif")):
            if os.path.basename(filename).split(".cif")[0].lower() in pdb_ids:
                yield filename


def pdbseq_with_missing_res(chains, res_names, all_res_names):
    """Returns a dictionary with the protein chain as keys and the protein sequence as values.
    The protein sequence has one letter code, in uppercase, except for the missing residues which is in lowercase.
    Coverts three letter amino acid codes into one letter code before processing.
    chains = List of protein chains
    res_names = List of all amino acid residues present in protein structure, where missing residues indicated as ?
    all_res_names = List of all amino acid residues present in protein structure
    """
    output = {}
    res_names = [res if res == "?" else seq1(res,undef_code='') for res in res_names]
    missing_res = [seq1(res,undef_code='').lower() for res in all_res_names]
    for i in range(len(chains)):
        if i == 0:
            if res_names[i] == "?":
                var = missing_res[i]
            else:
                var = res_names[i]
        else:
            if chains[i] == chains[i - 1]:
                if res_names[i] == "?":
                    var += missing_res[i]
                else:
                    var += res_names[i]
                output[chains[i]] = var
            else:
                var = ""
                if res_names[i] == "?":
                    var += missing_res[i]
                else:
                    var += res_names[i]
                output[chains[i]] = var
    output={k:v for k,v in output.items() if not v==""}
    return output


def missing_res_num(chains, res_names, res_num, all_res_names):
    """Returns a dictionary with the protein chain as keys and the position of missing amino acid residues as values.
    Coverts three letter amino acid codes into one letter code before processing.
    chains = List of protein chains
    res_names = List of residues present in protein structure, where missing residues indicated as ?
    res_num = List of numerical positions of the amino acid residues
    all_res_names = List of all amino acid residues present in protein structure
    """
    output = {}
    res_names = [res if res == "?" else seq1(res,undef_code='') for res in res_names]
    var = ""
    for i in range(len(chains)):
        if i == 0:
            if res_names[i] == "?" and seq1(all_res_names[i]) in list("ARNDCQEGHILKMFPSTWYVBJXUO"):
                var = res_num[i] + ","
        else:
            if chains[i] == chains[i - 1]:
                if res_names[i] == "?" and seq1(all_res_names[i]) in list("ARNDCQEGHILKMFPSTWYVBJXUO"):
                    var += res_num[i] + ","
                output[chains[i]] = var.strip(",")
            else:
                var = ""
                if res_names[i] == "?" and seq1(all_res_names[i]) in list("ARNDCQEGHILKMFPSTWYVBJXUO"):
                    var += res_num[i] + ","
                output[chains[i]] = var.strip(",")
    output={k:v for k,v in output.items() if not v==""}
    return output


def missing_res_name(chains, res_names, all_res_names):
    """Returns a dictionary with the protein chain as keys and the missing amino acid residues as values.
    Coverts three letter amino acid codes into one letter code before processing.
    chains = List of protein chains
    res_names = List of all amino acid residues present in protein structure, where missing residues indicated as ?
    all_res_names = List of all amino acid residues present in protein structure
    """
    output = {}
    res_names = [res if res == "?" else seq1(res) for res in res_names]
    missing_res = [seq1(res,undef_code="").lower() for res in all_res_names]
    for i in range(len(chains)):
        if i == 0:
            if res_names[i] == "?":
                var = missing_res[i]
            else:
                var = ""
        else:
            if chains[i] == chains[i - 1]:
                if res_names[i] == "?":
                    var += missing_res[i]
                else:
                    if not var.endswith(" "):
                        var += " "
                output[chains[i]] = var.strip()
            else:
                var = ""
                if res_names[i] == "?":
                    var += missing_res[i]
                else:
                    if not var.endswith(" "):
                        var += " "
                output[chains[i]] = var.strip()
    output={k:v for k,v in output.items() if not v==""}
    return output


def convert_to_ranges(num_dict):
    """Returns a tuple containing {'protein chain':'range of amino acid position'},{'protein chain':'length of the range of amino acid positions'}.
    num_dict = Dictionary containing keys as chains and values as position of amino acid residues
    """
    range_dict = {}
    len_dict = {}
    for key, value in num_dict.items():
        try:
            numbers = list(map(int, value.split(",")))
        except ValueError:
            range_str = value
            range_dict[key] = range_str
            len_dict[key] = ""
        else:
            numbers.sort()
            ranges = []
            start = end = numbers[0]
            for num in numbers[1:]:
                if num != end + 1:
                    ranges.append((start, end))
                    start = end = num
                else:
                    end = num
            ranges.append((start, end))

            range_str = ",".join(f"{s}-{e}" if s != e else str(s) for s, e in ranges)
            range_dict[key] = range_str

            len_str = ",".join(str(e - s + 1) for s, e in ranges)
            len_dict[key] = len_str

    return range_dict, len_dict


if __name__ == "__main__":
    main()
