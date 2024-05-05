from PDB_Missing_Residue import pdbseq_with_missing_res, get_input, get_PDBx
import argparse
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import os
import sys


def main():
    args = get_args()

    if args.path:
        path = args.path
    else:
        path = r"temp_files"

    if os.path.exists(args.o):
        sys.exit("Given output file already exists!")
    else:
        output_file = open(args.o, "w", newline="")

    pdb_ids = get_input(args.i, args.all_in_path)
    for filename in get_PDBx(path, pdb_ids, args.all_in_path):
        pdb_info = MMCIF2Dict(filename)
        chains = pdb_info["_pdbx_poly_seq_scheme.asym_id"]
        res_names = pdb_info["_pdbx_poly_seq_scheme.pdb_mon_id"]
        all_res_names = pdb_info["_pdbx_poly_seq_scheme.mon_id"]
        id = filename.split("\\")[-1].split(".cif")[0].upper()
        for chain, sequence in pdbseq_with_missing_res(
            chains, res_names, all_res_names
        ).items():
            output_file.write(">" + id + "_" + chain + "\n")
            output_file.write(sequence + "\n")

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
        help="Enter input filepath containg pdb_ids in each line",
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
        default="output.fasta",
        type=str,
        help="Enter output filepath, if ouput to be saved in a file",
    )

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main()
