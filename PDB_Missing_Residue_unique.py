from PDB_Missing_Residue import *


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
            output_file = open(args.o, "w")
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
                f"Sequence:{unique_chain_missing_res(pdbseq_with_missing_res(chains, res_names, all_res_names))}"
            )
        if args.res:
            results.append(
                f"Missing residues:{unique_chain_missing_res(missing_res_name(chains, res_names, all_res_names))}"
            )
        if args.num:
            results.append(
                f"Missing residue position:{unique_chain_missing_res(missing_res_num(chains, res_names, res_num))}"
            )
        if args.num_range:
            results.append(
                f"Missing residue position range:{unique_chain_missing_res(convert_to_ranges(missing_res_num(chains, res_names, res_num))[0])}"
            )
        if args.s_num:
            results.append(
                f"Missing residue serial number:{unique_chain_missing_res(missing_res_num(chains, res_names, seqid_res_num))}"
            )
        if args.len:
            results.append(
                f"Missing residue range length:{unique_chain_missing_res(convert_to_ranges(missing_res_num(chains, res_names, res_num))[1])}"
            )

        for result in results:
            if output_file:
                if args.save_as_csv:
                    writer = csv.writer(output_file)
                    if header == 0:
                        filednames = ["ID", "Chain", "residue", "range", "length"]
                        writer.writerow(filednames)
                        header = 1
                    id = [os.path.basename(filename).split('.cif')[0].upper()]
                    res = unique_chain_missing_res(
                        missing_res_name(chains, res_names, all_res_names)
                    )
                    res_position_range = unique_chain_missing_res(
                        convert_to_ranges(missing_res_num(chains, res_names, res_num))[
                            0
                        ]
                    )
                    res_length = unique_chain_missing_res(
                        convert_to_ranges(missing_res_num(chains, res_names, res_num))[
                            1
                        ]
                    )
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


def unique_chain_missing_res(missing_res_dict):
    """Returns a dictionary of [Chains:missing_res_information], where the chains having same missing_res and sequence are combined.
    missing_res_dict= [Chains:missing_res_information] as input
    """
    combined_dict = {}
    for chain, missing_res in missing_res_dict.items():
        if missing_res not in combined_dict.values():
            combined_dict[chain] = missing_res
        else:
            key = [k for k, v in combined_dict.items() if v == missing_res][0]
            combined_dict[key + "," + chain] = combined_dict.pop(key)
    return combined_dict


if __name__ == "__main__":
    main()
