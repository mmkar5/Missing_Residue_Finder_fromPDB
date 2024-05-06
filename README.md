# Missing Residue Finder: Finds the missing amino acid residue information from a given PDB id
## Contents
1. **PDB_Missing_Residue.py**
2. **Missing_Residue_Sequence.py**
3. **PDB_Missing_Residue_unique.py**

The program requires biophython to be installed.

## Description for PDB_Missing_Residue
The program **PDB_Missing_Residue.py** takes as input PDB_IDs or looks for all the mmCIF file present in the given directory, and for each such file, can provide the following information -
- Protein sequence with missing residues in lower case
- Position of the missing residues
- One letter code of the missing residues
- Range of the position of missing residues and their length

`usage: PDB_Missing_Residue.py [-h] [-i I] [-all_in_path] [-path PATH] [-o O] [-seq] [-res] [-num] [-num_range] [-s_num] [-len] [-save_as_csv]`

The program takes the following command line arguments to determine missing residues:
```
  -h, --help    show this help message and exit
  -i I          Enter input filepath containg pdb_ids in each line
  -all_in_path  uses all the mmCIF files in the given path as input, ignores -i
  -path PATH    Enter the path where the mmCIF input files are present, or will be dowloaded in
  -o O          Enter output filepath, if output to be saved in a file
  -seq          obtain sequence of pdb coordinates in single line with missing residues as small letters
  -res          obtain the missing residues
  -num          obtain the positions of the missing residues in the pdb(mmCIF) file
  -num_range    obtain the range of missing residue position in the pdb(mmCIF) file
  -s_num        obtain the missing residue position assuming the first residue position is one
  -len          obtain the length of ranges of missing residue position in the pdb(mmCIF) file
  -save_as_csv  saves the output in the output filepath provided as csv [ID,Chain,residue,range,length] ,only works if output filepath provided;ignores aruments other than -i or
                -all_in_path, -o, -path
```
All of the arguments are optional.

The following ways of providing inputs are accepted:
- An input file containing pdb ids in each line
- Path of a diectory where the mmCIF files are present
- User can type the pdb id, when prompted to do so

The outputs are provided in the following ways:
- If output filepath is provided, output will be saved in the file
- The output can be printed at the command line interface
- The output can be saved in a csv format - id,chain,residue,range,length
 
 If the given output file already exists, then the program will exit, indicating "Given output file already exists!", preventing any overwriting on an existing file. By default all files will be downloaded to "temp_files" folder in current working directory.

## Description for Missing_Residue_Sequence

The program **Missing_Residue_Sequence.py** takes the same kind of inputs as **Missing_Residue_Sequence.py**. It outputs the protein sequence, in the given file in a fasta format, where the missing residues in the protein sequence are in lower case, while the rest is in upper case. First, the identifier is given as >pdb_chain, then the next line contains the sequence and so on. If -unique is used, then chains having same sequence as well as missing residue information are removed and such chain names are concatenated.

If no ouput filepath (or filename) is provided, then uses default "output.fasta" as output file. If the output file already exists, then the program will exit, indicating "Given output file already exists!", preventing any overwriting on an existing file. By default all files will be downloaded to "temp_files" folder in current working directory.


`usage: Missing_Residue_Sequence.py [-h] [-i I] [-all_in_path] [-path PATH] [-o O] [-unique] `

The program takes the following command line arguments
```
  -h, --help    show this help message and exit
  -i I          Enter input filepath containg pdb_ids in each line
  -all_in_path  uses all the mmCIF files in the given path as input, ignores -pdb_id_list
  -path PATH    Enter the path where the mmCIF input files are present, or will be dowloaded in
  -o O          Enter output filepath, if output to be saved in a file
  -unique       Remove chains having same sequence and mssing residue
```
## Description for Missing_Residue_Sequence_unique

Same as **Missing_Residue_Sequence**, except that the chains having same sequence as well as missing residue information are removed and such chain names are concatenated.



Cock, P.J. et al., 2009. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), pp.1422–1423.
