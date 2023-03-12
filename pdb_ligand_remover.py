from Bio.PDB import PDBParser, PDBIO, Select, PDBList

class ProteinOnly(Select):
    def accept_residue(self, residue):
        hetero_flag = residue.get_id()[0].strip()
        if hetero_flag in ["W", "H", "L"]:
            return False
        if hetero_flag == "M":
            return False
        if residue.get_resname() not in ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]:
            return False
        return True

def remove_ligands(pdb_code, output_file):
#      Download PDB file
    pdbl = PDBList()
    pdb_file = pdbl.retrieve_pdb_file(pdb_code, file_format='pdb')

    # Parse the input PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_code, pdb_file)

    # Remove hetero atoms and non-amino acid residues
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file, select=ProteinOnly())
    print(f"Ligands removed from {pdb_code}. Output saved as {output_file}.")

# Example usage
remove_ligands("7pom", "7pom.pdb")