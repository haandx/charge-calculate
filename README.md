Summary
User uploads a PDB file with an optional charge value.
The file is saved with a unique timestamp.
The RDKit library processes the molecule:
Loads the structure.
Adds hydrogens.
Assigns a user-defined charge.
Computes Gasteiger charges.
Saves the output in MOL2 format.
The user receives a download link for the processed file.
This Flask app provides an easy-to-use web interface for modifying molecular files, making it useful for computational chemistry and bioinformatics applications. 🚀
