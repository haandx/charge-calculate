from flask import Flask, request, render_template, send_from_directory, jsonify
import os
import time
from rdkit import Chem
from rdkit.Chem import AllChem

app = Flask(__name__)

UPLOAD_FOLDER = "uploads"
OUTPUT_FOLDER = "outputs"
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['OUTPUT_FOLDER'] = OUTPUT_FOLDER

def add_hydrogens_and_calculate_charges(input_pdb, output_mol2):
    """
    Add hydrogens, calculate charges, and write charges to the MOL2 file in the required format.
    """
    try:
        # Load the molecule from the PDB file
        mol = Chem.MolFromPDBFile(input_pdb, removeHs=False)
        if mol is None:
            raise ValueError("Failed to load molecule from PDB file.")

        # Add hydrogens
        mol = Chem.AddHs(mol, addCoords=True)

        # Calculate Gasteiger charges
        AllChem.ComputeGasteigerCharges(mol)

        # Write the molecule to a new MOL2 file with charges in the correct format
        with open(output_mol2, 'w') as f:
            # Write molecule header
            f.write("@<TRIPOS>MOLECULE\n")
            f.write("4WI\n")  # Placeholder molecule name (can be replaced with actual name)
            f.write(f"   {mol.GetNumAtoms()}    {mol.GetNumBonds()}     1     0     0\n")
            f.write("SMALL\n")
            f.write("bcc\n\n")

            # Write atom section header
            f.write("@<TRIPOS>ATOM\n")

            # Write atom data, including charges
            for atom in mol.GetAtoms():
                atom_idx = atom.GetIdx() + 1  # MOL2 atom indices start at 1
                atom_name = atom.GetSymbol()  # Atom name
                atom_pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
                charge = atom.GetDoubleProp("_GasteigerCharge")
                f.write(f"{atom_idx:7} {atom_name:<4}    {atom_pos.x: .4f}    {atom_pos.y: .4f}   {atom_pos.z: .4f} "
                        f"{atom.GetSymbol().lower():<2}   401 4WI       {charge: .6f}\n")

        print(f"File saved successfully: {output_mol2}")
        return output_mol2
    except Exception as e:
        print(f"Error in add_hydrogens_and_calculate_charges: {e}")
        return None

@app.route("/", methods=["GET", "POST"])
def upload_file():
    if request.method == "POST":
        file = request.files["file"]
        if file:
            timestamp = int(time.time())
            base_filename = os.path.splitext(file.filename)[0]
            input_path = os.path.join(UPLOAD_FOLDER, f"{base_filename}_{timestamp}.pdb")
            output_mol2 = os.path.join(OUTPUT_FOLDER, f"{base_filename}_H_charged_{timestamp}.mol2")

            # Save the uploaded file
            file.save(input_path)

            # Add hydrogens and calculate charges
            try:
                add_hydrogens_and_calculate_charges(input_path, output_mol2)
                return render_template("index.html", mol2_file=os.path.basename(output_mol2))
            except Exception as e:
                return render_template("index.html", error=str(e))

    return render_template("index.html", mol2_file=None)

@app.route('/outputs/<filename>')
def output_file(filename):
    return send_from_directory(OUTPUT_FOLDER, filename)

@app.route('/save', methods=['POST'])
def save_file():
    data = request.json
    filename = os.path.join(OUTPUT_FOLDER, data['filename'])

    if os.path.exists(filename):  
        with open(filename, 'w') as file:
            file.write(data['content'])
        return jsonify({"message": "File saved successfully!"})
    else:
        return jsonify({"error": "File not found"}), 404

@app.route('/download/<filename>')
def download_file(filename):
    return send_from_directory(OUTPUT_FOLDER, filename, as_attachment=True)

if __name__ == "__main__":
    app.run(debug=True)
