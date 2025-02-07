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

def add_hydrogens_and_calculate_charges(input_pdb, output_mol2, user_charge):
    """
    Add hydrogens, set user-defined charge, and calculate Gasteiger charges.
    """
    try:
        mol = Chem.MolFromPDBFile(input_pdb, removeHs=False)
        if mol is None:
            raise ValueError("Failed to load molecule from PDB file.")

        mol = Chem.AddHs(mol, addCoords=True)
        
        # Assign user-defined charge
        for atom in mol.GetAtoms():
            atom.SetFormalCharge(user_charge)
        
        AllChem.ComputeGasteigerCharges(mol)

        with open(output_mol2, 'w') as f:
            f.write("@<TRIPOS>MOLECULE\n")
            f.write("Molecule\n")
            f.write(f"   {mol.GetNumAtoms()}    {mol.GetNumBonds()}     1     0     0\n")
            f.write("SMALL\n")
            f.write("bcc\n\n")
            f.write("@<TRIPOS>ATOM\n")
            
            for atom in mol.GetAtoms():
                atom_idx = atom.GetIdx() + 1
                atom_name = atom.GetSymbol()
                atom_pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
                charge = atom.GetDoubleProp("_GasteigerCharge")
                f.write(f"{atom_idx:7} {atom_name:<4} {atom_pos.x:.4f} {atom_pos.y:.4f} {atom_pos.z:.4f} {atom.GetSymbol().lower()} 401 MOL {charge:.6f}\n")

        return output_mol2
    except Exception as e:
        print(f"Error: {e}")
        return None

@app.route("/", methods=["GET", "POST"])
def upload_file():
    if request.method == "POST":
        file = request.files["file"]
        user_charge = int(request.form.get("charge", 0))  # Get charge input, default is 0

        if file:
            timestamp = int(time.time())
            base_filename = os.path.splitext(file.filename)[0]
            input_path = os.path.join(UPLOAD_FOLDER, f"{base_filename}_{timestamp}.pdb")
            output_mol2 = os.path.join(OUTPUT_FOLDER, f"{base_filename}_H_charged_{timestamp}.mol2")

            file.save(input_path)
            
            try:
                add_hydrogens_and_calculate_charges(input_path, output_mol2, user_charge)
                return render_template("index.html", mol2_file=os.path.basename(output_mol2))
            except Exception as e:
                return render_template("index.html", error=str(e))

    return render_template("index.html", mol2_file=None)

@app.route('/outputs/<filename>')
def output_file(filename):
    return send_from_directory(OUTPUT_FOLDER, filename)

@app.route('/download/<filename>')
def download_file(filename):
    return send_from_directory(OUTPUT_FOLDER, filename, as_attachment=True)

if __name__ == "__main__":
    app.run(debug=True)
