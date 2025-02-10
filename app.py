from flask import Flask, request, render_template, send_from_directory, jsonify
import os
import time
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem

app = Flask(__name__)

UPLOAD_FOLDER = "uploads"
OUTPUT_FOLDER = "outputs"
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['OUTPUT_FOLDER'] = OUTPUT_FOLDER

def convert_pdb_to_mol2(input_pdb, output_mol2):
    """Convert PDB file to MOL2 format."""
    mol = Chem.MolFromPDBFile(input_pdb, removeHs=False)
    if mol is None:
        raise ValueError("Error: Could not parse PDB file.")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    
    # Generate MOL2 content manually
    mol2_content = Chem.MolToMolBlock(mol)
    with open(output_mol2, 'w') as f:
        f.write(mol2_content)

def create_mopac_input(mol2_file, mopac_input_file):
    """Create MOPAC input file for AM1 charge calculation."""
    with open(mopac_input_file, 'w') as f:
        f.write(f"PM7\n")  # Using PM7, which is the default for MOPAC
        f.write(f"Molecule input\n")
        f.write(f"AM1\n")  # Using AM1 method for charge calculation

def run_mopac(mopac_input_file):
    """Run MOPAC using subprocess."""
    mopac_executable = r"C:\Program Files\MOPAC\bin\mopac.exe"  # Full path to mopac.exe
    result = subprocess.run([mopac_executable, mopac_input_file], capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"MOPAC run failed: {result.stderr}")

def extract_charges(mopac_output_file):
    """Extract AM1 charges from MOPAC output file."""
    charges = []
    with open(mopac_output_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if "FINAL CHARGES" in line:  # Adjust this pattern based on MOPAC output
                charges.append(line.strip())
    return charges


@app.route("/", methods=["GET", "POST"])
def upload_file():
    if request.method == "POST":
        file = request.files["file"]
        user_charge = int(request.form.get("charge", 0))  # Get charge input, default is 0

        if file:  # Creates unique filenames based on the current timestamp to prevent overwriting
            timestamp = int(time.time())
            base_filename = os.path.splitext(file.filename)[0]
            input_path = os.path.join(UPLOAD_FOLDER, f"{base_filename}_{timestamp}.pdb")
            output_mol2 = os.path.join(OUTPUT_FOLDER, f"{base_filename}_{timestamp}.mol2")
            mopac_input = os.path.join(OUTPUT_FOLDER, f"{base_filename}_{timestamp}.mop")
            mopac_output = os.path.join(OUTPUT_FOLDER, f"{base_filename}_{timestamp}.out")

            file.save(input_path)

            try:
                # Convert PDB to MOL2
                convert_pdb_to_mol2(input_path, output_mol2)

                # Create MOPAC input file for AM1 charge calculation
                create_mopac_input(output_mol2, mopac_input)

                # Run MOPAC to calculate charges
                run_mopac(mopac_input)

                # Extract AM1 charges from MOPAC output
                charges = extract_charges(mopac_output)

                # Render the result
                return render_template("index.html", charges=charges, mol2_file=os.path.basename(output_mol2))

            except Exception as e:
                return render_template("index.html", error=str(e))

    return render_template("index.html", charges=None)

@app.route('/outputs/<filename>')
def output_file(filename):
    """Serve the output MOL2 or MOPAC file."""
    return send_from_directory(OUTPUT_FOLDER, filename)

@app.route('/download/<filename>')
def download_file(filename):
    """Download the output file."""
    return send_from_directory(OUTPUT_FOLDER, filename, as_attachment=True)

if __name__ == "__main__":
    app.run(debug=True)
