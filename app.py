import os
import subprocess
from pathlib import Path
from flask import Flask, request, send_from_directory, jsonify, render_template
from werkzeug.utils import secure_filename

app = Flask(__name__)

# Configuration
UPLOAD_FOLDER = 'results'
OUTPUT_FOLDER = 'outputs'
ALLOWED_EXTENSIONS = {'pdb'}

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['OUTPUT_FOLDER'] = OUTPUT_FOLDER

# Make sure the directories exist
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def add_hydrogens(input_pdb, output_pdb):
    """Add hydrogens using Chimera"""
    script = f"""open {input_pdb}
addh
write format pdb #0 {output_pdb}
quit
"""
    with open("chimera_script.cmd", "w") as f:
        f.write(script)
    
    subprocess.run(['/mnt/c/Program Files/Chimera 1.18/bin/chimera.exe', '--nogui', 'chimera_script.cmd'])
    os.remove("chimera_script.cmd")

    return output_pdb

def calculate_charge(input_pdb, charge, output_mol2):
    """Calculate charge using Antechamber"""
    command = f"antechamber -i {input_pdb} -fi pdb -o {output_mol2} -fo mol2 -c bcc -at gaff2 -nc {charge}"
    subprocess.run(command, shell=True)
    return output_mol2

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/add_hydrogens', methods=['POST'])
def upload_file():
    if 'file' not in request.files:
        return jsonify({"error": "No file part"})
    file = request.files['file']
    
    if file.filename == '':
        return jsonify({"error": "No selected file"})
    
    if file and allowed_file(file.filename):
        filename = secure_filename(file.filename)
        input_pdb = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        file.save(input_pdb)

        # Add hydrogens
        output_filename = f"modified_{filename}"
        output_pdb = os.path.join(app.config['OUTPUT_FOLDER'], output_filename)
        add_hydrogens(input_pdb, output_pdb)
        
        # Get the charge input from the user
        charge = request.form.get('charge')
        if charge is None:
            return jsonify({"error": "Charge value is required."})

        # Calculate charge using Antechamber
        output_mol2 = output_filename.replace(".pdb", "_charge.mol2")
        output_mol2_path = os.path.join(app.config['OUTPUT_FOLDER'], output_mol2)
        calculate_charge(output_pdb, charge, output_mol2_path)

        # Provide the download link for the processed file
        return send_from_directory(app.config['OUTPUT_FOLDER'], output_mol2, as_attachment=True)

    return jsonify({"error": "Invalid file format. Please upload a .pdb file."})

if __name__ == '__main__':
    app.run(debug=True)
