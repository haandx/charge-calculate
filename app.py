import os
import subprocess
import time
from pathlib import Path
from flask import Flask, request, send_from_directory, jsonify, render_template
from werkzeug.utils import secure_filename

app = Flask(__name__)

# Configuration
UPLOAD_FOLDER = 'results'
OUTPUT_FOLDER = 'outputs'
ALLOWED_EXTENSIONS = {'pdb', 'mol2'}

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

def optimize_geometry(input_mol2, output_mol2):
    """Optimize geometry using Antechamber"""
    command = f"antechamber -i {input_mol2} -fi mol2 -o {output_mol2} -fo mol2 -c bcc -at gaff2"
    subprocess.run(command, shell=True)
    return output_mol2

def generate_frcmod(input_mol2, output_frcmod):
    """Generate force field parameters using parmchk2"""
    command = f"parmchk2 -i {input_mol2} -f mol2 -o {output_frcmod}"
    subprocess.run(command, shell=True)
    return output_frcmod

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
        timestamp = time.strftime("%Y%m%d%H%M%S")
        
        # Modify filenames to include a timestamp
        input_pdb = os.path.join(app.config['UPLOAD_FOLDER'], f"{timestamp}_{filename}")
        file.save(input_pdb)

        # Add hydrogens
        output_filename = f"{timestamp}_modified_{filename}"
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

        # Optimize geometry
        optimized_mol2 = output_mol2.replace("_charge.mol2", "_optimized.mol2")
        optimized_mol2_path = os.path.join(app.config['OUTPUT_FOLDER'], optimized_mol2)
        optimize_geometry(output_mol2_path, optimized_mol2_path)

        # Generate force field parameters
        frcmod_file = optimized_mol2.replace(".mol2", ".frcmod")
        frcmod_file_path = os.path.join(app.config['OUTPUT_FOLDER'], frcmod_file)
        generate_frcmod(optimized_mol2_path, frcmod_file_path)

        # Provide the download links for processed files
        return jsonify({
            "mol2_download": f"/download/{optimized_mol2}",
            "frcmod_download": f"/download/{frcmod_file}"
        })

    return jsonify({"error": "Invalid file format. Please upload a .pdb or .mol2 file."})

@app.route('/download/<filename>')
def download_file(filename):
    return send_from_directory(app.config['OUTPUT_FOLDER'], filename, as_attachment=True)

if __name__ == '__main__':
    app.run(debug=True)
