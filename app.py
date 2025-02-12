import os
import subprocess
import time
from flask import Flask, request, send_from_directory, jsonify, render_template
from werkzeug.utils import secure_filename

app = Flask(__name__)

# Configuration
UPLOAD_FOLDER = 'results'
OUTPUT_FOLDER = 'outputs'
ALLOWED_EXTENSIONS = {'pdb', 'mol2'}

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['OUTPUT_FOLDER'] = OUTPUT_FOLDER

os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def add_hydrogens(input_pdb, output_pdb):
    script = f"""open {input_pdb}\naddh\nwrite format pdb #0 {output_pdb}\nquit\n"""
    with open("chimera_script.cmd", "w") as f:
        f.write(script)
    subprocess.run(['/mnt/c/Program Files/Chimera 1.18/bin/chimera.exe', '--nogui', 'chimera_script.cmd'])
    os.remove("chimera_script.cmd")
    return output_pdb

def calculate_charge(input_pdb, charge, output_mol2):
    command = f"antechamber -i {input_pdb} -fi pdb -o {output_mol2} -fo mol2 -c bcc -at gaff2 -nc {charge}"
    subprocess.run(command, shell=True)
    return output_mol2

def generate_topology(mol2_file):
    prmtop = mol2_file.replace('.mol2', '.prmtop')
    inpcrd = mol2_file.replace('.mol2', '.inpcrd')
    final_pdb = mol2_file.replace('.mol2', '_final.pdb')
    gmx_gro = mol2_file.replace('.mol2', 'GMX.gro')
    gmx_top = mol2_file.replace('.mol2', 'GMX.top')
    
    # Placeholder commands, replace with actual ones
    subprocess.run(f'tleap -s -f {mol2_file}', shell=True)
    subprocess.run(f'parmchk2 -i {mol2_file} -f mol2 -o {prmtop}', shell=True)
    
    return prmtop, inpcrd, final_pdb, gmx_gro, gmx_top

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/step1', methods=['POST'])
def step1():
    if 'file' not in request.files:
        return jsonify({"error": "No file uploaded"})
    file = request.files['file']
    charge = request.form.get('charge')
    if not charge:
        return jsonify({"error": "Charge value is required"})
    
    if file and allowed_file(file.filename):
        filename = secure_filename(file.filename)
        timestamp = time.strftime("%Y%m%d%H%M%S")
        input_pdb = os.path.join(UPLOAD_FOLDER, f"{timestamp}_{filename}")
        file.save(input_pdb)
        
        output_pdb = os.path.join(OUTPUT_FOLDER, f"{timestamp}_modified_{filename}")
        add_hydrogens(input_pdb, output_pdb)
        
        output_mol2 = os.path.join(OUTPUT_FOLDER, f"{timestamp}_charge.mol2")
        calculate_charge(output_pdb, charge, output_mol2)
        
        return jsonify({
            "pdb_download": f"/download/{os.path.basename(output_pdb)}",
            "mol2_download": f"/download/{os.path.basename(output_mol2)}"
        })

@app.route('/step2', methods=['POST'])
def step2():
    mol2_file = request.form.get('mol2_file')
    if not mol2_file:
        return jsonify({"error": "MOL2 file is required"})
    
    mol2_path = os.path.join(OUTPUT_FOLDER, mol2_file)
    if not os.path.exists(mol2_path):
        return jsonify({"error": "MOL2 file not found"})
    
    prmtop, inpcrd, final_pdb, gmx_gro, gmx_top = generate_topology(mol2_path)
    
    return jsonify({
        "prmtop_download": f"/download/{os.path.basename(prmtop)}",
        "inpcrd_download": f"/download/{os.path.basename(inpcrd)}",
        "final_pdb_download": f"/download/{os.path.basename(final_pdb)}",
        "gmx_gro_download": f"/download/{os.path.basename(gmx_gro)}",
        "gmx_top_download": f"/download/{os.path.basename(gmx_top)}"
    })

@app.route('/download/<filename>')
def download_file(filename):
    return send_from_directory(OUTPUT_FOLDER, filename, as_attachment=True)

if __name__ == '__main__':
    app.run(debug=True)
