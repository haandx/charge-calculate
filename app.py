import os
import subprocess
import time
import tempfile
from flask import Flask, request, jsonify, send_from_directory, render_template
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

def check_parameters(mol2_file, frcmod_file):
    """Check if all parameters are available using parmchk2"""
    cmd = ['parmchk2', '-i', mol2_file, '-f', 'mol2', '-o', frcmod_file]
    subprocess.run(cmd)
    print("Force field parameters checked ✓")

def generate_amber_params(mol2_file, frcmod_file, prmtop_file, inpcrd_file, pdb_output):
    """Generate topology and coordinate files using TLeap"""
    tleap_input = "pfr.tleap"
    
    with open(tleap_input, "w") as f:
        f.write(f"""
source leaprc.gaff2
loadamberparams {frcmod_file}
PFR = loadmol2 {mol2_file}
saveamberparm PFR {prmtop_file} {inpcrd_file}
savepdb PFR {pdb_output}
quit
""")
    
    subprocess.run(['tleap', '-f', tleap_input])
    os.remove(tleap_input)
    print("Amber parameter files generated ✓")

def convert_to_gromacs(prmtop_file, inpcrd_file):
    """Convert Amber topology to GROMACS format using ACPYPE"""
    cmd = ['acpype', '-p', prmtop_file, '-x', inpcrd_file, '-b', "4wi_gro"]
    subprocess.run(cmd)
    print("GROMACS topology files generated ✓")

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
        return jsonify({"error": "No file uploaded"}), 400
    file = request.files['file']
    charge = request.form.get('charge')
    if not charge:
        return jsonify({"error": "Charge value is required"}), 400
    
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
        return jsonify({"error": "MOL2 file is required"}), 400
    
    mol2_path = os.path.join(OUTPUT_FOLDER, mol2_file)
    if not os.path.exists(mol2_path):
        return jsonify({"error": "MOL2 file not found"}), 404
    
    # Generate Amber Parameters
    frcmod_file = os.path.join(tempfile.gettempdir(), 'output.frcmod')
    prmtop_file = os.path.join(tempfile.gettempdir(), 'output.prmtop')
    inpcrd_file = os.path.join(tempfile.gettempdir(), 'output.inpcrd')
    pdb_output = os.path.join(tempfile.gettempdir(), 'output.pdb')

    try:
        check_parameters(mol2_path, frcmod_file)
        generate_amber_params(mol2_path, frcmod_file, prmtop_file, inpcrd_file, pdb_output)
        convert_to_gromacs(prmtop_file, inpcrd_file)

        return jsonify({
            'message': 'Files processed successfully!',
            'files': {
                'frcmod': frcmod_file,
                'prmtop': prmtop_file,
                'inpcrd': inpcrd_file,
                'pdb': pdb_output
            }
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/download/<filename>')
def download_file(filename):
    return send_from_directory(OUTPUT_FOLDER, filename, as_attachment=True)

if __name__ == '__main__':
    app.run(debug=True)
