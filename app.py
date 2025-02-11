import os
import subprocess
from pathlib import Path
from flask import Flask, request, send_from_directory, jsonify
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

@app.route('/')
def index():
    return '''
    <h1>Upload PDB file</h1>
    <form method="POST" enctype="multipart/form-data" action="/add_hydrogens">
        <input type="file" name="file">
        <input type="submit" value="Upload">
    </form>
    '''

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
        
        # Call the function to add hydrogens
        add_hydrogens(input_pdb, output_pdb)
        
        # Provide the download link for the processed file
        return send_from_directory(app.config['OUTPUT_FOLDER'], output_filename, as_attachment=True)

    return jsonify({"error": "Invalid file format. Please upload a .pdb file."})

if __name__ == '__main__':
    app.run(debug=True)
