from flask import Flask, request, render_template, send_from_directory
import os
import subprocess

app = Flask(__name__)

# Set the path where the user uploads files and stores results
UPLOAD_FOLDER = 'uploads'
RESULT_FOLDER = 'results'

# Make sure the folders exist
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(RESULT_FOLDER, exist_ok=True)

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

# Route to display the upload form
@app.route('/')
def index():
    return render_template('index.html')

# Route to handle the file upload and processing
@app.route('/upload', methods=['POST'])
def upload_file():
    if 'file' not in request.files:
        return 'No file part', 400
    file = request.files['file']
    
    if file.filename == '':
        return 'No selected file', 400
    
    # Save the uploaded file
    file_path = os.path.join(app.config['UPLOAD_FOLDER'], file.filename)
    file.save(file_path)

    # Call Antechamber to add hydrogen and generate mol2
    try:
        output_file = os.path.join(RESULT_FOLDER, file.filename.split('.')[0] + '_am1.mol2')
        command = [
            'antechamber',
            '-i', file_path,
            '-fi', 'pdb',
            '-o', output_file,
            '-fo', 'mol2',
            '-c', 'bcc',
            '-at', 'gaff2',
            '-nc', '0'
        ]
        
        subprocess.run(command, check=True)
        
        # Send the processed file back to the user
        return send_from_directory(RESULT_FOLDER, os.path.basename(output_file), as_attachment=True)

    except subprocess.CalledProcessError as e:
        return f"Error occurred: {e}", 500


if __name__ == '__main__':
    app.run(debug=True)
