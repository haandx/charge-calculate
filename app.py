from flask import Flask, request, render_template, send_from_directory
import os
from logic import process_file  # Import the function from logic.py

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

    # Process the file using logic.py's process_file function
    output_file = process_file(file_path, RESULT_FOLDER)

    if output_file:
        # Send the processed file back to the user
        return send_from_directory(RESULT_FOLDER, os.path.basename(output_file), as_attachment=True)
    else:
        return 'Error occurred during processing.', 500


if __name__ == '__main__':
    app.run(debug=True)
