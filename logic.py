import subprocess
import os




def process_file(file_path, result_folder):
    """
    Runs the Antechamber command to process the uploaded file and generate the result.

    Args:
        file_path (str): The path to the uploaded file.
        result_folder (str): The folder where the results should be saved.

    Returns:
        str: The path to the resulting file if successful, else None.
    """
    try:
        # Generate output file path
        output_file = os.path.join(result_folder, os.path.splitext(os.path.basename(file_path))[0] + '_am1.mol2')

        # Run the antechamber command
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
        
        return output_file
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")
        return None
