import os
import subprocess
from pathlib import Path

def add_hydrogens(input_pdb, output_pdb):
    """Add hydrogens using Chimera"""
    script = f"""open {input_pdb}
addh
write format pdb #0 {output_pdb}
quit
"""
    with open("chimera_script.cmd", "w") as f:
        f.write(script)
    
    subprocess.run(['chimera', '--nogui', 'chimera_script.cmd'])
    os.remove("chimera_script.cmd")

    return output_pdb  