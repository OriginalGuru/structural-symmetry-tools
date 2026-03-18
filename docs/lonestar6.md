# Deploying on TACC Lonestar6

## One-time setup

```bash
# Load the default Python module
module load python3/3.9.7

# Create a persistent virtual environment in your home directory
python3 -m venv $HOME/venvs/vasp-phonon-tools
source $HOME/venvs/vasp-phonon-tools/bin/activate

# Install the only dependency
pip install numpy

# Clone the repository to your work directory
cd $WORK
git clone https://github.com/yourname/vasp-phonon-tools.git
```

The virtual environment persists across login sessions and batch jobs.
You only need to do this once.

## Running interactively

```bash
module load python3/3.9.7
source $HOME/venvs/vasp-phonon-tools/bin/activate

cd /path/to/your/phonon/calculation
python $WORK/vasp-phonon-tools/scripts/get_phonon_mode_charges.py > output.txt
```

## Running in a batch job

```bash
#!/bin/bash
#SBATCH -J phonon_analysis
#SBATCH -A your_allocation
#SBATCH -p vm-small
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH -o phonon_%j.out
#SBATCH -e phonon_%j.err

module load python3/3.9.7
source $HOME/venvs/vasp-phonon-tools/bin/activate

cd $SLURM_SUBMIT_DIR

python $WORK/vasp-phonon-tools/scripts/get_phonon_mode_charges.py > phonon_output.txt
```

Save this as `run_phonon.sh` in your calculation directory and submit with:

```bash
sbatch run_phonon.sh
```

## Notes

- The analysis scripts are single-threaded and lightweight. The `vm-small`
  partition (or `development` for testing) is sufficient.
- `vasprun.xml` can be several hundred MB for large supercells. Ensure the
  job runs from `$SCRATCH` or `$WORK`, not `$HOME` (quota is small).
- The virtual environment must be in `$HOME` (persistent) not `$SCRATCH`
  (purged periodically).
- If you update the repository (`git pull`), no reinstallation is needed
  since there are no compiled components.
