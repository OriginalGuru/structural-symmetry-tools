# Deploying on TACC Lonestar6

## One-time setup

```bash
# Load the default Python module
module load python3/3.9.7

# Create a persistent virtual environment in your home directory
python3 -m venv $HOME/venvs/structural-symmetry-tools
source $HOME/venvs/structural-symmetry-tools/bin/activate

# Install the only dependency
pip install numpy

# Clone the repository to your work directory
cd $WORK
git clone https://github.com/OriginalGuru/structural-symmetry-tools.git
```

The virtual environment persists across login sessions and batch jobs.
You only need to do this once.

### ISOTROPY Software Suite (FINDSYM)

```bash
# Download and install
mkdir -p $WORK/software/isotropy
cd $WORK/software/isotropy
wget https://stokes.byu.edu/iso/iso.zip
unzip iso.zip
rm iso.zip
```

Add the following to `~/.bashrc` **above** the `[ -z "$PS1" ] && return` line
so that both interactive sessions and batch jobs pick it up:

```bash
export PATH=$PATH:$WORK/software/isotropy
export ISODATA=$WORK/software/isotropy/
```

Then reload and verify:
```bash
source ~/.bashrc
findsym findsym_sample.in   # should identify NaCl as Space Group 225 Fm-3m
```

## Running interactively

```bash
module load python3/3.9.7
source $HOME/venvs/structural-symmetry-tools/bin/activate

cd /path/to/your/phonon/calculation
python $WORK/structural-symmetry-tools/scripts/get_phonon_mode_charges.py > output.txt
```

## Running in a batch job

### Phonon analysis

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
source $HOME/venvs/structural-symmetry-tools/bin/activate

cd $SLURM_SUBMIT_DIR

python $WORK/structural-symmetry-tools/scripts/get_phonon_mode_charges.py > phonon_output.txt
```

Save as `run_phonon.sh` and submit with:
```bash
sbatch run_phonon.sh
```

### Space group identification with FINDSYM

```bash
#!/bin/bash
#SBATCH -J findsym
#SBATCH -A your_allocation
#SBATCH -p vm-small
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:10:00
#SBATCH -o findsym_%j.out
#SBATCH -e findsym_%j.err

module load python3/3.9.7
source $HOME/venvs/structural-symmetry-tools/bin/activate

cd $SLURM_SUBMIT_DIR

python $WORK/structural-symmetry-tools/scripts/run_findsym.py
```

Save as `run_findsym.sh` and submit with:
```bash
sbatch run_findsym.sh
```

## Notes

- All scripts are single-threaded and lightweight. The `vm-small` partition
  (or `development` for testing) is sufficient for all tasks.
- `vasprun.xml` can be several hundred MB for large supercells. Run from
  `$SCRATCH` or `$WORK`, not `$HOME` (quota is small).
- The virtual environment must live in `$HOME` (persistent). Do not put it
  in `$SCRATCH` (purged periodically).
- ISOTROPY/FINDSYM must be in `$WORK` or `$HOME`, not `$SCRATCH`, for the
  same reason.
- If you update the repository (`git pull`), no reinstallation is needed
  since there are no compiled components.
