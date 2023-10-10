# setup

## miniconda set up on linux / ubuntu
```bash
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
```

## install python packages
```bash
# use either conda or miniconda, both work
conda env create -f environment.yml
conda activate moler-env
pip install molecule-generation

# install plip
git clone https://github.com/pharmai/plip.git
cd ./plip/
# go to setup.py then comment out the line that says "openbabel" under install_requires (should be around line 23)
python3 setup.py install
```

NOTE: this assumes you already have GNINA working as a binary executable. please install GNINA separately.


# run the scaffold-constrained optimisation
Please first run `export TF_CPP_MIN_LOG_LEVEL=3` to silence tensorflow's bazillion error messages. 

The entrypoint is at: `python3 optimise/run.py`. See `script/optimise.sh` & the python script itself for arguments.

Run with:
```bash
bash script/optimise.sh
```

Modify the code & parameters as you desire.

# Outputs
With each epoch in the optimisation, the code will save several csv files into the output folder. by default, this folder is `./results/{run_name}_{time_hash}`.

## scores
The `best_solutions.csv` contains all SMILES generated during the entire run, with the residue interaction score, sorted from best SMILES at the top (1st row) in decreasing order.
```
smiles,fitness
COc1ccccc1Nc1cc(NC(=O)C2(CO)CC2)nnc1C(=O)NCc1ccc(-c2cc(C)c(Br)c(C)n2)nc1,0.7777777777777779
Cc1ccc(-c2ccc(Nc3cc(NC(=O)C4CC4)nnc3C(=O)NCc3ccc(C)nc3C)c(OC(CO)C(=O)O)c2)cn1,0.7777777777777779
O=C(NCc1ccc(-c2ccc(Br)cn2)nc1)c1nnc(NC(=O)C2CC2)cc1Nc1ccccc1OCCCCCCO,0.6666666666666665
COc1ccccc1Nc1cc(NC(=O)C2(C)CC2)nnc1C(=O)NCc1cc(-c2cc(C)cc(C)n2)ccn1,0.6666666666666665
COc1ccccc1Nc1cc(NC(=O)C2CC2)nnc1C(=O)NCc1cc(Oc2ccccn2)ccn1,0.6666666666666665
COc1ccccc1Nc1cc(NC(=O)C2CC2)nnc1C(=O)NCc1cc(OCc2ccc(Br)cn2)ccn1,0.6666666666666665
COc1ccccc1Nc1cc(NC(=O)C2(C)CC2)nnc1C(=O)NCc1ccc(-c2cc(C)cc(C)n2)cn1,0.6666666666666665
```

## docked SDFs

At the same time, the code will also save all the GNINA docked ligand SDFs into the temp dir, by default it is `./temp_dirs/moler_oracle_temp_{run_name}_{time_hash}`. in this folder, you will see SDFs named `best_score{score}_pose{pose}_{hash}.sdf`. 

The score is the residue interaction score ranging from 0 (worst) to 1 (best). Pose is the pose index that gave the best score (and corresponds to the conformer inside this SDF). Hash is just a random hash.

Inside the SDF, you will see the usual GNINA properties, followed by 2 more properties, `residues` and `smiles`.
`residues` refers to all residues that interacted with this docked pose, as returned by `plip`.

```bash
>  <minimizedAffinity>  (1) 
-8.33303

>  <CNNscore>  (1) 
0.8733225465

>  <CNNaffinity>  (1) 
7.5456223488

>  <CNN_VS>  (1) 
6.5897622108

>  <CNNaffinity_variance>  (1) 
0.0177119523

>  <residues>  (1) 
LEU551,GLN553,LEU579,LYS581,GLU627,PHE628,VAL629,ASP635,THR636,LYS677

>  <smiles>  (1) 
CNC(=O)c1nnc(NC(=O)C2CC2)cc1Nc1cccc(-c2cnn(C3CN(Cc4ccn(C)n4)C3)c2)c1OC

```
