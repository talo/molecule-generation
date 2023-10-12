#!/bin/bash
time_hash=$(date +%s)

# RUN THIS FROM ROOT OF GITHUB REPO
export TF_CPP_MIN_LOG_LEVEL=3

# store STDOUT
mkdir -p logs/

# note: please run within a "screen" session so that the run continues even if SSH connection is cut off
run_name="50parts_20steps_8res_5poses_8exh_seed8198_top10smis_v3"
python3 optimise/run.py \
    --run_name ${run_name} \
    --num_particles 50 \
    --num_epochs 20 \
    --seed 8198 \
    --num_top_smiles_to_init 10 \
    --num_poses 5 \
    --out_folder "results" \
    --scaffold "CNC(=O)C1=C(NC2=CC=CC=C2OC)C=C(NC(=O)C2CC2)N=N1" \
    --residues_of_interest "VAL629" "ARG715" "TRP718" "THR555" "CYS675" "PHE594" "PHE595" "GLU627" \
    --necessary_residues "VAL629" \
    --ref_smi_excel_path "data/JAK2_IC50_data_Sep18.xlsx" \
    --protein_pdb_path "data/mutant_PRT1008596.align.fixed.addH.pdb" \
    --autobox_ligand_path "data/PRT1008596.sdf" \
    2>&1 | tee logs/${run_name}_${time_hash}.log

#    "VAL629", "ARG715", "TRP718", "THR555", "CYS675", 

# hydrophobic pocket (on top of ARG715 and CYS675)
# PHE594, PHE595

# hinge h bonds
# VAL629, GLU627 

# other residues
#    "LYS581", "GLN626", "LYS677", "SER698",
