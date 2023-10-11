#!/bin/bash
time_hash=$(date +%s)

# RUN THIS FROM ROOT OF GITHUB REPO
export TF_CPP_MIN_LOG_LEVEL=3

# store STDOUT
mkdir -p logs/

# fixed hyperball clipping of latent vectors
run_name="100parts_10steps_5res_10poses_8exh_seed42_top20smis_fixed"
python3 optimise/run.py \
    --run_name ${run_name} \
    --num_particles 100 \
    --num_epochs 10 \
    --seed 42 \
    --num_top_smiles_to_init 20 \
    --num_poses 10 \
    --out_folder "results" \
    --scaffold "CNC(=O)C1=C(NC2=CC=CC=C2OC)C=C(NC(=O)C2CC2)N=N1" \
    --residues_of_interest "VAL629" "ARG715" "TRP718" "THR555" "CYS675" \
    --necessary_residues "VAL629" \
    --ref_smi_excel_path "data/JAK2_IC50_data_Sep18.xlsx" \
    --protein_pdb_path "data/mutant_PRT1008596.align.fixed.addH.pdb" \
    --autobox_ligand_path "data/PRT1008596.sdf" \
    2>&1 | tee logs/${run_name}_${time_hash}.log

# residues_of_interest
# {
#    "VAL629", "ARG715", "TRP718", "THR555", "CYS675", 
#    "LYS581", "GLN626", "LYS677", "SER698",
# },
