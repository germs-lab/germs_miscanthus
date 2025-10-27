#!/bin/bash

#SBATCH --nodes=1   # Number of nodes to use
#SBATCH --ntasks-per-node=32   # Use 32 processor cores per node 
#SBATCH --time=3-0:0:0   # Walltime limit (DD-HH:MM:SS)
#SBATCH --mem=32G   # Maximum memory per node
#SBATCH --job-name="mxg_BLAST"   # Job name to display in squeue
#SBATCH --mail-user=bolivar@iastate.edu   # Email address
#SBATCH --mail-type=ALL   # Send an email when the job starts
#SBATCH --chdir="/work/adina/bolivar/germs_miscanthus" # Set the working directory of the batch script to directory before it is executed. Absolute or relative paths
#SBATCH --output="/work/adina/bolivar/germs_miscanthus/hpc/slurm-%j-mxg_asv_blast.out"   # Job standard output file (%j will be replaced by the slurm job id)
#SBATCH --error="/work/adina/bolivar/germs_miscanthus/hpc/slurm-%j-mxg_asv_blast.out"   # Job standard error file (%j will be replaced by the slurm job id)


#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK # Set OMP_NUM_THREADS to the number of CPUs per task we asked for.

##Modules/Singularity
module purge
#module load micromamba/1.4.2-lcemqbe # Latest in Nova
module load blast-plus/2.13.0-py310-irl7uxo

# Basic session info
echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS
echo sbatch invoked from: $SLURM_SUBMIT_DIR

# Run BLAST job
# Paths
PROJ_DIR=/work/adina/bolivar/germs_miscanthus
SCRIPT_PATH=$PROJ_DIR/scripts/blast_nt_remote.sh


# Define regions and their corresponding input files
declare -A REGIONS=(
    ["16S"]="$PROJ_DIR/data/output/processed/sequences/mxg_16S_combined_asv_renamed.fa"
    ["ITS"]="$PROJ_DIR/data/output/processed/sequences/mxg_ITS_combined_asv_renamed.fa"
    ["AMF"]="$PROJ_DIR/data/output/processed/sequences/mxg_AMF_combined_asv_renamed.fa"
)

# Loop through each region
for REGION in "${!REGIONS[@]}"; do
    echo "Starting BLAST analysis for $REGION..."
    
    # Set input and output paths
    IN="${REGIONS[$REGION]}"
    OUT="$PROJ_DIR/data/output/processed/sequences/mxg_${REGION}_blast_results.tsv"
    
    # Check if input file exists
    if [[ -f "$IN" ]]; then
        echo "Processing $REGION: $IN -> $OUT"
        
        # Execute BLAST script with region-specific parameters
        $SCRIPT_PATH "$REGION" "$IN" "$OUT"
        
        echo "Completed BLAST analysis for $REGION"
        echo "----------------------------------------"
    else
        echo "Warning: Input file not found for $REGION: $IN"
        echo "Skipping $REGION analysis"
        echo "----------------------------------------"
    fi
done

# End job
module purge
echo "All BLAST jobs completed"

