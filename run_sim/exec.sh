#PBS -N DrugSim_Qui
#PBS -l nodes=1:ppn=10
#PBS -l walltime=20000:00:00
#PBS -e stderr.log
#PBS -o stdout.log
#Specific the shell types
#PBS -S /bin/bash
#Specific the queue type
#PBS -q dque

cd $PBS_O_WORKDIR
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS nodes

# Use this to export the library path
export LD_LIBRARY_PATH=/opt/prog/sundial/lib64:$LD_LIBRARY_PATH

find . -name "*.plt" -type f -delete
rm -rf *.log result
mpirun -machinefile $PBS_NODEFILE -np $NPROCS /home/cml/marcell/DrugSimulation_ALI_single/bin/drug_sim \
	-input_deck EDISON_INPUT_DECK_SINGLE.txt \
	-hill_file /home/cml/marcell/DrugSimulation_ALI_single/bin/drugs/quinidine/IC50_samples10.csv > logfile
