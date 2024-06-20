## I did this in the uni cluster, wihtout using slurm, still fairly fast.
## I can also use slrum in interactive mode if more computer power is needed
# activate where msprime is installed in python
# srun --pty --nodes=1 --cpus-per-task=4 --mem=16G --time=04:00:00 /bin/bash

source $HOME/my_python3_env/my_env/bin/activate
cd $WORK2/cod_sim/00_burnin
# export pymsprimeBurnin=/Users/khan/Desktop/WORK/Project_Fisheries_induced_evol/cod_sim/00_burnin/
export pymsprimeBurnin=/gxfs_work2/geomar/smomw426/cod_sim/scripts/00_msprime_burnin.py
export mutationRate=3.5e-9
export recombinationRate=3.11e-8

export chrLen=30000000
export ndiploid=1000
# run a set of burn-in sims
#  -n <Ne diploid> -L <length in bp> -m <mutation rate per site> -r <recombination rate per site> -o <outputfile>
#for i in {1..20}
for i in 1
do
	$pymsprimeBurnin -n $ndiploid -L $chrLen -m $mutationRate -r $recombinationRate -o msprime_burnin_n${ndiploid}_${chrLen}_$i.vcf
#	$pymsprimeBurnin -n 50000 -L $chrLen -m $mutationRate -r $recombinationRate -o msprime_burnin_n${ndiploid}_${chrLen}_$i.vcf
done

# output as tree sequence
$pymsprimeBurnin -n $ndiploid -L $chrLen -m $mutationRate -r $recombinationRate -o msprime_burnin_n${ndiploid}_${chrLen}_1.tree
