#!/home/rcf-40/hli465/bin/Rscript
#PBS -q cmb
#PBS -l walltime=100:00:00
#PBS -e /home/cmb-11/plr/hli465/POPRES/output
#PBS -o /home/cmb-11/plr/hli465/POPRES/output
#PBS -l nodes=1:ppn=1
#PBS -l mem=100gb,pmem=100gb,vmem=100gb

coded <- as.matrix(read.table("coded_data_chr1.txt"))
M <- rowMeans(coded,na.rm=TRUE)
coded <- coded-M
cov <- cov(coded,use="pairwise")
write.table(cov,"cov_data_for_chr1.txt",sep="\t")

