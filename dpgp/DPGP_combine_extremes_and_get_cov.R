#!/home/rcf-40/hli465/bin/Rscript
#PBS -q cmb
#PBS -l walltime=200:00:00
#PBS -e /home/cmb-11/plr/hli465/Dpgp/all_samples_Chr3L/output
#PBS -o /home/cmb-11/plr/hli465/Dpgp/all_samples_Chr3L/output
#PBS -l nodes=1:ppn=1
#PBS -l mem=100gb,pmem=100gb,vmem=100gb
chr=read.table("/home/cmb-11/plr/hli465/Dpgp/all_samples_Chr3L/coded_data_for_all_samples_seqs_both_low_NAs_Chr3L_with_SNP_Pos.txt",header=TRUE,row.names=1)
neibors= read.table("/home/cmb-11/plr/hli465/Dpgp/all_samples_Chr3L/neibors_win103_Chr3L.txt")
win=10^3
getchunk<-function(x){
        chunk=(win*(x-1)+1):(win*x)
        return(chunk)
}
chr_chunck_cov<-function(k){
        temp1<-sapply(c(sort(neibors[,k])),getchunk)
        chr_temp<-chr[c(temp1),]
        temp=data.matrix(chr_temp)
        data=temp
    M=rowMeans(data,na.rm=TRUE)
    M=rep(M,times=ncol(data))
    M=matrix(M,nrow=nrow(data),ncol=ncol(data),byrow=FALSE)
    data=data-M
    cov=cov(data,use="pairwise")
    return(cov)
}
cov1=chr_chunck_cov(1)
write.table(cov1,"/home/cmb-11/plr/hli465/Dpgp/all_samples_Chr3L/Chr3L_recomchunk_win103_cov1.txt",sep="\t")
cov2=chr_chunck_cov(2)
write.table(cov2,"/home/cmb-11/plr/hli465/Dpgp/all_samples_Chr3L/Chr3L_recomchunk_win103_cov2.txt",sep="\t")
cov3=chr_chunck_cov(3)
write.table(cov3,"/home/cmb-11/plr/hli465/Dpgp/all_samples_Chr3L/Chr3L_recomchunk_win103_cov3.txt",sep="\t")
