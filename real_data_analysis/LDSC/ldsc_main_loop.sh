BIMPATH='BimPath'
ANNOTPATH='AnnotPath'
OUTPATH='OutPath'
SNPPATH='SnpPath'
for j in proposed SCENT Signac;
do
        for i in {1..22}
        do
                BIMFILE=${BIMPATH}/1000G.EUR.hg38.$i
                ANNOT=${ANNOTPATH}/$j/annot/Microglia.$i.annot.gz
                OUT=${OUTPATH}/$j/annot_overlap/Microglia.$i
                SNPS=${SNPPATH}/hm3_no_MHC.list

        ldsc.py --l2 --bfile $BIMFILE --ld-wind-cm 1 --annot $ANNOT --thin-annot --out $OUT --print-snps $SNPS

        done
done

