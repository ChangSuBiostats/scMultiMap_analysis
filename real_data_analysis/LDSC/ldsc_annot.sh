BEDPATH='YourPath'
BIMPATH='BimPath'
ANNOTPATH='AnnotPath'
for j in proposed SCENT Signac;
do
        BEDFILE=${BEDPATH}/update_Microglia_"${j}_sig_sets.bed"
        for i in {1..22}
        do
                BIMFILE=${BIMPATH}/1000G.EUR.hg38.$i.bim
                ANNOT=${ANNOTPATH}/$j/annot/Microglia.$i.annot.gz
                make_annot_ver1.py --bed-file $BEDFILE --bimfile $BIMFILE --annot-file $ANNOT

        done
done
