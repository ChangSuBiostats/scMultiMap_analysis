H2PATH='H2Path'
CELLPATH='CellPath'
BASEPATH='BasePath'
WEIGHTSPATH='WeightsPath'
FREQPATH='FreqPath'
OUTDIRPATH='OutdirPath'
for k in Bellenguez Jansen Kunkle;
do
        H2=${H2PATH}/${k}_AD.sumstats.gz
        for j in proposed SCENT Signac;
        do
                CELL=${CELLPATH}/$j/annot_overlap/Microglia.
                BASE=${BASEPATH}/baselineLD.
                WEIGHTS=${WEIGHTSPATH}/weights.hm3_noMHC.
                FREQ=${FREQPATH}/1000G.EUR.hg38.
                OUTDIR=${OUTDIRPATH}/$j/heritability/$k/${k}_AD_Heritability_B22
                ldsc.py --h2 $H2 --ref-ld-chr $CELL,$BASE --w-ld-chr $WEIGHTS --frqfile-chr $FREQ --overlap-annot --print-coefficients --out $OUTDIR
        done
done

