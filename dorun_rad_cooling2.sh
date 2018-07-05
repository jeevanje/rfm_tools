# Args
case=$1
atm=$2
ctm=$3

zlevs=51
nqad=1
exec="rfm_NQAD=$nqad"
atmfile=${atm}_atm.txt
levfile=z${zlevs}.lev
gas=$(echo $case | cut -d'_' -f 1)

# Derived config vars
HDR="$atm atm, $GAS only, $ctm"
optFLG="flx zen opt vrt $ctm"
cooFLG="rad flx coo sfc $ctm"
tabFLG="tab $ctm"
ATM="../atm/$atmfile"
LEV="../lev/$levfile"
DIM=PLV
if [ $gas = h2o ]; then
    SPC="10 1500 1"
    GAS="H2O"
    HIT="../hit/h2o_1500cm-1.hit"
elif [ $gas = co2 ]; then
    SPC="500 850 1"
    GAS=CO2
    HIT="../hit/co2_500-850cm-1.hit"
fi

# Directories
mkdir -p $case 
for outdir in tab opt coo flx; do 
    if [ -d ${case}/$outdir ]; then
	rm ${case}/${outdir}/*
    elif [ ! -d ${case}/$outdir ]; then
	mkdir ${case}/${outdir}
    fi
done

# configure and run RFM 
for field in tab opt coo; do 
#field=tab  # for debugging 
    FLG=${field}FLG
    sed -e "/*HDR/ a\ $HDR"\
	-e "/*FLG/ a\ ${!FLG}"\
	-e "/*SPC/ a\ $SPC"\
	-e "/*GAS/ a\ $GAS"\
	-e "/*ATM/ a\ $ATM"\
	-e "/*DIM/ a\ $DIM"\
	-e "/*LEV/ a\ $LEV"\
	-e "/*HIT/ a\ $HIT"\
	< rfm_${field}.drv_pre > $case/rfm_${field}.drv

    cd $case
    cp rfm_${field}.drv rfm.drv
    ../src/$exec
    mv rfm.log rfm_${field}.log
    cd ../
done

Rscript asc2ncdf.R $case 