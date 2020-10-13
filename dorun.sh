# Args
case=$1
gas=$2
atm=$3
SFC=$4
dk=$5
save_sigma=$6

ctm=ctm
exec="rfm_NQAD=1"
atmfile=${atm}.atm
levfile=${atm}.lev
#zlevs=$(sed -n "3p" atm/$atmfile | awk '{print($1)}')
if [ ! -f lev/${levfile} ]; then
    Rscript make_lev.R $atm
fi

# Derived config vars
HDR="$atm atm, $ctm"
optFLG="flx zen opt vrt $ctm"
cooFLG="rad flx coo sfc $ctm"
tabFLG="tab $ctm"
rfmdir=/home/nadirj/17rad_cooling2/rfm
scratchdir=/tigress/nadirj/17rad_cooling2_data
ATM="${rfmdir}/atm/$atmfile"
LEV="${rfmdir}/lev/$levfile"
DIM=PLV

if [ $gas = h2o ]; then
    SPC="10 1500 $dk"
    GAS="H2O"
    HIT="${rfmdir}/hit/h2o_1500cm-1.hit"
elif [ $gas = co2 ]; then
    SPC="500 850 $dk"
    GAS=CO2
    HIT="${rfmdir}/hit/co2_500-850cm-1.hit"
elif [ $gas = both ]; then
    SPC="10 1500 $dk"
    GAS="H2O CO2"
    HIT="${rfmdir}/hit/h2o_co2_1500cm-1.hit"
fi
HDR="$atm atm, $gas, $ctm"


# Directories
mkdir -p ${scratchdir}/$case
ln -s ${scratchdir}/$case $case 
for outdir in tab opt coo flx; do 
    if [ -d ${case}/$outdir ]; then
	rm ${case}/${outdir}/*
    elif [ ! -d ${case}/$outdir ]; then
	mkdir -p ${case}/${outdir}
    fi
done

# configure and run RFM 
for field in tab opt coo; do 
    echo "RFM ${field} run"
    FLG=${field}FLG
    sed -e "/*HDR/ a\ $HDR"\
	-e "/*FLG/ a\ ${!FLG}"\
	-e "/*SPC/ a\ $SPC"\
	-e "/*GAS/ a\ $GAS"\
	-e "/*ATM/ a\ $ATM"\
	-e "/*DIM/ a\ $DIM"\
	-e "/*LEV/ a\ $LEV"\
        -e "/*SFC/ a\ $SFC"\
	-e "/*HIT/ a\ $HIT"\
	< rfm_${field}.drv_pre > $case/rfm_${field}.drv

    cd $case
    cp rfm_${field}.drv rfm.drv
    ${rfmdir}/src/$exec > rfm_${field}.out
    mv rfm.log rfm_${field}.log
    cd ../
done

Rscript asc2ncdf.R $case $gas $save_sigma
