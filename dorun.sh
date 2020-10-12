# Args
case=co2_Ts300_gamma7_Ttp150_dk0.1
gas=co2
atm=Ts300_rh0.75_gamma7_Ttp150
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
ATM="${rfmdir}/atm/$atmfile"
LEV="${rfmdir}/lev/$levfile"
DIM=PLV
SFC=$(grep -A 1 "*TEM" atm/${atmfile} | sed -n "2p" | awk '{print($1)}' | rev | cut -c 2- | rev)
#SFC=260
if [ $gas = h2o ]; then
    SPC="10 1500 1"
    GAS="H2O"
    HIT="${rfmdir}/hit/h2o_1500cm-1.hit"
elif [ $gas = co2 ]; then
    SPC="500 850 0.1"
    GAS=CO2
    HIT="${rfmdir}/hit/co2_500-850cm-1.hit"
elif [ $gas = both ]; then
    SPC="10 1500 1"
    GAS="H2O CO2"
    HIT="${rfmdir}/hit/h2o_co2_1500cm-1.hit"
fi
HDR="$atm atm, $gas, $ctm"


# Directories
scratchdir=/tigress/nadirj/17rad_cooling2_data
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
#field=tab  # for debugging 
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

Rscript asc2ncdf.R $case $gas
