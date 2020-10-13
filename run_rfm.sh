# Args
gas=$1
atm=$2
dk=$3   #cm-1
save_sigma=F
ctm=ctm
exec="rfm"
case=${gas}_${atm}_dk${dk}

# Derived script vars
atmfile=${atm}.atm
zlevs=$(sed -n "3p" atm/$atmfile | awk '{print($1)}')
levfile=z${zlevs}.lev

# Derived config vars
HDR="$atm atm, $ctm"
optFLG="flx zen opt vrt $ctm"
cooFLG="rad flx coo sfc $ctm"
tabFLG="tab $ctm"
ATM="${PWD}/atm/$atmfile"
LEV="${PWD}/lev/$levfile"
DIM=PLV
SFC=$(grep -A 1 "*TEM" atm/${atmfile} | sed -n "2p" | awk '{print($1)}' | rev | cut -c 2- | rev)

if [ $gas = h2o ]; then
    SPC="10 1500 $dk"
    GAS="H2O"
    HIT="${PWD}/hit/h2o_1500cm-1.hit"
elif [ $gas = co2 ]; then
    SPC="500 850 $dk"
    GAS=CO2
    HIT="${PWD}/hit/co2_500-850cm-1.hit"
elif [ $gas = both ]; then
    SPC="10 1500 1"
    GAS="H2O CO2"
    HIT="${PWD}/hit/h2o_co2_1500cm-1.hit"
fi
HDR="$atm atm, $gas, $ctm"


# Directories
echo "Making directories"
rfmdir=$PWD
projectdir=${rfmdir%/rfm}
project=${projectdir##*/}
datadir=/tigress/nadirj/${project}_data
mkdir -p $datadir
mkdir -p $datadir/${case}
ln -s ${datadir}/$case $case 
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
        -e "/*SFC/ a\ $SFC"\
	-e "/*HIT/ a\ $HIT"\
	< rfm_${field}.drv_pre > $case/rfm_${field}.drv

    cd $case
    cp rfm_${field}.drv rfm.drv
    $rfmdir/src/$exec
    mv rfm.log rfm_${field}.log
    cd $rfmdir
done

Rscript asc2ncdf.R $case $gas $save_sigma
