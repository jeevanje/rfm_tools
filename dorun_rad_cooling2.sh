# Args
gas=$1    # "h2o", "co2", or "both"
Ts=$2     # K
rh=$3     # absolute
gamma=$4  # K/km
ctm=$5    # "ctm" or ""

# Derived vars
atm="${Ts}_rh${rh}_gamma${gamma}"
case="${gas}_${atm}${ctm}"
atmfile=${atm}_atm.txt
if [ ! -f atm/${atmfile} ];then
    printf "Creating atm file"\n
    cd atm
    Rscript make_simple_atm.R $Ts $rh $gamma
    cd ../
else
    printf "atm file already exists"\n
fi
 
# Fixed RFM params
#zlevs=51
zlevs=351
levfile=z${zlevs}.lev
nqad=1
exec="rfm_NQAD=$nqad"

# Derived RFM config vars
optFLG="flx zen opt vrt $ctm"
cooFLG="rad flx coo sfc $ctm"
tabFLG="tab $ctm"
ATM="../atm/$atmfile"
LEV="../lev/$levfile"
DIM=PLV
#DIM1="1 1000"
#DIM2="1 $Ts"
SFC=$Ts
if [ $gas = h2o ]; then
    SPC="10 1500 1"
    GAS="H2O"
    HIT="../hit/h2o_1500cm-1.hit"
elif [ $gas = co2 ]; then
    SPC="500 850 1"
    GAS=CO2
    HIT="../hit/co2_500-850cm-1.hit"
elif [ $gas = both ]; then
    SPC="10 1500 1"
    GAS="H2O CO2"
    HIT="../hit/h2o_co2_1500cm-1.hit"
fi
HDR="$atm atm, $GAS only, $ctm"

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
#field=tab  # for debugging 
for field in tab opt coo; do 
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
    ../src/$exec
    mv rfm.log rfm_${field}.log
    cd ../
done

Rscript asc2ncdf.R $case $gas
