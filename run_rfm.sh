# Args
gas=$1
atm=$2
ctm=$3
dk=$4
save_sigma=$5
case=$6
exec="rfm_new"
mix=mix
kmin=500  # cm-1
kmax=850  # cm-1
#dk=0.1   # cm-1

if [ $ctm == "noctm" ] ; then
   ctm=""
fi

if [ $case == "" ] ; then
    echo "Case is blank. Exiting."
    exit
fi

echo "Begin run_rfm.sh"

# Derived script vars
atmfile=${atm}.atm
zlevs=$(sed -n "3p" atm/$atmfile | awk '{print($1)}')
levfile=z${zlevs}.lev

# Derived config vars
HDR="$atm atm, $ctm"
optFLG="flx zen opt vrt $ctm $mix"
cooFLG="rad flx coo sfc $ctm $mix"
tabFLG="tab $ctm $mix"
ATM="${PWD}/atm/$atmfile"
LEV="${PWD}/lev/$levfile"
DIM=PLV
SFC=$(grep -A 1 "*TEM" atm/${atmfile} | sed -n "2p" | awk '{print($1)}' | rev | cut -c 2- | rev)

if [ $gas = h2o ]; then
    SPC="$kmin $kmax $dk"
    GAS="H2O"
    HIT="${PWD}/hit/h2o_1500cm-1.hit"
elif [ $gas = co2 ]; then
    SPC="$kmin $kmax $dk"
    GAS=CO2
    HIT="${PWD}/hit/co2_500-850cm-1.hit"
elif [ $gas = both ]; then
    SPC="$kmin $kmax $dk"
    GAS="H2O CO2"
    HIT="${PWD}/hit/h2o_co2_1500cm-1.hit"
fi
HDR="$atm atm, $gas, $ctm"


# Directories
rfmdir=$PWD
projectdir=${rfmdir%/rfm}
project=${projectdir##*/}
datadir=/tigress/nadirj/${project}_data
mkdir -p $datadir
if [ -d $datadir/${case} ] ; then
    rm -rf $datadir/${case}
    rm -rf ${case}
    echo "Removing old case directories"
fi
echo "Making directories"
mkdir -p $datadir/${case}
ln -s ${datadir}/$case $case 
for outdir in tab opt coo flx; do 
    mkdir -p ${case}/${outdir}
done

# configure and run RFM 
for field in tab opt coo; do 
#field=tab  # for debugging 
    echo "${field} run"
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

echo "END run_rfm.sh"
