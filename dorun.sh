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

# Derived vars
if [ $gas = h2o ]; then
    OPTFLG="flx zen opt vrt $ctm"
    COOFLG="rad flx coo sfc $ctm"
    SPC="10 1500 1"
    GAS=H2O
    ATM="../atm/$atmfile"
    LEV="../lev/$levfile"
    HIT="../hit/h2o_1500cm-1.hit"
elif [ $gas = co2 ]; then
    OPTFLG="flx zen opt vrt"
    COOFLG="rad flx coo sfc"
    SPC="500 850 1"
    GAS=CO2
    ATM="../atm/$atmfile"
    LEV="../lev/$levfile"
    HIT="../hit/co2_500-850cm-1.hit"
fi

HDR="$atm atm, $GAS only, $ctm"

sed -e "/*HDR/ a\ $HDR"\
    -e "/*FLG/ a\ $OPTFLG"\
    -e "/*SPC/ a\ $SPC"\
    -e "/*GAS/ a\ $GAS"\
    -e "/*ATM/ a\ $ATM"\
    -e "/*LEV/ a\ $LEV"\
    -e "/*HIT/ a\ $HIT"\
    < rfm_opt.drv_pre > $case/rfm_opt.drv

sed -e "/*HDR/ a\ $HDR"\
    -e "/*FLG/ a\ $COOFLG"\
    -e "/*SPC/ a\ $SPC"\
    -e "/*GAS/ a\ $GAS"\
    -e "/*ATM/ a\ $ATM"\
    -e "/*LEV/ a\ $LEV"\
    -e "/*HIT/ a\ $HIT"\
    < rfm_coo.drv_pre > $case/rfm_coo.drv

cd $case

# opt calculation
cp rfm_opt.drv rfm.drv
../src/$exec
mv rfm.log rfm_opt.log

# coo and flx calculation
cp rfm_coo.drv rfm.drv
../src/$exec
mv rfm.log rfm_coo.log

cd ../
Rscript asc2ncdf.R $case 