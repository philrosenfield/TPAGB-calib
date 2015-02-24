#!/bin/bash
# should be in code directory
# usage: ./firstrun_vary_sfh.sh nprocs

#ipcluster start -n=$1 &
CMD="nice -n 19 python -m tpagb_calibration.sfhs.vary_sfh"
LOC="/home/rosenfield/research/TP-AGBcalib/SNAP/varysfh"
# only one filter in SNAP/tables/snap_galaxies.dat
for galaxy in eso540-030 scl-de1 ugc-04459 ugc-4305-2 ugc-5139 ugc8508 ddo78 kdg73 kkh37 ngc3741 hs117 ngc2403-deep
do
    echo $galaxy
    #python /home/rosenfield/research/TP-AGBcalib/code/TPAGB-calib/tpagb_calibration/prepare_data.py -n $1 -d $galaxy
    echo "$CMD $LOC/$galaxy/$galaxy.vsfhinp > $loc/$galaxy/$galaxy.vsfhout"
    #> $galaxy_vsfh.sh
    #python -m tpagb_calibration.plotting.plotting /home/rosenfield/research/TP-AGBcalib/SNAP/varysfh/$galaxy/$galaxy.plotinp
done
# more than one filter in SNAP/tables/snap_galaxies.dat
echo "ic2574-sgs"
#python /home/rosenfield/research/TP-AGBcalib/code/TPAGB-calib/tpagb_calibration/prepare_data.py -n $1 -df 555 ic2574-sgs
echo "$CMD $LOC/ic2574-sgs/ic2574-sgs_555.vsfhinp > $loc/ic2574-sgs/ic2574-sgs_555.vsfhout"
#> ic2574-sgs_555_vsf.sh

for galaxy in ugca292 ngc300-wide1
do
    echo $galaxy
    #python /home/rosenfield/research/TP-AGBcalib/code/TPAGB-calib/tpagb_calibration/prepare_data.py -n $1 -df 475 $galaxy
    echo "$CMD $LOC/$galaxy/$galaxy_475.vsfhinp > $loc/$galaxy/$galaxy_475.vsfhout &"
    #> $galaxy_475_vsfh.sh
done

for galaxy in ugca292 ngc300-wide1 ddo82 ngc4163
do
    echo $galaxy
    #python /home/rosenfield/research/TP-AGBcalib/code/TPAGB-calib/tpagb_calibration/prepare_data.py -n $1 -df 606 $galaxy
    echo "$CMD $LOC/$galaxy/$galaxy_606.vsfhinp > $loc/$galaxy/$galaxy_606.vsfhout &"
    #> $galaxy_606_vsfh.sh
done

#ipcluster stop