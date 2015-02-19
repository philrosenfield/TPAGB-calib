#!/bin/bash
# should be in code directory
# usage: ./firstrun_vary_sfh.sh

# only one filter in SNAP/tables/snap_galaxies.dat
for galaxy in eso540-030 scl-de1 ugc-04459 ugc-4305-2 ugc-5139 ugc8508 ddo78 kdg73 kkh37 ngc3741 hs117 ngc2403-deep
do
    echo $galaxy
    #python /home/rosenfield/research/TP-AGBcalib/code/TPAGB-calib/tpagb_calibration/prepare_data.py -n $1 -d $galaxy
    python -m tpagb_calibration.sfhs.vary_sfh --dry_run /home/rosenfield/research/TP-AGBcalib/SNAP/varysfh/$galaxy/$galaxy.vsfhinp > /home/rosenfield/research/TP-AGBcalib/SNAP/varysfh/$galaxy/$galaxy.vsfhout
    wait
    #python -m tpagb_calibration.plotting.plotting /home/rosenfield/research/TP-AGBcalib/SNAP/varysfh/$galaxy/$galaxy.plotinp
done
wait
# more than one filter in SNAP/tables/snap_galaxies.dat
echo "ic2574-sgs"
#python /home/rosenfield/research/TP-AGBcalib/code/TPAGB-calib/tpagb_calibration/prepare_data.py -n $1 -df 555 ic2574-sgs
python -m tpagb_calibration.sfhs.vary_sfh --dry_run /home/rosenfield/research/TP-AGBcalib/SNAP/varysfh/c2574-sgs/c2574-sgs_555.vsfhinp > /home/rosenfield/research/TP-AGBcalib/SNAP/varysfh/c2574-sgs/c2574-sgs_555.vsfhout

for galaxy in ugca292 ngc300-wide1
do
    echo $galaxy
    #python /home/rosenfield/research/TP-AGBcalib/code/TPAGB-calib/tpagb_calibration/prepare_data.py -n $1 -df 475 $galaxy
    python -m tpagb_calibration.sfhs.vary_sfh --dry_run /home/rosenfield/research/TP-AGBcalib/SNAP/varysfh/$galaxy/$galaxy_475.vsfhinp > /home/rosenfield/research/TP-AGBcalib/SNAP/varysfh/$galaxy/$galaxy_475.vsfhout &
    wait
done

wait

for galaxy in ugca292 ngc300-wide1 ddo82 ngc4163
do
    echo $galaxy
    #python /home/rosenfield/research/TP-AGBcalib/code/TPAGB-calib/tpagb_calibration/prepare_data.py -n $1 -df 606 $galaxy
    python -m tpagb_calibration.sfhs.vary_sfh --dry_run /home/rosenfield/research/TP-AGBcalib/SNAP/varysfh/$galaxy/$galaxy_606.vsfhinp > /home/rosenfield/research/TP-AGBcalib/SNAP/varysfh/$galaxy/$galaxy_606.vsfhout &
    wait
done


