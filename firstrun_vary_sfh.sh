#!/bin/bash
# should be in code directory

CMD="nice -n 19 python tpagb_calibration/sfhs/vary_sfh -n 8"
LOC="/home/rosenfield/research/TP-AGBcalib/SNAP/varysfh"
EXT="_vsfh.sh"
# only one filter in SNAP/tables/snap_galaxies.dat
for galaxy in eso540-030 scl-de1 ugc-04459 ugc-4305-2 ugc-5139 ugc8508 ddo78 kdg73 kkh37 ngc3741 hs117 ngc2403-deep
do
    echo $galaxy
    #python /home/rosenfield/research/TP-AGBcalib/code/TPAGB-calib/tpagb_calibration/prepare_data.py -n $1 -d $galaxy
    echo "$CMD $LOC/$galaxy/$galaxy.vsfhinp > $LOC/$galaxy/$galaxy.vsfhout" > $galaxy$EXT
    #python -m tpagb_calibration.plotting.plotting /home/rosenfield/research/TP-AGBcalib/SNAP/varysfh/$galaxy/$galaxy.plotinp
done

# more than one filter in SNAP/tables/snap_galaxies.dat
FILT="_555"
GALAXY="ic2574-sgs"
echo $GALAXY
#python /home/rosenfield/research/TP-AGBcalib/code/TPAGB-calib/tpagb_calibration/prepare_data.py -n $1 -df 555 ic2574-sgs
echo "$CMD $LOC/$GALAXY/$GALAXY$FILT.vsfhinp > $LOC/$GALAXY/$GALAXY$FILT.vsfhout" > $GALAXY$FILT$EXT

FILT="_475"
for galaxy in ugca292 ngc300-wide1
do
    echo $galaxy
    #python /home/rosenfield/research/TP-AGBcalib/code/TPAGB-calib/tpagb_calibration/prepare_data.py -n $1 -df 475 $galaxy
    echo "$CMD $LOC/$galaxy/$galaxy$FILT.vsfhinp > $LOC/$galaxy/$galaxy$FILT.vsfhout &" > $galaxy$FILT$EXT
done

FILT="_606"
for galaxy in ugca292 ngc300-wide1 ddo82 ngc4163
do
    echo $galaxy
    #python /home/rosenfield/research/TP-AGBcalib/code/TPAGB-calib/tpagb_calibration/prepare_data.py -n $1 -df 606 $galaxy
    echo "$CMD $LOC/$galaxy/$galaxy$FILT.vsfhinp > $LOC/$galaxy/$galaxy$FILT.vsfhout &" > $galaxy$FILT$EXT
done
