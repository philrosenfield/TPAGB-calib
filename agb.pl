#!/usr/bin/perl
# Script PERL to fit the RGB+AGB LF
use CGI qw/:standard/;
use File::Basename;

$runitreally=1 ; # set to 0 to skip trilegal and spread
$magabovetip = 2 ; # mag cut above tip for star counts
$magbelowtip = 2 ; # mag cut below tip for star counts
$meanage = 14 ; # mean age of population in Gyr
$deltage = 1 ; # delta age of population in Gyr
$factdeltalogz = 0.5 ; # multiplicative factor for delta[M/H] of population
$nsigma = 1 ; # how far from 1 sigma level of RGB Poisson distribution

$argc = @ARGV;          # get the number of arguments
if ($argc==0) {
  print "Usage: agb.pl <galdataname>\n";
  print "Example: agb.pl 9884_M81K61\n";
  die;
} else {
  $galdataname = @ARGV[0];
  print "Doing it for $photfile\n";
}

$galname = $galdataname ;
$galname =~ s/[0-9]+_// ; #eliminates initial number and _

print "$galname\n";
if (-e "$galdataname/proc_new_dol/${galdataname}_F606W_F814W.gst") {
  $photfile = "$galdataname/proc_new_dol/${galdataname}_F606W_F814W.gst" ;
  $filter = "F606W F814W";
}
elsif (-e "$galdataname/proc/${galdataname}_F606W_F814W.gst") {
  $photfile = "$galdataname/proc/${galdataname}_F606W_F814W.gst" ;
  $filter = "F606W F814W";
}
elsif (-e "$galdataname/proc_new_dol/${galdataname}_F475W_F814W.gst") {
  $photfile = "$galdataname/proc_new_dol/${galdataname}_F475W_F814W.gst" ;
  $filter = "F475W F814W";
}
elsif (-e "$galdataname/proc/${galdataname}_F475W_F814W.gst") {
  $photfile = "$galdataname/proc/${galdataname}_F475W_F814W.gst" ;
  $filter = "F475W F814W";
}
$filedir = dirname $photfile;
#previous $sfhfile = $filedir . "/sfh_full_grid.dat" ;
$astfile = $photfile ;
$astfile =~ s/.gst/_gst/ ;
$pezzo = basename $astfile ;
###############
$type = 'noAGB';
$zinc = 'yes';
$cmd_inputfile = 'cmd_input_ma08.dat' ;
#$cmd_inputfile = 'cmd_input_ma08.dat' ;
###############
if ($type eq 'noAGB') {
  if ($zinc eq 'no') {
    $sfhfile = "noAGB/".$pezzo.".sfh.dat";
  } else {
    $sfhfile = "noAGB/".$pezzo.".zinc.sfh.dat";
  }
#9884_M81F6D1_F606W_F814W_gst.sfh.dat
} else {
  $pezzo =~ s/\_gst/.gst.fits/ ; #eliminates initial number and _
  if ($zinc eq 'no') {
    $sfhfile = "withAGB/".$pezzo.".sfh.dat";
  } else {
    $sfhfile = "withAGB/".$pezzo.".zinc.sfh.dat";
  }
#9884_M81K61_F606W_F814W.gst.fits.sfh.dat
#9884_M81K61_F606W_F814W.gst.fits.zinc.sfh.dat
}
$astfile = $astfile.".matchfake";
$filedir = dirname $photfile;
#############################################
#print "pezzo is found to be $pezzo\n";
print "galname is found to be $galname\n";
print "photfile is found to be $photfile\n";
print "astfile is found to be $astfile\n";
print "sfhfile is found to be $sfhfile\n";
print "filters are found to be $filter\n";
#die;

open (IN, "<low_sfr_trgb.dat");
while (<IN>) {
  chomp; ($n1, $av, $trgb, $errtrgb, $dm0, $mtot) = split;
  if (index($n1,$galname)>=0) { 
    print "Found $n1: Av=$av, TRGB=$trgb+-$errtrgb, (m-M)0=$dm0\n";
    last;
  }
}
close IN ;
if ($trgb>10.0 && $trgb<28.0) {
  $camera="ACS";
}
if ($camera eq "ACS") {
  $photsysname="acs_wfc"; 
  $maglimfilter=10;
  $maglimnumber=29;
} elsif ($camera eq "WFPC2") {
  $photsysname="wfpc2vega"; 
  $maglimfilter=11;
  $maglimnumber=27;
}

###################################
#then get a distance
if ($trgb>10.0 && $trgb<28.0) {
  $mTRGB=$trgb;
  $dist = 10.0**(0.2*$dm0+1.0);
  $distMpc = $dist/1.0e6;
  print "Put $galname at d=$distMpc Mpc, Av=$av, mTRGB=$mTRGB\n";
}

&makesfhfile2;

##########################sets dirs#########################
chomp($homedir = `echo \$HOME`);
chomp($basedir = `pwd`); # contains perl macros and templates
$trilegaldir = "$homedir/eso/trilegal_1.3";

##########################counts stars in data #########################
print "Counting stars in $photfile\n";
$mTRGB2 = $mTRGB+$magbelowtip;
$mTRGB0 = $mTRGB-$magabovetip;
#$awkpattern = "BEGIN {rgb=0;agb=0}; NR>78 \$13\>$mTRGB && \$13\<$mTRGB2 {rgb=rgb+1}; NR>78 && \$13\>$mTRGB0 && \$13\<$mTRGB {agb=agb+1}; END {print rgb, agb}";
if  ($trgb>10.0 && $trgb<28.0) {
  $awkpattern = "BEGIN {rgb=0;agb=0}; \$13\>$mTRGB && \$13\<$mTRGB2 && \$5\<99 && \$13\<99 {rgb=rgb+1}; \$13\>$mTRGB0 && \$13\<$mTRGB && \$5\<99 && \$13\<99 {agb=agb+1}; END {print rgb, agb}";
} else {
  $awkpattern = "BEGIN {rgb=0;agb=0}; \$2\>$mTRGB && \$2\<$mTRGB2 && \$1\<99 && \$2\<99 {rgb=rgb+1}; \$2\>$mTRGB0 && \$2\<$mTRGB && \$1\<99 && \$2\<99 {agb=agb+1}; END {print rgb, agb}";
}
#print "awk \'$awkpattern\' <$photfile\n";
$_ = `awk \'$awkpattern\' <$photfile\n`;
($Nrgb, $Nagb) = split;
$ratio = $Nagb/$Nrgb ;
printf "Ratio = %.3f\n", $ratio;
$Nrgbinf = $Nrgb-int($nsigma*sqrt($Nrgb));
$Nrgbsup = $Nrgb+int($nsigma*sqrt($Nrgb));
print "Nrgb = $Nrgb\nNagb = $Nagb\n";

##########################runs trilegal#########################
#$mtot = 10000000;
#&trilegal4LMC; # runs TRILEGAL for LMC
while (1) {
  &trilegal4LMC; # runs TRILEGAL for LMC
  #################counts stars in model###############
  $awkpattern = "BEGIN {rgb=0;agb=0}; NR>1 && \$12+\$14\>$mTRGB && \$12+\$14\<$mTRGB2 && \$13\<99 && \$14\<99 {rgb=rgb+1}; NR>1 && \$12+\$14\>$mTRGB0 && \$12+\$14\<$mTRGB && \$13\<99 && \$14\<99 {agb=agb+1}; END {print rgb, agb}";
  #print "awk \'$awkpattern\' <$modfile\n";
  #$_ = `awk \'$awkpattern\' <$modfile\n`;
  $_ = `awk \'$awkpattern\' <$astmodfile\n`;
  ($Mrgb, $Magb) = split;
  print "Mrgb = $Mrgb\nMagb = $Magb\n";
  $ratio = $Magb/$Mrgb ;
  printf "Ratio = %.3f\n", $ratio;
  printf ("Modelled %.0f RGBs as compared to (%.0f,%.0f,%.0f)\n",
	  $Mrgb, $Nrgbinf, $Nrgb, $Nrgbsup );
  $fact_RGB = $Nrgb/$Mrgb;
  if ($Mrgb>$Nrgbinf && $Mrgb<$Nrgbsup) {
    print "that's within $nsigma sigma\n";
    printf("Mtot= %.3e\n",$mtot);
    last; #leaves if $N_RGB is within 1sigma of observed value
  } 
  $mtot = $mtot * $fact_RGB ; # else change mass and do it again
}

##########################counts stars in model#########################
$awkpattern = "BEGIN {rgb=0;agb=0}; NR>78 && \$21\>$mTRGB && \$21\<$mTRGB2 {rgb=rgb+1}; NR>78 && \$21\>$mTRGB0 && \$21\<$mTRGB {agb=agb+1}; END {print rgb, agb}";
print "awk \'$awkpattern\' <$modfile\n";
$_ = `awk \'$awkpattern\' <$modfile\n`;
($Mrgb, $Magb) = split;
print "Mrgb = $Mrgb\nMagb = $Magb\n";
$ratio = $Magb/$Mrgb ;
printf("Ratio = %.3f\n", $ratio);

########################makes final plot and statistics###################
$plotfile = $galname.".ps";
$plotfile2 = $galname.".eps";
open (SM, "|sm2 >lixo");
print SM "device postencap tmp.ps\n";
if ($filter eq "F606W F814W") {
  print SM "macro read cmd.sm cmd $photfile $astmodfile $galname $mTRGB 606\n";
  print "sm: macro read cmd.sm cmd $photfile $astmodfile $galname $mTRGB 606\n";
} elsif ($filter eq "F475W F814W") {
  print SM "macro read cmd.sm cmd $photfile $astmodfile $galname $mTRGB 475\n";
  print "sm: macro read cmd.sm cmd $photfile $astmodfile $galname $mTRGB 475\n";
}
print SM "hardcopy\n";
print SM "quit\n";
#`gv $plotfile \&`;
close SM;
`mv tmp.ps $plotfile`;
`ps1eps $galname`;

print "******* FINISHED, see $plotfile ************\n\n";
print "mv input*.dat out_*.dat output_*.dat *.ps *.eps fig_ma08????";
die;

##########################trilegal4LMC#####################################
sub  trilegal4LMC {
  ####### runs TRILEGAL 
  $inpfile = "$basedir/input_$galname.dat";
  $modfile = "$basedir/output_$galname.dat";
  $astmodfile = "$basedir/out_$galname.dat";
  $astmodfile = basename($astmodfile);
  printf("Preparing input TRILEGAL file at $inpfile\n");
  $tmp = $basedir; $tmp =~ s#\/#\\/#g; #print "$basedir $workdir \n";
  `sed  -e 's/PHOTSYS/$photsysname/' -e 's/MAGLIMFILTER MAGLIMNUMBER/$maglimfilter $maglimnumber/' -e 's/EXTINCTION/0.0/' -e 's/MTOT DIST/$mtot $dist/' -e 's/AV/$av/' -e 's/IMAGEDIR/$tmp/' -e 's/SFHFILE/$trilegalsfhfile/'  -e 's/FACA FACB/$faca $facb/' < $basedir/input_template.dat >$inpfile`; #prepares input file
  #runs trilegal:
  chdir("$trilegaldir");
  if ($runitreally==1) {
    $message = `code_1.3.1/main -a -f $cmd_inputfile $inpfile $modfile`;
  }
  print "code/main -a -f cmd_input_ma08.dat $inpfile $modfile\n";
  $message .= "$modfile has size " . `wc $modfile` . "\n";
  if ($debug==1) {print "TRILEGAL run for galaxy:\n $message\n";}
  chdir("$basedir");
  ####### applies ASTs
  $astmodfile = "$basedir/out_$galname.dat";
  $astmodfile = basename($astmodfile);
  if ($runitreally==1) {
    open(ASTCODE,"|code/spread_angst >lixo");
    print ASTCODE "$astfile $modfile $astmodfile\n";
    close ASTCODE;
  }
}

###################makes SFH file##########################
sub makesfhfile1
  {
    #log(age)_low
    #log(ag)_high
    #dmod_high_best_or_low
    #SFR
    #SFR_err
    #SFR_err
    #<[M/H]>
    #<[M/H]>_err
    #<[M/H]>_err
    #[M/H]_spread
    #[M/H]_spread_err
    #[M/H]_spread_err
    #Cumulative_Stellar_Mass
    #Cumulative_Stellar_Mass_err_high
    #Cumulative_Stellar_Mass_err_low
    $trilegalsfhfile = "sfh_".$galname.".dat";
    print "Writing SFH from $sfhfile\n to $trilegalsfhfile\n";
    
    open (SFHI,"<$sfhfile");
    open (SFHO,">$trilegalsfhfile");
    $sumcumsfr = 0.0;
    while (<SFHI>) {
      chomp; 
      if (index($_,"Best fit:")==0) {
	$_ = <SFHI>;
	($av,undef,$dmod) = split;
	$av =~ s/Av=//g ; #eliminates $
	$av =~ s/\+\/\-\d+.\d+,//g ; #eliminates errorbar
	print "Reset Av to $av\n"; 
      } else {
	($lage1,$lage2,$dm,$sfr,undef,undef,$mh,undef,undef,$sigmh) = split;
	if ($lage1>6.5 && $lage1<10.3) {
	  $age1 = 1.0e-9*10**$lage1;
	  $age2 = 1.0e-9*10**($lage2-0.01*($lage2-$lage1));
	  $sfr = $sfr;
	  $zeta = 0.019*10.0**$mh ;
	  $sigmh = ($sigmh>0.0) ? $sigmh : 0.01 ;
	  $cumsfr = $sfr*($age2-$age1)/0.366391312651252;
	  $sumcumsfr = $sumcumsfr+$cumsfr;
	  $insumcumsfr = 1.0-$sumcumsfr;
	  #    print SFHO "$age1 $sfr $zeta 0.01\n$age2 $sfr $zeta 0.01\n";
	  print SFHO "$age1 $sfr $zeta $sigmh\n$age2 $sfr $zeta $sigmh\n";
	  printf("%.3g %.3g %.3g %.3g %.3g %.3g\n", 
		 $age1, $age2, $sfr, $zeta, $sigmh, $cumsfr, $insumcumsfr);
	}
      }
    }
    $age2 = $age2*1.001;
    print SFHO "$age2 0.0 $zeta $sigmh\n";
    close SFHI;
    close SFHO;
    $faca=1.0e9; $facb=0.0;
  }

############################################
sub makesfhfile2
  {
    $trilegalsfhfile = "sfh_".$galname.".dat";
    print "Writing SFH from $sfhfile\n to $trilegalsfhfile\n";
    
    $sfh_young = 0.0; $sfh_int = 0.0; $sfh_old = 0.0;
    $fehmass_sum = 0.0; $mass_sum = 0.0;
    open (SFHI,"<$sfhfile");
    open (SFHO,">$trilegalsfhfile");
    while (<SFHI>) {
      chomp;
      ($feh1,$feh2,$logt1,$logt2,$sfr) = split; 
      next if ($sfr==0.0) ;
      $age1a=1.0e-9*10**$logt1;
      $age1p=1.0e-9*10**($logt1+0.0001);
      $age2a=1.0e-9*10**$logt2;
      $age2p=1.0e-9*10**($logt2+0.0001);
      $zeta1=0.019*10**$feh1;
      $zeta2=0.019*10**$feh2;
      $sfr=$sfr*1.0e6;
      if ($logt2<=9.0) {
	$sfh_young = $sfh_young + $sfr*($age2a-$age1a);	
      } elsif ($logt2<=9.477){
	$sfh_int = $sfh_int + $sfr*($age2a-$age1a);
      } else {
	$sfh_old = $sfh_old + $sfr*($age2a-$age1a);
      }
      $fehmass_sum = $fehmass_sum + 0.5*($feh1+$feh2)*$sfr*($age2a-$age1a);
      $mass_sum = $mass_sum + $sfr*($age2a-$age1a);
      print SFHO "$age1a 0.0 $zeta1\n";
      print SFHO "$age1p $sfr $zeta1\n";
      print SFHO "$age2a $sfr $zeta2\n";
      print SFHO "$age2p 0.0 $zeta2\n";
      print SFHO "$age1a 0.0 $zeta2\n";
      print SFHO "$age1p $sfr $zeta2\n";
      print SFHO "$age2a $sfr $zeta1\n";
      print SFHO "$age2p 0.0 $zeta1\n";
    }
    close SFHI;
    close SFHO;
    $sfh_tot = $sfh_young + $sfh_int + $sfh_old ;
    $feh_mean = $fehmass_sum / $mass_sum ;
    printf "%s SFH young= %.2f, int= %.2f old= %.2f, [Fe/H]mean= %.2f\n",
      $galname, $sfh_young/$sfh_tot, $sfh_int/$sfh_tot, $sfh_old/$sfh_tot, 
	$feh_mean ; 
    $faca=1.0e9; $facb=0.0;
  }
