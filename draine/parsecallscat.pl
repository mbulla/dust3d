#!/usr/bin/perl -w
#
# ABOUT
#   Script to parse the callscat_init_p* files from Bruce Draine into
#   one matrix for each dust type, which will be much easier to read
#   from a low level language.
#
#   Three files are produced for each "dust type", one containing the
#   table, called '.dat', and two files describing the angle
#   ('_angle.dat') and wavelength ('_lambda.dat') vectors.
#
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::Graphics::LUT;
use PDL::Graphics::PGPLOT::Window;

my $base = "callscat_init_p";

my @lambrange = ('uv','','ir');
my @dusttypes = ('LMC_avg','SMC_bar','MW_3.1');

foreach my $type (@dusttypes) {                  # Loop over dust type
  my $angles = undef;
  my (@lambda,@table);

  foreach my $range (@lambrange) {               # Loop over lambda range
    my $ifile = "${base}.out_${type}";
    $ifile .= "_${range}" unless $range eq '';

    my $reverse = 0;                 # Should the lambda vector be reversed
    
    # Retrieve lambda vector
    #
    my @foo = split /\=/, `grep wav ${ifile}`; #/
    my @thislambda = split /\s+/, $foo[1];
    shift @thislambda; shift @thislambda;

    my @thistable = rcols $ifile;
    if (defined $angles) {
      shift @thistable;
    } else {
      $angles = shift @thistable;
    }

    # Check if the vector is descending or ascending
    #
    my $thislambda = pdl \@thislambda;
    my $difflambda = $thislambda(1:) - $thislambda(:-2);
    $reverse = 1 if $difflambda->where($difflambda < 0)->nelem > 0;
    $difflambda *= -1 if $reverse;
    
    # Make sure the vector is sorted
    #
    die "Error: lambda vector does not appear to be sorted for '${type}_${range}'!"
      unless $difflambda->where($difflambda < 0)->nelem == 0;

    if ($reverse) {
      print "REVERSING VECTORS FOR: $type $range\n";
      @thislambda = reverse @thislambda;
      @thistable = reverse @thistable;
    }

    push @lambda, @thislambda;
    push @table, @thistable;

    $thislambda = pdl \@thislambda;
    my $thistable = pdl \@thistable;
    
    # Diagnostics plots that can be compared to Fig. 5 in
    # Draine et al. (2003)
    #
    my $w = pgwin(Dev=>"${base}_${type}_${range}.png/png");
    $w->env($angles->min,$angles->max,$thistable->min,$thistable->max,
	    {XTitle=>"angle (deg)",YTitle=>"P",Title=>"$type $range"});
    my (@legend,@colour);
    for (my $i=0; $i<$thislambda->nelem; $i++) {
      $w->line($angles,$thistable(:,($i)),{COLOUR=>$i+1});
      push @legend, $thislambda($i)->sclr;
      push @colour, $i+1;
    }
    $w->legend(\@legend,5,0.95*$thistable->max,{COLOUR=>\@colour,CharSize=>0.7});
    $w->release;
    $w->close;
  }

  my $lambda = pdl \@lambda;
  my $m = pdl \@table;
  $m = $m->transpose;
  my ($nx,$ny) = $m->dims;

  # Print the main table
  #
  open FP, ">${base}_${type}.dat" ||
    die "Error: failed to open file handle!";
  print FP "# Each row corresonponds to one angle in the angle vector\n".
    "# Each column corresponds to one wavelength in the lambda vctor\n";
  for (my $jj=0; $jj<$ny; $jj++) {
    for (my $ii=0; $ii<$nx; $ii++) {
      printf FP "%7.4f ", $m($ii,$jj)->sclr;
    }
    print FP "\n";
  }
  close FP;

  # Print the angle and lambda vectors
  #
  wcols "%8.4f", $angles, "${base}_${type}_angle.dat",
    { HEADER => "# Scattering angle in degrees" };
  wcols "%8.4f", $lambda, "${base}_${type}_lambda.dat",
    { HEADER => "# Wavelength in microns" };

  
  # Diagnostics plot
  #
  # This plot will not reflect reality since the PGPLOT does not support
  # surface plots with non-equidistant steps.
  #
#  $lambda *= 10000;
#  my $w = pgwin(Dev=>"${base}_${type}.png/png");
#
#  my ( $l, $r, $g, $b ) = lut_data('bgyrw');
#  $w->ctab($l,$r,$g,$b);
#
#  my $t = pdl [$lambda->min,($lambda->max-$lambda->min)/($lambda->nelem-1),0,
#	       $angles->min,0,($angles->max-$angles->min)/($angles->nelem-1)];
#  $w->env($lambda->min,$lambda->max,$angles->min,$angles->max,
#	  {CharSize=>0.7,XTitle=>'\gl (\A)',YTitle=>'angle (deg)',Title=>$type});
#  $w->imag($m,$m->min,$m->max,{TRANSFORM=>$t,DrawWedge=>1,ITF=>'SQRT'});
#  foreach my $minmax (@lambdaminmax) {
#    my $foo = pdl [$minmax,$minmax];
#    my $bar = pdl [$angles->min,$angles->max];
#    $w->lines(10000*$foo,$bar,{LineWidth=>3,LineStyle=>'solid',Colour=>'white'});
#  }
#  $w->release;
#  $w->close;

}

# EOF
#
