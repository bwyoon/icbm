 #!/usr/bin/perl

use Math::Trig;

use constant {
    G => 6.67408E-11,
    M => 5.972E+24,
    R => 6.371E+6
};

if (scalar(@ARGV) < 2) {
    print STDERR "USAGE: perl icbm.pl Vi Angle\n";
    exit(1);
}

my ($vi, $ai) = @ARGV;
$vi *= 1000.0;
my @v = ($vi*cos($ai*pi/180.0), $vi*sin($ai*pi/180.0));
my $dt = 1;
my @x = (0.0, R);
my $r =$rmax = $rprev = R;
my @a;
my $t = 0.0;

# initial acceleration
for (my $k = 0; $k < 2; $k++) {
    $r = sqrt($x[0]*$x[0]+$x[1]*$x[1]);
    $a[$k] = -G*M*$x[$k]/$r/$r/$r;
}

# time integration
my $n = 0;
while (1) {
    #print "$n : @x\n";
    $n++;
    # velocity verlet
    for ($k = 0; $k < 2; $k++) {
        $v[$k] += $a[$k]*0.5*$dt;
        $x[$k] += $v[$k]*$dt;
    }

    $rprev = $r;
    $r = sqrt($x[0]*$x[0]+$x[1]*$x[1]);
    $t += $dt;
    last if ($r < R);

    for ($k = 0; $k < 2; $k++) {
        $a[$k] = -G*M*$x[$k]/$r/$r/$r;
        $v[$k] += $a[$k]*0.5*$dt;
    }

    $rmax = $r if ($r > $rmax);
}

# correction of the last position
for ($k = 0; $k < 2; $k++) {
    $x[$k] -= $v[$k]*$dt;
}
$t -= $dt;

$dt *= (R-$rprev)/($r-$rprev);

for ($k = 0; $k < 2; $k++) {
    $x[$k] += $v[$k]*$dt;
}
$r = sqrt($x[0]*$x[0]+$x[1]*$x[1]);
$t += $dt;

print "$n : @x\n";

# maximum altitude
my $hmax = $rmax - R;

# maximum distance
my $angle = acos($x[1]/R);
my $dist = R*$angle;

print "max altitude = ", $hmax/1000.0, " km\n";
print "max range = ", $dist/1000.0, " km\n";