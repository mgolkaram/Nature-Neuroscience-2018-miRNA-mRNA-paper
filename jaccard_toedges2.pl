use warnings; 
use strict; 

open (NEIGHBORS, "<", $ARGV[0]); 

my %cellkey; 
my %dist; 
my %neigh; 

my $head = 0; 
while (<NEIGHBORS>){
    chomp; 
    if ($head == 0){
	$head++; 
	next; 
    }
    my @col = split /\t/; 
    $cellkey{$head} = $col[0]; 
    my @neigh = @col[2..$#col]; 
    $neigh{$head} = \@neigh; 
    $head++; 
}

close (NEIGHBORS); 

open (DISTANCE, "<", $ARGV[1]);

my $count = 0; 

while (<DISTANCE>){
    chomp;
    if ($count == 0){
	$count++; 
	next; 
    }
    my @col = split /\t/;
    my @dist = @col[2..$#col];
    $dist{$count} = \@dist;
    $count++; 
}

close (DISTANCE); 

my %edges; 

my @cells = keys (%cellkey); 

print "V1\tV2\tweight\n"; 

foreach my $c (@cells){
    my @neigh = @{$neigh{$c}}; 
    my %check; 
    foreach my $n (@neigh){
	$check{$n} = 1; 
    }
    my @dist = @{$dist{$c}}; 
    my $n = 0; 
    for ($n = 0; $n <= $#neigh; $n++){
	my $same = 0; 
	my $diff = 0; 
	my $edge = join ("_", $c, $neigh[$n]); 
	if (exists $edges{$edge}){
	    next; 
	}
	else{
	    $edges{$edge} = 1; 
	}
#	print "$c\t$n\t$neigh[$n]\n"; 
	my @neighbors = @{$neigh{$neigh[$n]}};
	foreach my $y (@neighbors){
	    if (exists $check{$y}){
		$same++; 
	    }
	    else{
		$diff++; 
	    }
	}
	my $jaccard = $same/(($#neigh+ 1 - $same) + $same + $diff); 
	my $weighted = $jaccard*$dist[$n]; 
	my $edge1 = $cellkey{$c}; 
	my $edge2 = $cellkey{$neigh[$n]}; 
	print "$edge1\t$edge2\t$weighted\n"; 
    }
}
