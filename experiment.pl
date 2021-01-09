#!/usr/bin/env perl

$experiment_name = "AnyAngle";
@algorithms = (
	"EXP_A_EUC", 
	"EXP_A_OCT", 
	"EXP_T", 
	"EXP_L", 
	"EXP_F", 
	"EXP_B", 
	"EXP_SUB_1_A", 
	"EXP_SUB_1_T", 
	"EXP_SUB_2_A", 
	"EXP_SUB_2_T", 
	"EXP_SUB_10000_A", 
	"EXP_SUB_10000_T", 
	"EXP_ANYA"
	);

# Read all the map names from the current directory
# Only look at .map.scen files
@maps = ();
opendir(DIR,".") || "Can't open dir";
@FILES= readdir(DIR);
foreach my $filename (@FILES) {
	if($filename =~ /^(.*)\.map\.scen$/){
		$mapname = $1;
		push(@maps, $mapname);
	}
}
closedir(DIR);

# Write the maplist 
open(OUT, ">maplist");

foreach my $map (@maps) {
	print OUT "$map \n";
}
close(OUT);

# Create the results folder
$results = "results/$experiment_name/";
$instances = "results/$experiment_name/instances/";
#system "rm -rf $results";
system "mkdir results/";
system "mkdir $results";
system "mkdir $instances";

# Create a temporary folder for storing processed scenario files
# This is done so that if a crash occurs, the experiment can resume
# from where it left off
system "mkdir temp/";

# For each map, run each algorithm
foreach my $map (@maps) {
	foreach my $algorithm (@algorithms){
		system "./$algorithm $map\n";
		system "mv Instances-* $instances";
		system "mv *-$map $results";
	}
	
	# Once the experiment is done, move the associated scenario
	# file to a temporary directory
	
	system "mv $map.map.scen temp/";
}

# Once the experiment is done, move all summary files and maplist
system "mv maplist $results/";

# move all the scen files back to their original location
system "mv temp/* .";
system "rm -rf temp/";
