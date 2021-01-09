#!/usr/bin/env perl

@experiments = (
"bg512",
"DAO",
"starcraft",
#"wc3maps512",
"random10",
"random15",
"random20",
"random25",
"random30",
"random35",
"random40",
"game",
"random",
#"room8",
#"room16",
#"room32",
#"room64",
#"maze1",
#"maze2",
#"maze4",
#"maze8",
#"maze16",
#"maze32"
);

@algorithms = (
"A_EUC",
"A_OCT", 
"T", 
"L",
"F",
"B",
"S_1_A", 
"S_2_A", 
"S_10000_A",
"S_1_T", 
"S_2_T", 
"S_10000_T",
"ANYA"
);

@ps_algorithms = (
"A_EUC_PS",
"A_OCT_PS", 
"T_PS", 
"L_PS",
"F_PS",
"B_PS",
"S_1_A_PS", 
"S_2_A_PS", 
"S_10000_A_PS",
"S_1_T_PS", 
"S_2_T_PS", 
"S_10000_T_PS",
"ANYA"
);

@all_algorithms = (@algorithms, @ps_algorithms);

sub getMaps
{
	my $algorithm = $_[0];
	my @maps = ();
	opendir(DIR,".") || "Can't open dir";
	@FILES= readdir(DIR);
	foreach my $filename (@FILES) {
		chomp $filename;
		if($filename =~ /^$algorithm-(.*)$/){
			$mapname = $1;
			if ($mapname ne "summary")
			{
				push(@maps, $mapname);
			}
		}
	}

	closedir(DIR);
	return @maps;
}

sub getVal
{
	my $algorithm = $_[0];
	my $mapname = $_[1];
	my $description = $_[2];
	
	open(IN, "$algorithm-$mapname");
	foreach $line (<IN>)  {
		chomp $line;
		if ($line =~ /$description.*:\s*(\S*)\s*$/)
		{
			#print "$1\n";
			return $1;
		}
	}
}

sub getData
{
	my $algorithm = shift;
	my $description = shift;
	my $average = shift;
	
	my @maps = getMaps($algorithm);
	
	if ($#maps == -1)
	{
		return "-";
	}

	my $sum = 0;
	my $count = 0;
	
	foreach my $map (@maps) {
		my $val = getVal($algorithm, $map, $description);
		$count ++;
		$sum = $sum + $val;
	}
	
	if($average == 0) {
		$count = 1;
	}
	return $sum/$count;
}

sub printData
{
	my $algorithm = shift;
	my $description = shift;
	my $average = shift;

	my $val = getData($algorithm, $description, $average);
	#print "($algorithm) $new_description : $ave \n";
	print OUT "\t$val";
}

sub getMapRow
{
	my $experiment = $_[0];
	chdir($experiment);
	print OUT "$experiment";

	my $average = 0;
	
	printData("A_OCT", "Number of solved instances", $average);
	
	foreach (@algorithms){
	#	printData($_, "Number of solved instances", $average);
	}	
	
	foreach (@algorithms){
	#	printData($_, "Number of unsolved instances", $average);
	}		
	
	foreach (@algorithms){
	#	printData($_, "Valid searches - Total search time", $average);
	}	
	
	foreach (@algorithms){
	#	printData($_, "Valid searches - Total solution cost", $average);
	}
	
	foreach (@ps_algorithms){
	#	printData($_, "Valid searches - Total search time", $average);
	}	
	
	foreach (@ps_algorithms){
	#	printData($_, "Valid searches - Total solution cost", $average);
	}	
	
	foreach (@algorithms){
	#	printData($_, "Valid searches - Total number of expansions", $average);
	}
		
	foreach (@algorithms){
	#	printData($_, "Valid searches - Total number of LOS checks", $average);
	}
	
	foreach (@ps_algorithms){
	#	printData($_, "Valid searches - Total number of direction changes", $average);
	}	
	
	foreach (@ps_algorithms){
	#	printData($_, "Valid searches - Total number of freespace direction changes", $average);
	}	
	
	printData("S_2_T", "Time to construct SSG", $average);
	printData("S_2_T", "Time to partition the SSG", $average);
	
	print OUT "\n";
	chdir("..");
}

open (OUT, ">results.txt");
foreach $experiment (@experiments)
{
	getMapRow($experiment);
}
close (OUT);