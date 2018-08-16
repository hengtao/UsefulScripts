#!/usr/bin/perl
#Name: Transcript Similarity Calc
#Author: Chaoyong Xie
my $newfile = $ARGV[0];
my $lastfile = $ARGV[1];
if(defined($newfile) && defined($lastfile)){
	open(NEWFH,$newfile) or die "Cannot read $newfile:$!";
	open(LASTFH,$lastfile) or die "Cannot read $lastfile:$!";
	my @records_new = ; chomp @records_new;
	my @records_last = ; chomp @records_last;
	#print "Both files Loaded!\n";
	close(NEWFH) or die "Cannot close $newfile:$!";
	close(LASTFH) or die "Cannot close $lastfile:$!";
	my @fields_new = ();
 	my @fields_last = ();
	foreach $record_new (@records_new){
		@fields_new = split(/\t/,$record_new);
		foreach $record_last (@records_last){
			@fields_last = split(/\t/,$record_last);
			if($fields_new[0] eq  $fields_last[0] && $fields_new[5] eq $fields_last[5]){
				my @offset_new = split(/,/,$fields_new[11]);
				my @length_new = split(/,/,$fields_new[10]);
				my @offset_last = split(/,/,$fields_last[11]);
				my @length_last = split(/,/,$fields_last[10]);
				my $len_sum = 0;
				for(my $i = 0; $i < $fields_new[9]; $i ++){
					for(my $j = 0; $j < $fields_last[9]; $j ++){
					if($fields_last[1] + $offset_last[$j] > $fields_new[1] + $offset_new[$i] + $length_new[$i]){
						next;
					}else{
						if($fields_last[1] + $offset_last[$j] + $length_last[$j] < $fields_new[1] + $offset_new[$i]){
							next;
						}else{
							if($fields_last[1] + $offset_last[$j] < $fields_new[1] + $offset_new[$i] ){
								if($fields_last[1] + $offset_last[$j] + $length_last[$j] < $fields_new[1] + $offset_new[$i] + $length_new[$i]){
									$len_overlap += $fields_last[1] + $offset_last[$j] + $length_last[$j] - $fields_new[1] -$offset_new[$i];
								}else{
									$len_overlap += $length_new[$i];
								}
							}else{
                     			if($fields_last[1] + $offset_last[$j] + $length_last[$j] < $fields_new[1] + $offset_new[$i] + $length_new[$i]){
									$len_overlap += $length_last[$j];
								}else{
									$len_overlap += $fields_new[1] + $offset_new[$i] + $length_new[$i] - $fields_last[1] - $offset_last[$j];
								}
							}
						}
					}
				}
			}
        	if($len_overlap > 0){
				$len_sum = 0;
				for(my $t = 0; $t < $fields_new[9]; $t++){
					$len_sum += $length_new[$t];
				}
				my $similarity = 0;
				$similarity = $len_overlap / ($len_sum + 0.0);
				print "$fields_new[0]\t$fields_new[1]\t$fields_new[2]\t$fields_new[3]\t",substr($lastfile,-13,9),"\t$fields_last[3]\t$similarity\n";
        	}
      		}else{
         		next;
			}
		}
	}
}else{
	die "There must be 2 parameters!\n";
}
