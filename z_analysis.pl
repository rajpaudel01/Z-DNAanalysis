#!/usr/bin/perl
#Rajan Paudel 
#pur-pyr
$chromo = 1;$pairs ='';
for $ko(0..9){
	%{h.$ko} = 0; 
	%{he.$ko} = 0;
	${tot.$ko} = 0;
	${ek.$ko} = 0;
}
for $chromo(1..22){
	%ekpos = 0;
	open(EK,"/home/rajan/the11/the1101data/z_80/z_r_$chromo")||die 'cannot open';
	while(<EK>){
		$_ =~ /chr\d+\t\d+\t(\d+)-(\d+)/;
		for $con($1..$2){
			$ekpos{$con}=1;
			}
	}
	close EK;
	open(CHR,"/home/rajan/the11/the1101data/chro/chr$chromo")|| die 'cannot open it';
	while(<CHR>){
		chomp $_;
		@a = split ("",$_);
	}
	close CHR;
	open(CHECK,">>check")|| die 'cannot create check';
	print CHECK "chromosone$chromo\n";
	open(IN,"/home/rajan/gthou/1000G_chr$chromo")||die 'cannot find the file';
	while(<IN>){
		#last if $. > 400000;
		if  ($_ =~ /^\d+\t(\d+)\t(\w)\t(\w)\t(\d\.?\d*)\t(\w)/){
			$anc = $5;
			if ($anc eq $2){$der = $3;$dfr = $4;}
			elsif($anc eq $3){$der = $2;$dfr = 1-$4;}
			else{next;}
			$bef = $1-2;
			#$ano = $1-1;
			$pairs = "$a[$bef]"."$anc"."$a[$1]"."\t"."$a[$bef]"."$der"."$a[$1]";
			#print"$pairs\n";
			#print"$a[$bef].$a[$ano].$2.$a[$1]\n";
			$pairs =~ tr/aaAAggGGccCCttTT/UUUUUUUUYYYYYYYY/;
			#print "$pairs\n";
			$tot0++;$h0{$pairs}++;
			if($dfr > 0 && $dfr <= 0.002){$tot1++;$h1{$pairs}++;}else{;}
			if($dfr > 0.002 && $dfr <= 0.10){$tot2++;$h2{$pairs}++;}else{;}
			if($dfr > 0 && $dfr <= 0.20){$tot3++;$h3{$pairs}++;}else{;}
			if($dfr > 0.20 && $dfr <= 0.40){$tot4++;$h4{$pairs}++;}else{;}
			if($dfr > 0.40 && $dfr <= 0.60){$tot5++;$h5{$pairs}++;}else{;}
			if($dfr > 0.60 && $dfr <= 0.80){$tot6++;$h6{$pairs}++;}else{;}
			if($dfr > 0.80 && $dfr <= 1){$tot7++;$h7{$pairs}++;}else{;}
			if($dfr > 0.90 && $dfr <= 0.998){$tot8++;$h8{$pairs}++;}else{;}
			if($dfr > 0.998 && $dfr <= 1){$tot9++;$h9{$pairs}++;}else{;}
			if ($ekpos{$1}){
				$ek0++;$he0{$pairs}++;
	                        if($dfr > 0 && $dfr <= 0.002){$ek1++;$he1{$pairs}++;}else{;}
	                        if($dfr > 0.002 && $dfr <= 0.10){$ek2++;$he2{$pairs}++;}else{;}
        	                if($dfr > 0 && $dfr <= 0.20){$ek3++;$he3{$pairs}++;}else{;}
                	        if($dfr > 0.20 && $dfr <= 0.40){$ek4++;$he4{$pairs}++;}else{;}
                        	if($dfr > 0.40 && $dfr <= 0.60){$ek5++;$he5{$pairs}++;}else{;}
                       		if($dfr > 0.60 && $dfr <= 0.80){$ek6++;$he6{$pairs}++;}else{;}
                       		if($dfr > 0.80 && $dfr <= 1){$ek7++;$he7{$pairs}++;}else{;}
                        	if($dfr > 0.90 && $dfr <= 0.998){$ek8++;$he8{$pairs}++;}else{;}
                        	if($dfr > 0.998 && $dfr <= 1){$ek9++;$he9{$pairs}++;}else{;}
			}else{;}

		}else{print"no match\n";}
		print "$.\n";
	}
$. = 0;		
}
open(ALL,">fulltable")|| die 'cannot create all';
$" = "\t";
foreach $am(keys %h0){
		push (@{a.$am},$am);
		for $ti(0..9){
        		push(@{a.$am},${h.$ti}{$am});
			}
		for $ti(0..9){
			push(@{a.$am},${he.$ti}{$am});
			}
	print ALL "@{a.$am}\n";
	}
open(DAT,"fulltable")|| die 'cannot open file';
open(OUT,">finaltable")|| die "cannot create"; 
while(<DAT>){
	chomp $_;
       	if  ($_ =~ /(\w)(\w)(\w)\t(\w)(\w)(\w)\t(.+)/){
		if ($2 ne $5){
			@ray = split("\t",$7);
			if ($1 eq $2 && $2 eq $3){
				for $yo(0..19){
 				${a.$yo} +=  $ray[$yo];
				$incre += $ray[$yo];
				}
			}elsif($4 eq $5 && $5 eq $6){
                		for $yo(0..19){
 				${b.$yo} +=  $ray[$yo]; 
				$decre += $ray[$yo];
				}                      
                        }else{;}
 		}else{;}
	}else{;}
}

print "@ray\n";
foreach $ya(0..19){
	push(@arr,${a.$ya});
	push(@brr,${b.$ya});
	$rati = $arr[$ya]/$brr[$ya];
	$rati = sprintf("%.2f", $rati);
	push (@crr,$rati);
}
print"@arr\n";
print"@brr\n";
$" = "\t";
print OUT "increasing z content >> $incre\n decreasing z content >> $decre\n\n$mutation\tall\t(<0.2%)\t(>0.2% & <=10%)\t(>0%&<=20%)\t(>20%&<=40%)\t(>40%&<=60%)\t(>60%&<=80%)\t(>80%&<=100%)\t(>90%&<=99.8%)\t(>99.8%&<=100%)\tz_all\tz_(<0.2%)\tz_(0.2% & <= 10%)\tz_(>0%&<=20%)\tz_(>20%&<=40%)\tz_(>40%&<=60%)\tz_(>60%&<=80%)\tz_(>80%&<=100%)\tz_(>90%&<=99.8%)\tz_(>99.8%&<=100%)\nnonz>z\t@arr\nz> non-z\t@brr\n\t@crr\n";




	
