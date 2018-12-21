#!/usr/bin/perl
$ch = 1;$p = 0;
for $ch(1..22){
	#system("tail -n +2 chr$ch.fa | paste -s -d \"\" >chr$ch");
	open(OUT,">z_r_80_$ch")||die 'cannot do open';
	open(IN,"/home/rajan/the1120/chr/chr$ch")|| die "Couldn't open file: $!";
	$ini = '';$tot = 0;$n=1;$se = '';$seq ='';
	while (<IN>) {
		$whole = $_;
		$leng = () = $_ =~ /[GgCcAaTtN]/g;
	}
	print "$leng\n";
	close IN;
	while(){last if ($n > $leng);
		$whole =~ /\G(.{100})/gc;$ini = $1;$n = $n+100;
		$nb = () = $ini =~ /[aAgG][cCtT]/g;$nc = () = $ini =~ /[cCtT][aAgG]/g;$na = $nb + $nc;
		if ($na > 80){$m = $n-100;
               		$seq = $ini;
			START:
			$whole =~ /\G(.{10})/gc;
        		$str = $1;$pseq = $seq;$pn = $n;$nall = $pn - $m;
			$seq = $seq.$str;
        		$nama = () = $seq =~ /[AaGg][cCtT]/g;$namb = () = $seq =~ /[cCtT][AaGg]/g;
			$nam = $nama + $namb;
			$n = $n + 10;$at = $nam/$nall;
        			if ($at > 0.80){goto START;
				}else{$eachreg = $pn - $m; $p += $eachreg;$ztot += $nam; print OUT "chr$ch\t$nall\t$m-$pn\t$pseq\n";}
		}else{$ztot += $na;}
	}
	$wtot +=$n;
	close OUT;
}
$seqocc = $p*100/$wtot;
print "percentage occupied by rich regions = $seqocc\n";
$perce = $ztot*100/$wtot;
print "percentage of total U and Y alternations = $perce\n";


