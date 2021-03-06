#!/bin/bash

#R. Berlemont - renoberlemont@gmail.com
#This Script requires the HMMer package (see http://hmmer.org)
#This Script requires Pfam-A.hmm... if you need it, run the folowinf command in the terminal to dowload the db in you current directory.
#wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
#gunzip  Pfam-A.hmm.gz
###########################################################
TS=$(date +"%y_%m_%d_%H%M")
clear
echo Loading the data\:
echo 1\-List of PFam is in\: 
read PFDom
if [[ ! -f $PFDom ]] ; 
then
    echo File  $PFDom is not there, aborting.
    exit
fi
echo "###################################################"
echo creating the custom DB and \"pressing\" Pfam\-A.hmm \(if required\)
echo "..."
gtr '\n' '@' < Pfam-A.hmm | sed 's/\/\//\/\/\n/g' | grep -F -f $PFDom |  tr '@' '\n'  >.db1
hmmpress Pfam-A.hmm
echo done
echo "###################################################"
echo 2\-Sequence file \(FASTA\)\:
read SEQFILE
	if [[ ! -f $SEQFILE ]] ; 
		then
    	echo File  $SEQFILE is not there, aborting.
    	exit
		else
		grep -v "^$" $SEQFILE | cut -d " " -f 1 | tr '\n' '@'  | sed 's/>/\n>/g' | tr -d "-" >.SEQFILE 
		#grep -v "^$" $SEQFILE | cut -d " " -f 1 | tr '\n' '@'  | sed 's/>/\n>/g'  >.SEQFILE 
		echo "###################################################"
		echo 3\-Output file\:
		read OUTFILE
		OUTFILE2=$(echo $OUTFILE | tr ' ' '_' | sed 's/^/_STAT_/'  |sed "s/^/$TS/g"  )
		OUTFILE3=$(echo $OUTFILE | tr ' ' '_' | sed 's/^/_SEQ_/'  |sed "s/^/$TS/g" )
		echo "###################################################"
		echo Stats will be saved in $OUTFILE2 \- Sequences will be saved in $OUTFILE3
		echo "###################################################"
		echo -e \#$TS >$OUTFILE2.txt
		echo -e \#Imput file\: $SEQFILE >>$OUTFILE2.txt
		echo -e \#Dom. of interest\: >>$OUTFILE2.txt
		#cat $PFDom | grep -v "#" | sed 's/^/#/g' >>$OUTFILE2.txt
		gtr '\n' '@' < Pfam-A.hmm | sed 's/\/\//\/\/\n/g'  >.mPfam
	 	grep -F -f $PFDom .mPfam | cut -d "@" -f 3,4 | tr '@' '\t' | awk '{print $4"\t"$2}' | sed 's/^/#/g' >>$OUTFILE2.txt
		echo "######################################################################################################" >>$OUTFILE2.txt
		echo -e PFam"\t"PFlen"\t"Query"\t"Qlen"\t"DomEval"\t"PFfrom"\t"PFto"\t"Qfrom"\t"Qto"\t"AliRatio >>$OUTFILE2.txt
		hmmsearch --domtblout .POTHit.txt -E 10e-3 .db1 $SEQFILE
		grep -v "#" .POTHit.txt | cut -d " " -f 1 | sort -k 1 | uniq | grep -v "^$" | sed 's/$/@/g' > .POTHit2.txt
			if [[ -s .POTHit2.txt ]]; then
			grep -F -f .POTHit2.txt .SEQFILE | tr '@' '\n' >.POTHit3.txt
			hmmscan --domtblout .temp1.4 Pfam-A.hmm .POTHit3.txt
			cat .temp1.4 | grep -v '^#' | awk '{print $2,$3,$4,$6,$13,$16,$17,$18,$19}' | sed 's/ /\t/g' | sort -k 3,3 -k 8n -k 9n | perl -e 'while(<>){chomp;@a=split;next if $a[-1]==$a[-2];push(@{$b{$a[2]}},$_);}foreach(sort keys %b){@a=@{$b{$_}};for($i=0;$i<$#a;$i++){@b=split(/\t/,$a[$i]);@c=split(/\t/,$a[$i+1]);$len1=$b[-1]-$b[-2];$len2=$c[-1]-$c[-2];$len3=$b[-1]-$c[-2];if($len3>0 and ($len3/$len1>0.5 or $len3/$len2>0.5)){if($b[4]<$c[4]){splice(@a,$i+1,1);}else{splice(@a,$i,1);}$i=$i-1;}}foreach(@a){print $_."\n";}}' | perl -e 'while(<>){chomp;@a=split(/\t/,$_);if(($a[-1]-$a[-2])>80){print $_,"\t",($a[-3]-$a[-4])/$a[1],"\n" if $a[4]<=1e-5;}else{print $_,"\t",($a[-3]-$a[-4])/$a[1],"\n" if $a[4]<=1e-3;}}' | awk '$NF>0.3' |sort -k 3 -k 8,9g | grep -F -f $PFDom > .temp1.5
		#set the evalue (by default 1e-3) and the HMM converage (by default >30%)

			cut -f 3 .temp1.5  | sort -k 1 | uniq > .POTHit4.txt
				if [[ -s .temp1.5 ]]; then
				cat .temp1.4 | grep -v '^#' | awk '{print $2,$3,$4,$6,$13,$16,$17,$18,$19}' | sed 's/ /\t/g' | sort -k 3,3 -k 8n -k 9n | perl -e 'while(<>){chomp;@a=split;next if $a[-1]==$a[-2];push(@{$b{$a[2]}},$_);}foreach(sort keys %b){@a=@{$b{$_}};for($i=0;$i<$#a;$i++){@b=split(/\t/,$a[$i]);@c=split(/\t/,$a[$i+1]);$len1=$b[-1]-$b[-2];$len2=$c[-1]-$c[-2];$len3=$b[-1]-$c[-2];if($len3>0 and ($len3/$len1>0.5 or $len3/$len2>0.5)){if($b[4]<$c[4]){splice(@a,$i+1,1);}else{splice(@a,$i,1);}$i=$i-1;}}foreach(@a){print $_."\n";}}' | perl -e 'while(<>){chomp;@a=split(/\t/,$_);if(($a[-1]-$a[-2])>80){print $_,"\t",($a[-3]-$a[-4])/$a[1],"\n" if $a[4]<=1e-5;}else{print $_,"\t",($a[-3]-$a[-4])/$a[1],"\n" if $a[4]<=1e-3;}}' | awk '$NF>0.3' |sort -k 3 -k 8,9g  >>$OUTFILE2.txt
		#set the evalue (by default 1e-3) and the HMM converage (by default >30%)
				grep -F -f .POTHit4.txt .SEQFILE | tr '@' '\n' > $OUTFILE3.txt
				else
				echo NO HIT\!
				fi
			fi
	fi
echo "######################################################################################################" >>$OUTFILE2.txt
echo \#Domain Identified\: >>$OUTFILE2.txt
cut -f 1 $OUTFILE2.txt | sort -k 1 | uniq > .allPF

grep -F -f .allPF .mPfam | cut -d "@" -f 3,4 | tr '@' '\t' | awk '{print $4"\t"$2}'  | sort -k 2| sed 's/^/#/g'>>$OUTFILE2.txt


rm -r  .POTHit*  .SEQFILE .db1 .temp* .allPF


