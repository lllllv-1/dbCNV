#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long qw/GetOptions/;
use List::Util qw/sum min max/;
use List::MoreUtils qw/uniq/;
use Statistics::Sequences::Runs 0.10;
use Statistics::Distributions;
use Cwd;
my $dir = getcwd;
my $usage=<<'USAGE';
Usage:

    perl dbCNV.pl  [Options]
   
    -i File    CNV region file (Chr Start End Type)
    -n Number  the number of disease type (2 or 5)
    -r File    the file of Rscript path
    -h Help    help information

Example
   perl dbCNV.pl -i ./cnv_file.txt -n 2 -r /usr/bin/Rscript
USAGE

my ($in,$n,$help,$Rscript);
GetOptions(
    'i=s' =>\$in,
    'n=i'=>\$n,
    'r=s'=>\$Rscript,
    'help|?' => \$help,
);

die "$usage\n" if ($help || !$in || !$n || !$Rscript);
 
my %hash;
my $gene_feature="$dir/Database/gene_feature.txt";
open GENE,"<$gene_feature"or die "Can't open $gene_feature: $!";
while(<GENE>){
    chomp;
    my @nn=split;
    $hash{$nn[0]}{"$nn[1]\_$nn[2]"}=$_;
}
close GENE;
my %gene;
my $gene_types="$dir/Database/chr_region.txt";
open TYPE,"<$gene_types" or die "Can't open $gene_types: $!";
while(<TYPE>){
    chomp;
    my @line=split;
    $gene{$line[0]}{"$line[1]\_$line[2]"}=$line[3] ;
}
close TYPE;
my %snp;
my $snp="$dir/Database/snp_feature.txt";
open SNP,"<$snp" or die "Can't open $snp: $!" ;
while(<SNP>){
    chomp;
    my @snp_line=split;
    next if($snp_line[1]=~/Start/);
    $snp{$snp_line[0]}{"$snp_line[1]\_$snp_line[2]"}=$_; 
}
close SNP;
my %element;
my $element="$dir/Database/hg19-cCREs.txt";
open CRE,"<$element" or die "Can't open $element: $!";
while(<CRE>){
    chomp;
    my @CRE_line=split;
    $element{$CRE_line[0]}{"$CRE_line[1]\_$CRE_line[2]"}=$CRE_line[-1];
}
close CRE;

my @Clingen_value_HI=(0,1,2,3,30,40,"NA");
my @Clingen_value_TS=(0,1,2,3,40,"NA");
open IN,"<$in" or die "Can't open $in: $!";
open LOSS,">>$dir/loss\_$n\_test.txt";
open GAIN,">>$dir/gain\_$n\_test.txt";
my $gain_file="N",my $loss_file="N",
my $header="Chr\tStart\tEnd\tLength\tGene_sum\tHI_0_sum\tHI_1_sum\tHI_2_sum\tHI_3_sum\tHI_30_sum\tHI_40_sum\tHI_NA_sum\tTS_0_sum\tTS_1_sum\tTS_2_sum\tTS_3_sum\tTS_40_sum\tTS_NA_sum\tGHIS_mean\tEpiscore_max\tLOEUF_min\tpLI_max\tMorbid_gene_sum\tOMIM_gene_sum\tZscore_gene_max\tthree_prime_UTR_num\tfive_prime_UTR_num\tlincRNA_gene_num\tmiRNA_gene_num\tpseudogene_num\trRNA_gene_num\tsnoRNA_gene_num\tsnRNA_gene_num\tSIFT_score_min\tSIFT_pred_D\tSIFT_pred_T\tPolyphen2_HDIV_score_max\tPolyphen2_HDIV_pred_D\tPolyphen2_HDIV_pred_pred_P\tPolyphen2_HDIV_pred_B\tPolyphen2_HVAR_score_max\tPolyphen2_HVAR_pred_D\tPolyphen2_HVAR_pred_P\tPolyphen2_HVAR_pred_B\tLRT_score_max\tLRT_pred_D\tLRT_pred_N\tLRT_pred_U\tMutationTaster_score_max\tMutationTaster_pred_A\tMutationTaster_pred_D\tMutationTaster_pred_N\tMutationTaster_pred_P\tMutationAssessor_score_max\tMutationAssessor_pred_H\tMutationAssessor_pred_L\tMutationAssessor_pred_M\tMutationAssessor_pred_N\tFATHMM_score_min\tFATHMM_pred_D\tFATHMM_pred_T\tRadialSVM_score_max\tRadialSVM_pred_D\tRadialSVM_pred_T\tLR_score_max\tLR_pred_D\tLR_pred_T\tVEST3_score_max\tCADD_raw_max\tCADD_phred_max\tGERP_RS_max\tphyloP46way_placental_max\tphyloP100way_vertebrate_max\tSiPhy_29way_logOdds_max\tCDTS_min\trevel_max\tCTCF-bound_num\tCTCF-only_num\tdELS_num\tDNase-H3K4me3_num\tpELS_num\tPLS_num\tType";
   ##Define the median value
my %median;
$median{"19"}=0.52;$median{"20"}=0.61;$median{"21"}=0.36;$median{"22"}=0.9;
$median{"25"}=1.79;$median{"34"}=0;$median{"37"}=1;$median{"41"}=1;
$median{"45"}=0.982;$median{"49"}=1;$median{"54"}=4.31;$median{"59"}=-4.68;
$median{"62"}=1.098;$median{"65"}=0.968;$median{"68"}=0.996;$median{"69"}=9.052;
$median{"70"}=41;$median{"71"}=5.8;$median{"72"}=2.773;$median{"73"}=9.434;
$median{"74"}=20.256;$median{"75"}=-8.82154;$median{"76"}=0.967;
print LOSS "$header\n";
print GAIN "$header\n";
while(<IN>){
   chomp;
   my @string=();
   my @region=split/\s+/,$_;
   next if ($region[-1]=~/Type/);
   my $len=$region[2]-$region[1]+1;
   push @string,"$region[0]\t$region[1]\t$region[2]\t$len";
   my (%HI,%TS,@GHIS,@Episcores,@LOEUF,@pLI,@Zscore_gene);
   my ($gene_num,$Morbid_gene,$OMIM_gene)=(0,0,0);
   ##Gene feature calculation
   if (exists$hash{$region[0]}){
      for my $st_ed(keys%{$hash{$region[0]}}){ 
         my @location=split /\_/,$st_ed; 
         if(($location[0]>=$region[1])and($location[1]<=$region[2])){
           push @GHIS,"NA";push @Episcores,"NA";push @LOEUF,"NA";push @pLI,"NA";push @Zscore_gene,"NA";
           $gene_num++;
           my $gene_value=$hash{$region[0]}{$st_ed};
           my @value=split /\t/,$gene_value;
           $HI{$value[5]}+=1;$TS{$value[6]}+=1;
             if($value[7] ne "NA"){
             pop @GHIS;
             push @GHIS,$value[7];
             }else{pop @GHIS;}
             if($value[8] ne "NA"){
             pop @Episcores;
             push @Episcores,$value[8];
             } if($value[8] eq "NA"){ pop @Episcores;} 
             if($value[9] ne "NA"){
             pop @LOEUF;
             push @LOEUF,$value[9];
             }else {pop @LOEUF;}
             if($value[10] ne "NA"){
             pop @pLI;
             push @pLI,$value[10];
             }else {pop @pLI;}
             $Morbid_gene+=$value[11];
             $OMIM_gene+=$value[12];
             if($value[13] !~/NA/){
             pop @Zscore_gene; 
             push @Zscore_gene,$value[13];
             }else{pop @Zscore_gene;}
         }
      }
  }
    push @string, "$gene_num";
    for my $i(@Clingen_value_HI){
        if(exists$HI{$i}){push @string, "$HI{$i}";}
        if(!exists$HI{$i}){push @string, "0";}
  }
    for my $j(@Clingen_value_TS){
       if(exists$TS{$j}){push @string, "$TS{$j}";} 
       if(!exists$TS{$j}){push @string, "0";}
  }   
   if (!defined $GHIS[0]){push @string, "NA";}
    else{ 
       my $GHIS_mean=&mean(@GHIS);
       $GHIS_mean=sprintf "%.2f",$GHIS_mean;
       push @string, "$GHIS_mean";
  }
   if (!defined $Episcores[0]){push @string, "NA";}
    else{ my $Episcores_max=max@Episcores;
        $Episcores_max=sprintf "%.2f",$Episcores_max;
        push @string, "$Episcores_max";
  }  
    if (!defined $LOEUF[0]){push @string, "NA";}  
    else{ my $LOEUF_min=min@LOEUF;
       push @string, "$LOEUF_min";
 }
    if (!defined $pLI[0]){push @string, "NA";}
     else{my $pLI_max=max@pLI;
       push @string, "$pLI_max";
 }
       push @string, "$Morbid_gene";
       push @string, "$OMIM_gene";
    if (!defined $Zscore_gene[0]){push @string, "NA";}
       else{my $Zscore_gene_max=max@Zscore_gene;
       push @string, "$Zscore_gene_max";
 }
     
     my %gene_type_num;
    if (exists$gene{$region[0]}){
         for my $start_end(keys%{$gene{$region[0]}}){
         my @loc=split /\_/,$start_end;
           if(($loc[0]>=$region[1])and($loc[1]<=$region[2])){
           my $gene_types=$gene{$region[0]}{$start_end};
           $gene_type_num{$gene_types}+=1;
       }
     }
   }  
     #three_prime_UTR_num\tfive_prime_UTR_num\tlincRNA_gene\tmiRNA_gene\tpseudogene\trRNA_gene\tsnoRNA_gene\tsnRNA_gene\t
    if (exists$gene_type_num{"three_prime_UTR"}){push @string, "$gene_type_num{'three_prime_UTR'}";}
    if (!exists$gene_type_num{"three_prime_UTR"}){push @string, "0";}
    if (exists$gene_type_num{"five_prime_UTR"}){push @string, "$gene_type_num{'five_prime_UTR'}";}
    if (!exists$gene_type_num{"five_prime_UTR"}){push @string, "0";}
    if (exists$gene_type_num{"lincRNA_gene"}){push @string, "$gene_type_num{'lincRNA_gene'}";}
    if (!exists$gene_type_num{"lincRNA_gene"}){push @string, "0";}
    if (exists$gene_type_num{"miRNA_gene"}){push @string, "$gene_type_num{'miRNA_gene'}";}
    if (!exists$gene_type_num{"miRNA_gene"}){push @string, "0";}
    if (exists$gene_type_num{"pseudogene"}){push @string, "$gene_type_num{'pseudogene'}";}
    if (!exists$gene_type_num{"pseudogene"}){push @string, "0";}
    if (exists$gene_type_num{"rRNA_gene"}){push @string, "$gene_type_num{'rRNA_gene'}";}
    if (!exists$gene_type_num{"rRNA_gene"}){push @string, "0";}
    if (exists$gene_type_num{"snoRNA_gene"}){push @string, "$gene_type_num{'snoRNA_gene'}";}
    if (!exists$gene_type_num{"snoRNA_gene"}){push @string, "0";}
    if (exists$gene_type_num{"snRNA_gene"}){push @string, "$gene_type_num{'snRNA_gene'}";}
    if (!exists$gene_type_num{"snRNA_gene"}){push @string,  "0";}
    ## SNP  feature calculation 
      my ($SIFT_pred_D,$SIFT_pred_T,$Polyphen2_HDIV_pred_D,$Polyphen2_HDIV_pred_P,$Polyphen2_HDIV_pred_B,$Polyphen2_HVAR_pred_D,$Polyphen2_HVAR_pred_P,$Polyphen2_HVAR_pred_B,$LRT_pred_D,$LRT_pred_N,$LRT_pred_U,$MutationTaster_pred_A,$MutationTaster_pred_D,$MutationTaster_pred_N,$MutationTaster_pred_P,$MutationAssessor_pred_H,$MutationAssessor_pred_L,$MutationAssessor_pred_M,$MutationAssessor_pred_N,$FATHMM_pred_D,$FATHMM_pred_T,$RadialSVM_pred_D,$RadialSVM_pred_T,$LR_pred_D,$LR_pred_T)=(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
      my (@SIFT_score,@Polyphen2_HDIV_score,@Polyphen2_HVAR_score,@LRT_score,@MutationTaster_score,@MutationAssessor_score,@FATHMM_score,@RadialSVM_score,@LR_score,@VEST3_score,@CADD_raw,@CADD_phred,@GERP_RS,@phyloP46way_placental,@phyloP100way_vertebrate,@SiPhy_29way_logOdds,@CDTS,@revel);
    if (exists$snp{$region[0]}){ 
      for my $st_ed(keys%{$snp{$region[0]}}){
        my @location2=split /\_/,$st_ed;
         if(($location2[0]>=$region[1])and($location2[1]<=$region[2])){
         my $snp_value=$snp{$region[0]}{$st_ed};
         my @value=split /\s+/,$snp_value;
         if($value[4] ne "NA"){push @SIFT_score,$value[4];}
         $SIFT_pred_D+=$value[5];$SIFT_pred_T+=$value[6];
         if($value[7] ne "NA"){ push @Polyphen2_HDIV_score,$value[7];}
         $Polyphen2_HDIV_pred_D+=$value[8];$Polyphen2_HDIV_pred_P+=$value[9];$Polyphen2_HDIV_pred_B+=$value[10];
         if($value[11] ne "NA"){ push @Polyphen2_HVAR_score,$value[11];}
         $Polyphen2_HVAR_pred_D+=$value[12];$Polyphen2_HVAR_pred_P+=$value[13];$Polyphen2_HVAR_pred_B+=$value[14];
         if($value[15] ne "NA"){  push @LRT_score,$value[15];}
         $LRT_pred_D+=$value[16];$LRT_pred_N+=$value[17];$LRT_pred_U+=$value[18];
         if($value[19] ne "NA"){ push @MutationTaster_score,$value[19];}
         $MutationTaster_pred_A+=$value[20];$MutationTaster_pred_D+=$value[21];$MutationTaster_pred_N+=$value[22];$MutationTaster_pred_P+=$value[23];
         if($value[24] ne "NA"){  push @MutationAssessor_score,$value[24];}
          $MutationAssessor_pred_H+=$value[25];$MutationAssessor_pred_L+=$value[26];$MutationAssessor_pred_M+=$value[27];$MutationAssessor_pred_N+=$value[28];
         if($value[29] ne "NA"){  push @FATHMM_score,$value[29];}
         $FATHMM_pred_D+=$value[30]; $FATHMM_pred_T+=$value[31];
         if($value[32] ne "NA"){  push @RadialSVM_score,$value[32];}
         $RadialSVM_pred_D+=$value[33];$RadialSVM_pred_T+=$value[34];
         if($value[35] ne "NA"){ push @LR_score,$value[35];}
         $LR_pred_D+=$value[36];$LR_pred_T+=$value[37];
         if($value[38] ne "NA"){ push @VEST3_score,$value[38];}
         if($value[39] ne "NA"){ push @CADD_raw,$value[39];}
         if($value[40] ne "NA"){ push @CADD_phred,$value[40];}
         if($value[41] ne "NA"){ push @GERP_RS,$value[41];}
         if($value[42] ne "NA"){ push @phyloP46way_placental,$value[42];}
         if($value[43] ne "NA"){ push @phyloP100way_vertebrate,$value[43];}
         if($value[44] ne "NA"){ push @SiPhy_29way_logOdds,$value[44];}
         if($value[45] ne "NA"){ push @CDTS,$value[45];}
         if($value[46] ne "NA"){ push @revel,$value[46];}
      }
    }
  }  
         if (!defined $SIFT_score[0]){push @string, "NA";}
         else{ my $SIFT_score_min=min@SIFT_score;
         push @string, "$SIFT_score_min";}
         push @string, "$SIFT_pred_D\t$SIFT_pred_T";
         if (!defined $Polyphen2_HDIV_score[0]){push @string, "NA";}
         else{ my $Polyphen2_HDIV_score_max=max@Polyphen2_HDIV_score;
         push @string, "$Polyphen2_HDIV_score_max";}
         push @string, "$Polyphen2_HDIV_pred_D\t$Polyphen2_HDIV_pred_P\t$Polyphen2_HDIV_pred_B";
         if (!defined $Polyphen2_HVAR_score[0]){push @string, "NA";}
         else{ my $Polyphen2_HVAR_score_max=max@Polyphen2_HVAR_score;
         push @string, "$Polyphen2_HVAR_score_max";} 
         push @string, "$Polyphen2_HVAR_pred_D\t$Polyphen2_HVAR_pred_P\t$Polyphen2_HVAR_pred_B";
         if (!defined $LRT_score[0]){push @string, "NA";}
          else {my $LRT_score_max=max@LRT_score;
         push @string, "$LRT_score_max";}
         push @string, "$LRT_pred_D\t$LRT_pred_N\t$LRT_pred_U";
         if (!defined $MutationTaster_score[0]){push @string, "NA";}
           else {my $MutationTaster_score_max=max@MutationTaster_score;
         push @string, "$MutationTaster_score_max";}
         push @string, "$MutationTaster_pred_A\t$MutationTaster_pred_D\t$MutationTaster_pred_N\t$MutationTaster_pred_P";
         if (!defined $MutationAssessor_score[0]){push @string, "NA";}         
            else {my $MutationAssessor_score_max=max@MutationAssessor_score;
         push @string, "$MutationAssessor_score_max";}
         push @string, "$MutationAssessor_pred_H\t$MutationAssessor_pred_L\t$MutationAssessor_pred_M\t$MutationAssessor_pred_N";
         if (!defined $FATHMM_score[0]){push @string, "NA";}                    
             else {my $FATHMM_score_min=min@FATHMM_score;
         push @string, "$FATHMM_score_min";}
         push @string, "$FATHMM_pred_D\t$FATHMM_pred_T";
         if (!defined $RadialSVM_score[0]){push @string, "NA";}
           else {my  $RadialSVM_score_max=max@RadialSVM_score;
         push @string, "$RadialSVM_score_max";}
         push @string, "$RadialSVM_pred_D\t$RadialSVM_pred_T";
         if (!defined $LR_score[0]){push @string, "NA";}
            else {my $LR_score_max=max@LR_score;
         push @string, "$LR_score_max";}
         push @string, "$LR_pred_D\t$LR_pred_T";

         if (!defined $VEST3_score[0]){push @string,"NA";}
                else {my  $VEST3_score_max=max@VEST3_score;
                    push @string, "$VEST3_score_max";}
         if (!defined $CADD_raw[0]){push @string, "NA";}
                  else {my $CADD_raw_max=max(@CADD_raw);
                    push @string, "$CADD_raw_max";}
         if (!defined $CADD_phred[0]){push @string, "NA";}
                  else {my $CADD_phred_max=max@CADD_phred;
                     push @string, "$CADD_phred_max";}
         if (!defined $GERP_RS[0]){push @string, "NA";}
                   else {my $GERP_RS_max=max@GERP_RS;
                     push @string, "$GERP_RS_max";}
         if (!defined $phyloP46way_placental[0]){push @string, "NA";}
                 else {my $phyloP46way_placental_max=max@phyloP46way_placental;
                push @string, "$phyloP46way_placental_max";}
         if (!defined $phyloP100way_vertebrate[0]){push @string, "NA";}
               else {my $phyloP100way_vertebrate_max=max@phyloP100way_vertebrate;
               push @string, "$phyloP100way_vertebrate_max";}
         if (!defined $SiPhy_29way_logOdds[0]){push @string, "NA";}
               else {my $SiPhy_29way_logOdds_max=max@SiPhy_29way_logOdds;
                push @string, "$SiPhy_29way_logOdds_max";}
         if (!defined $CDTS[0]){push @string, "NA";}
                 else {my $CDTS_min=min@CDTS; 
                push @string, "$CDTS_min";}
         if (!defined $revel[0]){push @string, "NA";}
           else {my $revel_max=max@revel;
           push @string, "$revel_max";}
      ##Hg19 cCREs element number calculation
          my %element_num;
         if (exists$element{$region[0]}){
            for my $s_e(keys%{$element{$region[0]}}){
              my @location3=split /\_/,$s_e;
               if(($location3[0]>=$region[1])and($location3[1]<=$region[2])){
                my $element_value=$element{$region[0]}{$s_e};
                my @ele=split /\,/,$element_value;
            for my $i(@ele){
                $element_num{$i}++;
           }
          }
        }
      }
      for my $j('CTCF-bound','CTCF-only','dELS','DNase-H3K4me3','pELS','PLS'){
       if (exists$element_num{$j}){push @string, "$element_num{$j}";}else {push @string, "0";} 
     } 
      my $loss_n=0;
      my $gain_n=0;
      if($region[3] eq "loss") {
          $loss_file="Y";
         for my $loss(@string){
           my @loss_line=split /\t/,$loss;
           for my $loss_i(@loss_line){   
           $loss_n++;
           if ($loss_i eq "NA"){print LOSS "$median{$loss_n}\t";}else{ print LOSS "$loss_i\t";}
         }    
       }
           print LOSS "$region[3]\n";  
      }
       if($region[3] eq "gain") {
          $gain_file="Y";
        for my $gain(@string){ 
           my @gain_line=split /\t/,$gain;
           for my $gain_i(@gain_line){
           $gain_n++;
           if ($gain_i eq "NA"){print GAIN "$median{$gain_n}\t";}else{ print GAIN "$gain_i\t";}
          }
        }
           print GAIN "$region[3]\n";
      } 
   }
close IN;
close LOSS;
close GAIN;
  ###Output the R script
my $trainset_2_gain_file="$dir/Database/trainset_2_gain.txt";
my $trainset_2_loss_file="$dir/Database/trainset_2_loss.txt";
my $trainset_5_gain_file="$dir/Database/trainset_5_gain.txt";
my $trainset_5_loss_file="$dir/Database/trainset_5_loss.txt";
my $head="library(xgboost)\nlibrary(caret)\nlibrary(ggplot2)\nlibrary(lattice)\nlibrary(Matrix)\nset.seed(1234)\n";
my $train_code="trainset1 <- data.matrix(trainset[,c(4:82)])\ntrainset2 <- Matrix(trainset1,sparse = T)\ntrain_y <- factor(trainset[,84])\nf <- train_y\ntrain_y <- as.numeric(f)-1\ntraindata <- list(data=trainset2,label=train_y)\ndtrain <- xgb.DMatrix(data=traindata\$data,label=traindata\$label)\n";
my $test_code="testset1 <- data.matrix(testset[,c(4:82)])\ntestset2 <- Matrix(testset1,sparse = T)\ndtest <- xgb.DMatrix(data=testset2)\npre <- predict(fit,newdata = dtest)\n";

if(($n==5)and($gain_file eq "Y")){
   open  G5,">$dir/predication_5_gain.R";
   print G5 "$head";
   print G5 "trainset <- read.table(\"$trainset_5_gain_file\",sep=\"\\t\",header=T)\n";
   print G5 "$train_code";
   print G5 "params <- list(eta = 0.6801917,gamma = 0,max_depth = 10, min_child_weight = 1,subsample = 1, nfold = 10, objective = \"multi:softmax\",num_class=5)\n";  
   print G5 "fit <- xgboost(params = params,data = trainset1, label = train_y,nrounds = 34, eval_metric = \"mlogloss\")\n";
   print G5 "testset <- read.table(\"./gain_5_test.txt\",sep=\"\\t\",header=T)\n";
   print G5 "$test_code";
   print G5 "write.table(pre,\"$dir/gain_5_predication.txt\",sep=\"\\t\")\n";
   `$Rscript $dir/predication_5_gain.R`;
   &trans("$dir/gain_5_predication.txt","$dir/gain_5_test.txt");
  `rm $dir/gain_5_predication.txt`;
  `rm $dir/predication_5_gain.R`; 
}
if (-e "$dir/gain_5_test.txt"){`rm $dir/gain_5_test.txt`;}
close G5;

if(($n==5)and($loss_file eq "Y")){
   open  L5,">$dir/predication_5_loss.R";
   print L5 "$head";
   print L5 "trainset <- read.table(\"$trainset_5_loss_file\",sep=\"\\t\",header=T)\n";
   print L5 "$train_code";
   print L5 "params <- list(eta = 0.4662765,gamma = 0,max_depth = 10, min_child_weight = 25,subsample = 0.25, nfold = 10, objective = \"multi:softmax\",num_class=5)\n";
   print L5 "fit <- xgboost(params = params,data = trainset1, label = train_y,nrounds = 94, eval_metric = \"mlogloss\")\n";
   print L5 "testset <- read.table(\"./loss_5_test.txt\",sep=\"\\t\",header=T)\n";
   print L5 "$test_code";
   print L5 "write.table(pre,\"$dir/loss_5_predication.txt\",sep=\"\\t\")\n";
  `$Rscript predication_5_loss.R`;
  &trans("$dir/loss_5_predication.txt","$dir/loss_5_test.txt");
  `rm $dir/loss_5_predication.txt`;
  `rm $dir/predication_5_loss.R`;
}
if (-e "$dir/loss_5_test.txt"){`rm $dir/loss_5_test.txt`;}
close L5;

if(($n==2)and($gain_file eq "Y")){
   open G2,">$dir/predication_2_gain.R";
   print G2 "$head";
   print G2 "trainset <- read.table(\"$trainset_2_gain_file\",sep=\"\\t\",header=T)\n";
   print G2 "$train_code";
   print G2 "params <- list(eta = 1,gamma = 0,max_depth = 10, min_child_weight = 1,subsample = 1, nfold = 10, objective = \"binary:logistic\")\n";
   print G2 "fit <- xgboost(params = params,data = trainset1, label = train_y,nrounds = 72, eval_metric = \"auc\")\n";
   print G2 "testset <- read.table(\"./gain_2_test.txt\",sep=\"\\t\",header=T)\n";
   print G2 "$test_code";
   print G2 "pre <- as.numeric(pre > 0.50)\n";
   print G2 "write.table(pre,\"$dir/gain_2_predication.txt\",sep=\"\\t\")\n";
   `$Rscript predication_2_gain.R`;
   &trans("$dir/gain_2_predication.txt","$dir/gain_2_test.txt");
   `rm $dir/gain_2_predication.txt`;
   `rm $dir/predication_2_gain.R`;
}
if (-e "$dir/gain_2_test.txt"){`rm $dir/gain_2_test.txt`;}
close G2;

if(($n==2)and($loss_file eq "Y")){
   open L2,">$dir/predication_2_loss.R";
   print L2 "$head";
   print L2  "trainset <- read.table(\"$trainset_2_loss_file\",sep=\"\\t\",header=T)\n";
   print L2 "$train_code";
   print L2 "params <- list(eta = 0.5777781,gamma = 0,max_depth = 10, min_child_weight = 4.920737,subsample = 1, nfold = 10, objective = \"binary:logistic\")\n";
   print L2 "fit <- xgboost(params = params,data = trainset1, label = train_y,nrounds = 98, eval_metric = \"auc\")\n";
   print L2 "testset <- read.table(\"./loss_2_test.txt\",sep=\"\\t\",header=T)\n";
   print L2 "$test_code";
   print L2 "pre <- as.numeric(pre > 0.50)\n";
   print L2 "write.table(pre,\"$dir/loss_2_predication.txt\",sep=\"\\t\")\n";
   `$Rscript predication_2_loss.R`;
   &trans("$dir/loss_2_predication.txt","$dir/loss_2_test.txt"); 
  `rm $dir/loss_2_predication.txt`;
  `rm $dir/predication_2_loss.R`;
}
if (-e "$dir/loss_2_test.txt"){`rm $dir/loss_2_test.txt`;}
close L2;
   

sub mean {
   my (@data)=@_;
   my $mean_sum;
   foreach my $d(@data){
      $mean_sum+=$d;
   }
      $mean_sum/=($#data+1);
    return $mean_sum;
}

sub trans {
   my ($file1, $file2) = @_;
   open DN1,"<$file1" or die $!;
   open DN2,"<$file2" or die $!;
my %result;
my %disease;
if($n==5){
     $disease{"0"}="benign";$disease{"1"}="likely benign";$disease{"2"}="likely pathogenic";
     $disease{"3"}="pathogenic";$disease{"4"}="uncertain significance";
}else{
     $disease{"0"}="benign";$disease{"1"}="pathogenic";
     }
while(<DN1>){
    chomp;
    my @line=split /\s+/,$_;
    next if ($#line<1);
    my $num=$line[0];
    $num=~s/[""]//g;
    $result{$num}=$line[-1];
}
close DN1;
my @name=split /\_/,(split /\//,$file1)[-1];

open OU2,">$dir/$name[0]\_$n\_predication_result.txt" or die $!;
  my $cnv_num=-1;
  while(<DN2>){
  chomp;
  $cnv_num++;
  my @outline=split /\s+/,$_;
  if($outline[-1]=~/Type/) {print OU2 "$outline[0]\t$outline[1]\t$outline[2]\t$outline[-1]\tDisease\n";}
    else {
    my $state=$result{$cnv_num};
    print OU2 "$outline[0]\t$outline[1]\t$outline[2]\t$outline[-1]\t$disease{$state}\n";}
}
close DN2;  
close OU2;
}
