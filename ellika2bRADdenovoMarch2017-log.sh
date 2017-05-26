# Trimming - Cutting out restriction fragments	from the reads
../scripts/2bRAD_trim_launch.pl fastq sampleID=1 > ../scripts/trims.sh
# open file and add # !/bin/bash to show it's a bash file before running!
sh ../scripts/trims.sh

#Filtering - trims out low quality read ends
../scripts/2bRAD_filt_launch.sh
# open file and add # !/bin/bash to show it's a bash file before running!
sh ../scripts/filt.sh

# Uniquing - Removs duplicates, only keeping unique reads
../scripts/2bRAD_uni_launch.sh
# open file and add # !/bin/bash to show it's a bash file before running!
sh ../scripts/unii.sh

# Merges all the .uni files to a file called Merged.uniq, where count information is stored about how many reads each unique sequence represents for each individual (used for genotyping later)
#also creates file called “mergedUniqTags.fasta”, which only contains the sequences (used for clustering)
~/scripts/mergeUniq.pl uni >Merged.uniq

#Clustering with cd-hit-est:
#identity threshold -c 91 %, which translates to about 3 nucleotide substitutions allowed.
~/scripts/cd-hit-est -i mergedUniqTags.fasta -o cdh_alltags.fas -aL 1 -aS 1 -g 1 -c 0.91 -M 0 -T 0

#Merging the clustering output into alltags.ul:
#Merges clusters with the counts data for each individual (stored in Merged.uniq), in order to be able to genotype.
#alltags.ul contains information about loci, tags (alleles), the sequences, and counts for all individuals 
#revcom=reverse complement of sequence.
~/scripts/uniq2loci.pl Merged.uniq cdh_alltags.fas.clstr cdh_alltags.fas minDP=20 >alltags.ul

#Calling SNPs with haplocall_denovo, using filter settings recommended by Mischa Matz (cranked up as high_mem ran out of memory with default settings):
#Mindp= Minimum depth for calling a homozygote (5), 
#ind= Number of individuals should have the allele to call it (2)
# hetero=fraction heterozygotes allowed (0.8)
#This will output information about the number of loci and alleles 
#on the screen, and will save genotype data both for SNPs and 
#haplotypes in two separate files in the .vcf file format (variant call format).
#-----better to use qsub than qlogin!
~/scripts/haplocall_denovo.pl alltags.ul mindp=5 ind=2 hetero=0.8 aobs=20 strbias=20 found=0.25

# OUTPUT::
#	5447085	total tags
#	1509430	pass depth cutoff 20
#	1256149	pass strand bias cutoff 5
#	reading start: Mon Apr 24 22:33:17 2017	
#		
#	allele filters start: Mon Apr 24 22:57:57 2017	
#		
#	locus filters start: Mon Apr 24 23:36:17 2017	
#	VCF writing start: Mon Apr 24 23:40:15 2017	
#	Haplo-VCF writing start: Mon Apr 24 23:44:58 2017	
#	all done: Mon Apr 24 23:48:42 2017	
#	Input parameters (parameter names):	
#	clipping 0 bases from read ends (clip)	
#	must see an allele at least 20 times (aobs)	
#	strand bias cutoff (strbias): 20	
#	allele bias cutoff (abias): 10	
#	keep loci genotyped in at least 0.25 fraction of all samples (found)	
#	must see an allele in at least 2 individual(s) (ind)	
#	keep monomorphic tags? (mono): 0 	
#	maximum acceptable fracton of heterozygotes at a locus (hetero): 0.8	
#	minimum depth to call a homozygote (mindp): 5	
#		
#		file alltags.ul
#		
#	-----------------	
#		
#	Allele filters:	
#		
#	1256111	raw alleles
#	1256111	with 20 or more reads
#	1251336	in 2 or more samples
#	1123858	pass strand bias cutoff 20
#	323308	pass allele bias cutoff 10
#		
#	--------------------	
#		
#	Locus filters:	
#		
#	245594	total
#	245594	remain after applying allele filters
#	242280	have less than 0.8 fraction of heterozygotes
#	235430	with not more than 2 alleles per individual
#	79412	genotyped in 25% of samples
#	11929	polymorphic

# make a tab-delimited table of clone (technical replicate) sample pairs called clonepairs.tab
nano clonepairs.tab
#

# using technical replicate in order to optimize dataset using a machine-learning algorithm.
# Each of the following filtering steps creates a vcf output file that is a modified and renamed version of the input file.

# extracting "true snps" subset (reproducible across replicates)
~/scripts/replicatesMatch.pl vcf=alltags.ul_Variants_count20_ab10_sb20_clip0.vcf replicates=clonepairs.tab polyonly=1 > vqsr.denovo.vcf


# 12220 total SNPs 1180 pass hets and match filters 589 show non-reference alleles 197 have alterantive alleles in at least 2 replicate pair(s) 197 have matching heterozygotes in at least 0 replicate pair(s) 197 polymorphic 
# 197 written
 
#non-parametric recalibration (strand bias not used, -nosb) the goal is to generate a composite filter and determine the setting that gives maximum "gain"
~/scripts/recalibrateSNPs.pl vcf=alltags.ul_Variants_count20_ab10_sb20_clip0.vcf true=vqsr.denovo.vcf -nosb -notp >denovo.recal.vcf

# OUTPUT::
# 22356494        total tags
# 7051292 pass depth cutoff 20
# 6799207 pass strand bias cutoff 5
# reading start: Thu Aug 18 16:47:18 2016
# 
# allele filters start: Thu Aug 18 19:00:52 2016
# 
# locus filters start: Thu Aug 18 22:59:57 2016
# VCF writing start: Thu Aug 18 23:44:10 2016
# Haplo-VCF writing start: Thu Aug 18 23:58:15 2016
# all done: Fri Aug 19 00:07:39 2016
# Input parameters (parameter names):
# clipping 0 bases from read ends (clip)
# must see an allele at least 20 times (aobs)
# strand bias cutoff (strbias): 20
# allele bias cutoff (abias): 10
# keep loci genotyped in at least 0.5 fraction of all samples (found)
# must see an allele in at least 10 individual(s) (ind)
# keep monomorphic tags? (mono): 0
# maximum acceptable fracton of heterozygotes at a locus (hetero): 0.8
# minimum depth to call a homozygote (mindp): 5
# 
#         file alltags.ul
# 
# -----------------
# 
# Allele filters:
# 
# 6796332 raw alleles
# 6796332 with 20 or more reads
# 5470290 in 10 or more samples
# 5288794 pass strand bias cutoff 20
# 4724303 pass allele bias cutoff 10
# 
# --------------------
# 
# Locus filters:
# 
# 3422577 total
# 3422577 remain after applying allele filters
# 3405097 have less than 0.8 fraction of heterozygotes
# 3235182         with not more than 2 alleles per individual
# 37304   genotyped in 50% of samples

# make a tab-delimited table of clone (replicate) sample pairs
nano clonepairs.tab
#edit accordingly 
#	SMID17a	SMID17b
#	SMID18a	SMID18b
#	SMID19a	SMID19b
#	SMID20a	SMID20b
#	SMID37a	SMID37b
#	SMID38a	SMID38b
#	SMID39a	SMID39b
#	SMID40a	SMID40b
#	SM2014RI16a	SM2014RI16b
#	SM2014RI17a	SM2014RI17b
#	SM2014RI18a	SM2014RI18b
#	SM2014RI19a	SM2014RI19b
#	SM2014SF17a	SM2014SF17b
#	SM2014SF18a	SM2014SF18b
#	SM2014SF19a	SM2014SF19b
#	SM2014SF20a	SM2014SF20b
#	SMAF17a	SMAF17b
#	SMAF18a	SMAF18b
#	SMAF19a	SMAF19b
#	SMAF20a	SMAF20b
#	SMAF37a	SMAF37b
#	SMAF38a	SMAF38b
#	SMAF39a	SMAF39b
#	SMAF40a	SMAF40b
#	SMFB17a	SMFB17b
#	SMFB18a	SMFB18b
#	SMFB19a	SMFB19b
#	SMFB20a	SMFB20b
#	SMFB37a	SMFB37b
#	SMFB38a	SMFB38b
#	SMFB39a	SMFB39b
#	SMFB40a	SMFB40b
#	SMKB17a	SMKB17b
#	SMKB18a	SMKB18b
#	SMKB19a	SMKB19b
#	SMKB20a	SMKB20b
#	SMKB37a	SMKB37b
#	SMKB38a	SMKB38b
#	SMKB39a	SMKB39b
#	SMKB40a	SMKB40b
#	SMLB17a	SMLB17b
#	SMLB18a	SMLB18b
#	SMLB19a	SMLB19b
#	SMLB20a	SMLB20b
#	SMLB37a	SMLB37b
#	SMLB38a	SMLB38b
#	SMLB39a	SMLB39b
#	SMLB40a	SMLB40b
#	SMLYF17a	SMLYF17b
#	SMLYF18a	SMLYF18b
#	SMLYF19a	SMLYF19b
#	SMLYF20a	SMLYF20b
#	SMLYF37a	SMLYF37b
#	SMLYF38a	SMLYF38b
#	SMLYF39a	SMLYF39b
#	SMLYF40a	SMLYF40b
#	SMOS17a	SMOS17b
#	SMOS18a	SMOS18b
#	SMOS19a	SMOS19b
#	SMOS20a	SMOS20b
#	SMOS37a	SMOS37b
#	SMOS38a	SMOS38b
#	SMOS39a	SMOS39b
#	SMOS40a	SMOS40b
#	SMRI17a	SMRI17b
#	SMRI18a	SMRI18b
#	SMRI19a	SMRI19b
#	SMRI20a	SMRI20b
#	SMRI37a	SMRI37b
#	SMRI38a	SMRI38b
#	SMRI39a	SMRI39b
#	SMRI40a	SMRI40b
#	SMSF17a	SMSF17b
#	SMSF18a	SMSF18b
#	SMSF19a	SMSF19b
#	SMSF20a	SMSF20b
#	SMSF37a	SMSF37b
#	SMSF38a	SMSF38b
#	SMSF39a	SMSF39b
#	SMSF40a	SMSF40b
#	SMST17a	SMST17b
#	SMST18a	SMST18b
#	SMST19a	SMST19b
#	SMST20a	SMST20b
#	SMST37a	SMST37b
#	SMST38a	SMST38b
#	SMST39a	SMST39b
#	SMST40a	SMST40b
#	SMTF17a	SMTF17b
#	SMTF18a	SMTF18b
#	SMTF19a	SMTF19b
#	SMTF20a	SMTF20b
#	SMTF37a	SMTF37b
#	SMTF38a	SMTF38b
#	SMTF39a	SMTF39b
#	SMTF40a	SMTF40b


# extracting "true snps" subset (reproducible across replicates)
~/scripts/replicatesMatch.pl vcf=alltags.ul_Variants_count20_ab10_sb20_clip0.vcf replicates=clonepairs.tab polyonly=1 > vqsr.denovo.vcf
#OUTPUT::
# 13254 total SNPs
# 1433 pass hets and match filters
# 806 show non-reference alleles
# 308 have alterantive alleles in at least 2 replicate pair(s)
# 308 have matching heterozygotes in at least 0 replicate pair(s)
# 308 polymorphic
# 308 written

# non-parametric recalibration,  creates a composite filter and determine the setting that gives maximum "gain"	
~/scripts/recalibrateSNPs.pl vcf=alltags.ul_Variants_count20_ab10_sb20_clip0.vcf true=vqsr.denovo.vcf -nosb -notp >denovo.recal.vcf
#OUTPUT::
#	AB quantiles:		
#	1	14	
#	5	21	
#	10	24	
#	15	29	
#	20	31	
#	30	34	
#	40	38	
#	50	41	
#	60	44	
#	70	48	
#	80	51	
#	90	57	
#			
#	SB quantiles:		
#	1	42	
#	5	52	
#	10	55	
#	15	60	
#	20	64	
#	30	70	
#	40	74	
#	50	80	
#	60	82	
#	70	86	
#	80	89	
#	90	94	
#			
#	DP quantiles:		
#	1	5640	107848
#	5	6340	72504
#	10	6968	57655
#	15	7383	54062
#	20	7902	48405
#	30	9167	40279
#	40	9783	31257
#	50	10598	27042
#	60	11492	24714
#	70	12676	20907
#	80	13565	17948
#	90	14281	16869
#			
#	TP quantiles:		
#	1	1	36
#	5	1	36
#	10	2	35
#	15	3	34
#	20	4	33
#	30	5	32
#	40	7	30
#	50	9	29
#	60	10	28
#	70	11	26
#	80	16	25
#	90	17	20
#			
#	JOINT quantiles:		
#	1	0.005	
#	5	0.015	
#	10	0.03	
#	15	0.06	
#	20	0.075	
#	30	0.12	
#	40	0.18	
#	50	0.24	
#	60	0.32	
#	70	0.42	
#	80	0.56	
#	90	0.7	
#			
#	------------------------		
#	17.38%	at qual <1 (16.38% gain)	
#	33.82%	at qual <5 (28.82% gain)	
#	41.32%	at qual <10 (31.32% gain)	
#	52.82%	at qual <15 (37.82% gain)	
#	56.87%	at qual <20 (36.87% gain)	
#	62.65%	at qual <30 (32.65% gain)	
#	------------------------		


# Replicates-based recalibration showed best gain at 15, hence minQ 15:
# also, minimum genotype quality should be minGQ20, biallelic only, present in 50% of samples:	
~/scripts/vcftools --vcf denovo.recal.vcf --minQ 15 --minGQ 20 --min-alleles 2 --max-alleles 2 --max-missing 0.5 --recode --recode-INFO-all --out denovo.filt0
# OUTPUT::
# After filtering, kept 576 out of 576 Individuals
# After filtering, kept 6142 out of a possible 13254 Sites


# Discarding loci with too many heterozygotes, which are likely lumped paralogs
~/scripts/hetfilter.pl vcf=denovo.filt0.recode.vcf >denovo.hetfilt.vcf

# OUTPUT::
# 6142 total loci
# 0 dropped because fraction of missing genotypes exceeded 0.5
# 10 dropped because fraction of heterozygotes exceeded 0.75
# 6132 written

# Thinning SNP dataset - leaving one snp per tag. by default, it will leave the SNP with the highest minor allele frequency
~/scripts/thinner.pl vcf=denovo.hetfilt.vcf > thinDenov.vcf

# OUTPUT::
# 6132 total loci
# 0 loci skipped because they were closer than 40
# 5315 loci selected

# evaluating accuracy (the most important one is Heterozygote discovery rate, last column) based on replicates
~/scripts/repMatchStats.pl vcf=thinDenov.vcf replicates=clonepairs.tab 

# OUTPUT::
#	pair	gtyped	match	[ 00	1	11 ]	HetMatch	HomoHetMismatch	HetNoCall	HetsDiscoveryRate
#	SMID17a:SMID17b	5240	5013(95.7%)	 [65%	26%	9% ]	1315	111	6	0.96
#	SMID18a:SMID18b	5066	4770(94.2%)	 [65%	27%	9% ]	1281	118	25	0.95
#	SMID19a:SMID19b	5169	4092(79.2%)	 [62%	34%	4% ]	1396	864	52	0.75
#	SMID20a:SMID20b	5008	4784(95.5%)	 [65%	26%	9% ]	1237	163	17	0.93
#	SMID37a:SMID37b	5026	4706(93.6%)	 [64%	29%	7% ]	1361	225	16	0.92
#	SMID38a:SMID38b	5129	4744(92.5%)	 [64%	27%	8% ]	1301	221	16	0.92
#	SMID39a:SMID39b	4997	4499(90.0%)	 [63%	28%	8% ]	1267	279	20	0.89
#	SMID40a:SMID40b	4980	4633(93.0%)	 [64%	28%	8% ]	1282	234	18	0.91
#	SM2014RI16a:SM2014RI16b	4585	3799(82.9%)	 [70%	27%	3% ]	1011	354	93	0.82
#	SM2014RI17a:SM2014RI17b	4597	4023(87.5%)	 [71%	26%	3% ]	1057	384	67	0.82
#	SM2014RI18a:SM2014RI18b	4855	4049(83.4%)	 [71%	27%	3% ]	1089	451	66	0.81
#	SM2014RI19a:SM2014RI19b	4980	4298(86.3%)	 [71%	27%	2% ]	1150	409	30	0.84
#	SM2014SF17a:SM2014SF17b	5078	4298(84.6%)	 [71%	27%	2% ]	1152	574	43	0.79
#	SM2014SF18a:SM2014SF18b	4923	4193(85.2%)	 [71%	26%	3% ]	1101	396	47	0.83
#	SM2014SF19a:SM2014SF19b	5016	4255(84.8%)	 [71%	27%	2% ]	1154	630	28	0.78
#	SM2014SF20a:SM2014SF20b	4977	4445(89.3%)	 [72%	24%	4% ]	1047	219	26	0.9
#	SMAF17a:SMAF17b	5114	4577(89.5%)	 [74%	23%	3% ]	1055	391	12	0.84
#	SMAF18a:SMAF18b	4998	4212(84.3%)	 [72%	26%	2% ]	1080	370	57	0.83
#	SMAF19a:SMAF19b	4932	4244(86.1%)	 [72%	25%	3% ]	1082	513	34	0.8
#	SMAF20a:SMAF20b	4858	4290(88.3%)	 [72%	25%	2% ]	1093	506	33	0.8
#	SMAF37a:SMAF37b	5030	4708(93.6%)	 [74%	22%	4% ]	1040	193	13	0.91
#	SMAF38a:SMAF38b	4509	4047(89.8%)	 [72%	24%	3% ]	986	264	46	0.86
#	SMAF39a:SMAF39b	5035	4363(86.7%)	 [72%	25%	3% ]	1083	352	31	0.85
#	SMAF40a:SMAF40b	4948	4466(90.3%)	 [74%	22%	4% ]	978	202	34	0.89
#	SMFB17a:SMFB17b	4977	4620(92.8%)	 [74%	22%	5% ]	998	168	15	0.92
#	SMFB18a:SMFB18b	4979	4670(93.8%)	 [73%	23%	4% ]	1084	196	23	0.91
#	SMFB19a:SMFB19b	5182	4553(87.9%)	 [74%	22%	4% ]	1019	292	9	0.87
#	SMFB20a:SMFB20b	5125	4790(93.5%)	 [74%	22%	4% ]	1063	201	10	0.91
#	SMFB37a:SMFB37b	5137	4830(94.0%)	 [74%	22%	5% ]	1040	183	9	0.92
#	SMFB38a:SMFB38b	4907	4651(94.8%)	 [74%	22%	4% ]	1030	199	20	0.9
#	SMFB39a:SMFB39b	5004	4643(92.8%)	 [73%	23%	4% ]	1080	310	14	0.87
#	SMFB40a:SMFB40b	4906	4612(94.0%)	 [73%	22%	4% ]	1026	213	29	0.89
#	SMKB17a:SMKB17b	5149	4893(95.0%)	 [74%	21%	5% ]	1035	183	7	0.92
#	SMKB18a:SMKB18b	5021	4709(93.8%)	 [73%	22%	4% ]	1047	183	12	0.91
#	SMKB19a:SMKB19b	5114	4818(94.2%)	 [75%	21%	5% ]	1008	228	7	0.9
#	SMKB20a:SMKB20b	5097	4867(95.5%)	 [73%	22%	5% ]	1067	148	14	0.93
#	SMKB37a:SMKB37b	5189	4521(87.1%)	 [72%	22%	5% ]	1015	213	25	0.9
#	SMKB38a:SMKB38b	5219	4880(93.5%)	 [73%	22%	5% ]	1096	159	3	0.93
#	SMKB39a:SMKB39b	5165	4792(92.8%)	 [74%	22%	4% ]	1048	239	6	0.9
#	SMKB40a:SMKB40b	5183	4300(83.0%)	 [72%	25%	3% ]	1064	413	26	0.83
#	SMLB17a:SMLB17b	5220	4631(88.7%)	 [74%	22%	5% ]	1009	285	20	0.87
#	SMLB18a:SMLB18b	4554	4182(91.8%)	 [73%	23%	4% ]	969	211	49	0.88
#	SMLB19a:SMLB19b	5013	4372(87.2%)	 [73%	25%	3% ]	1073	419	23	0.83
#	SMLB20a:SMLB20b	5079	4658(91.7%)	 [74%	22%	4% ]	1023	318	19	0.86
#	SMLB37a:SMLB37b	5201	4758(91.5%)	 [73%	24%	3% ]	1136	438	2	0.84
#	SMLB38a:SMLB38b	5131	4697(91.5%)	 [74%	22%	4% ]	1019	310	13	0.86
#	SMLB39a:SMLB39b	5097	4652(91.3%)	 [74%	23%	3% ]	1072	344	17	0.86
#	SMLB40a:SMLB40b	5046	4707(93.3%)	 [73%	23%	4% ]	1084	204	16	0.91
#	SMLYF17a:SMLYF17b	5064	4753(93.9%)	 [74%	22%	4% ]	1038	248	9	0.89
#	SMLYF18a:SMLYF18b	4776	4469(93.6%)	 [73%	23%	4% ]	1021	277	31	0.87
#	SMLYF19a:SMLYF19b	5034	4134(82.1%)	 [72%	25%	3% ]	1047	352	46	0.84
#	SMLYF20a:SMLYF20b	5018	4719(94.0%)	 [73%	23%	4% ]	1083	251	11	0.89
#	SMLYF37a:SMLYF37b	4194	3962(94.5%)	 [71%	24%	5% ]	965	160	72	0.89
#	SMLYF38a:SMLYF38b	5042	4587(91.0%)	 [73%	23%	4% ]	1042	295	21	0.87
#	SMLYF39a:SMLYF39b	4913	4355(88.6%)	 [73%	24%	3% ]	1033	358	24	0.84
#	SMLYF40a:SMLYF40b	4920	4050(82.3%)	 [71%	25%	4% ]	1019	217	56	0.88
#	SMOS17a:SMOS17b	5029	4379(87.1%)	 [73%	23%	3% ]	1029	408	40	0.82
#	SMOS18a:SMOS18b	5143	4391(85.4%)	 [73%	24%	3% ]	1044	403	52	0.82
#	SMOS19a:SMOS19b	5050	4421(87.5%)	 [74%	22%	4% ]	985	373	41	0.83
#	SMOS20a:SMOS20b	5108	4477(87.6%)	 [75%	22%	3% ]	999	418	29	0.82
#	SMOS37a:SMOS37b	5211	4764(91.4%)	 [73%	24%	3% ]	1126	445	4	0.83
#	SMOS38a:SMOS38b	5223	4640(88.8%)	 [72%	25%	3% ]	1164	555	7	0.81
#	SMOS39a:SMOS39b	5194	4694(90.4%)	 [74%	23%	3% ]	1101	496	2	0.82
#	SMOS40a:SMOS40b	5237	4743(90.6%)	 [74%	23%	3% ]	1099	477	3	0.82
#	SMRI17a:SMRI17b	4983	4735(95.0%)	 [74%	22%	4% ]	1055	220	15	0.9
#	SMRI18a:SMRI18b	4913	4637(94.4%)	 [74%	22%	4% ]	1018	240	16	0.89
#	SMRI19a:SMRI19b	5193	4904(94.4%)	 [74%	21%	5% ]	1046	218	6	0.9
#	SMRI20a:SMRI20b	4971	4738(95.3%)	 [75%	21%	5% ]	976	183	8	0.91
#	SMRI37a:SMRI37b	5184	4894(94.4%)	 [74%	21%	5% ]	1024	227	5	0.9
#	SMRI38a:SMRI38b	5089	4582(90.0%)	 [73%	23%	4% ]	1045	296	11	0.87
#	SMRI39a:SMRI39b	5119	4720(92.2%)	 [74%	22%	4% ]	1046	385	8	0.84
#	SMRI40a:SMRI40b	5103	4754(93.2%)	 [74%	22%	4% ]	1045	320	5	0.87
#	SMSF17a:SMSF17b	5109	4737(92.7%)	 [74%	22%	5% ]	1028	199	16	0.91
#	SMSF18a:SMSF18b	5046	4713(93.4%)	 [74%	22%	4% ]	1029	202	12	0.91
#	SMSF19a:SMSF19b	4947	4434(89.6%)	 [74%	22%	4% ]	989	267	24	0.87
#	SMSF20a:SMSF20b	5032	4725(93.9%)	 [74%	21%	4% ]	1008	235	16	0.89
#	SMSF37a:SMSF37b	4959	4468(90.1%)	 [74%	23%	3% ]	1021	366	24	0.84
#	SMSF38a:SMSF38b	4915	4479(91.1%)	 [73%	23%	4% ]	1016	252	24	0.88
#	SMSF39a:SMSF39b	5072	4636(91.4%)	 [73%	22%	4% ]	1040	258	11	0.89
#	SMSF40a:SMSF40b	4871	4529(93.0%)	 [74%	22%	4% ]	993	253	17	0.88
#	SMST17a:SMST17b	5150	4677(90.8%)	 [73%	23%	3% ]	1082	366	12	0.85
#	SMST18a:SMST18b	5072	4571(90.1%)	 [74%	22%	4% ]	996	281	22	0.87
#	SMST19a:SMST19b	5022	4604(91.7%)	 [74%	23%	4% ]	1039	312	15	0.86
#	SMST20a:SMST20b	4994	4217(84.4%)	 [73%	26%	2% ]	1087	664	34	0.76
#	SMST37a:SMST37b	5040	4541(90.1%)	 [73%	24%	3% ]	1070	421	33	0.82
#	SMST38a:SMST38b	5142	4593(89.3%)	 [74%	22%	4% ]	1027	268	22	0.88
#	SMST39a:SMST39b	5135	4656(90.7%)	 [75%	22%	3% ]	1027	451	15	0.82
#	SMST40a:SMST40b	4959	4707(94.9%)	 [73%	22%	5% ]	1038	250	14	0.89
#	SMTF17a:SMTF17b	5152	4804(93.2%)	 [74%	21%	5% ]	994	208	12	0.9
#	SMTF18a:SMTF18b	5015	4733(94.4%)	 [74%	22%	5% ]	1029	109	17	0.94
#	SMTF19a:SMTF19b	5065	4761(94.0%)	 [75%	21%	4% ]	990	251	6	0.89
#	SMTF20a:SMTF20b	5029	4819(95.8%)	 [73%	22%	5% ]	1074	168	18	0.92
#	SMTF37a:SMTF37b	5189	4619(89.0%)	 [72%	26%	2% ]	1210	527	6	0.82
#	SMTF38a:SMTF38b	4970	4606(92.7%)	 [74%	23%	4% ]	1047	304	23	0.86
#	SMTF39a:SMTF39b	5156	4634(89.9%)	 [72%	24%	4% ]	1108	346	17	0.86
#	SMTF40a:SMTF40b	5161	4548(88.1%)	 [74%	22%	4% ]	1006	359	22	0.84
#										
#	------------------------									
#	hets called homos depth: 									
#	lower 25%	29								
#	median		51							
#	upper 75%	80								


# creating a file of replicates and other poorly sequenced samples (low sites number/high-homozygosity)
~/scripts/vcftools --vcf thinDenov.vcf --het --out denovo.out

#OPEN denovo.out.het AND NOTE ALL INDIVIDUALS WITH LOW N_SITES, ADD THESE TO clones2remove
nano clones2remove
#edit accordingly 
#	SM2014RI16b
#	SM2014RI17b
#	SM2014RI18b
#	SM2014RI19b
#	SM2014SF17b
#	SM2014SF18b
#	SM2014SF19b
#	SM2014SF20b
#	SMAF17b
#	SMAF18b
#	SMAF19b
#	SMAF20b
#	SMAF37b
#	SMAF38b
#	SMAF39b
#	SMAF40b
#	SMFB17b
#	SMFB18b
#	SMFB19b
#	SMFB20b
#	SMFB37b
#	SMFB38b
#	SMFB39b
#	SMFB40b
#	SMID17b
#	SMID18b
#	SMID19b
#	SMID20b
#	SMID37b
#	SMID38b
#	SMID39b
#	SMID40b
#	SMKB17b
#	SMKB18b
#	SMKB19b
#	SMKB20b
#	SMKB37b
#	SMKB38b
#	SMKB39b
#	SMKB40b
#	SMLB17b
#	SMLB18b
#	SMLB19b
#	SMLB20b
#	SMLB37b
#	SMLB38b
#	SMLB39b
#	SMLB40b
#	SMLYF17b
#	SMLYF18b
#	SMLYF19b
#	SMLYF20b
#	SMLYF37b
#	SMLYF38b
#	SMLYF39b
#	SMLYF40b
#	SMOS17b
#	SMOS18b
#	SMOS19b
#	SMOS20b
#	SMOS37b
#	SMOS38b
#	SMOS39b
#	SMOS40b
#	SMRI17b
#	SMRI18b
#	SMRI19b
#	SMRI20b
#	SMRI37b
#	SMRI38b
#	SMRI39b
#	SMRI40b
#	SMSF17b
#	SMSF18b
#	SMSF19b
#	SMSF20b
#	SMSF37b
#	SMSF38b
#	SMSF39b
#	SMSF40b
#	SMST17b
#	SMST18b
#	SMST19b
#	SMST20b
#	SMST37b
#	SMST38b
#	SMST39b
#	SMST40b
#	SMTF17b
#	SMTF18b
#	SMTF19b
#	SMTF20b
#	SMTF37b
#	SMTF38b
#	SMTF39b
#	SMTF40b

# creating thinned dataset without clones and poorly covered individuals
~/scripts/vcftools --vcf thinDenov.vcf --remove clones2remove --recode --recode-INFO-all --out skg
# OUTPUT
# Excluding individuals in 'exclude' list
# After filtering, kept 480 out of 576 Individuals
# After filtering, kept 5315 out of a possible 5315 Sites


# additional filtering on final datasets max misson 0.7 then minor allele frequency 0.01 
~/scripts/vcftools --vcf skg.recode.vcf --max-missing 0.7 --recode --recode-INFO-all --out skg.filt
# OUTPUT
# After filtering, kept 5132 out of a possible 5315 Sites
~/scripts/vcftools --vcf skg.filt.recode.vcf --maf 0.01 --recode --recode-INFO-all --out skg.filt.maf0.01
# OUTPUT
# After filtering, kept 4123 out of a possible 5132 Sites

