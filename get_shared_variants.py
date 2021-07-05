#python3

# script to read in VCF and manifest file detailing pedigrees, annotate variants by GnoMAD frequency, functional effect and sharing status
# Ben Kinnersley (ben.kinnersley@icr.ac.uk)
# python3 get_shared_variants.py <INPUT_VCF> <INPUT_MANIFEST> <OUTPUT_FILE> <OUTPUT_RECESSIVE>

import gzip
import sys
import os

if len(sys.argv) == 5:
	input_vcf = sys.argv[1]
	input_manifest = sys.argv[2]
	output_file = sys.argv[3]
	output_recessive = sys.argv[4]
else:
	print('incorrect arguments provided: python3 get_shared_variants.py <INPUT_VCF> <INPUT_MANIFEST> <OUTPUT_FILE> <OUTPUT_RECESSIVE>')
	sys.exit()

# set thresholds for outputting variants (todo: add as command line options)
gq_threshold = 30
depth_threshold = 5
GnoMAD_NFE_threshold = 0.005
family_sharing_flag = "AT_LEAST_TWO_AFFECTED" # "ALL_AFFECTED" - only output variants shared in all affecteds in a family, "AT_LEAST_TWO_AFFECTED" = at least two affected shared

if os.path.isfile(input_vcf):
	if input_vcf.endswith('gz'):
		opened_input_vcf = gzip.open(input_vcf, 'rt')
	else:
		opened_input_vcf = open(input_vcf)
	print('reading query variants from'+input_vcf)
else:
	print('could not open file '+input_vcf)
	sys.exit()

if os.path.isfile(input_manifest):
	opened_input_manifest = open(input_manifest)
	opened_input_manifest.readline()
	
	print('reading families from '+input_manifest)
else:
	print('could not open file '+input_manifest)
	sys.exit()
	
opened_output_file = open(output_file, 'w')

print('writing to output file '+output_file)

print('now getting families to query...')

total_vcf_id_list = []
total_family_id_list = []
total_affected_vcf_id_list = []
total_unaffected_vcf_id_list = []
all_carriers_per_family_dict = {}
affected_carriers_per_family_dict = {}
unaffected_carriers_per_family_dict = {}
family_id_dict = {}

for line in opened_input_manifest:
	fields = line.split('\t')
	vcf_id = fields[0].strip()
	family_id = fields[1].strip()
	affected_flag = fields[2].strip()
	
	family_id_dict[vcf_id] = family_id
	
	if family_id not in total_family_id_list:
		total_family_id_list.append(family_id)
	
	total_vcf_id_list.append(vcf_id)
	
	if affected_flag == "1":
		total_affected_vcf_id_list.append(vcf_id)
		
		if family_id in affected_carriers_per_family_dict:
			affected_carriers_per_family_dict[family_id] = affected_carriers_per_family_dict[family_id]+'/'+vcf_id
		else:
			affected_carriers_per_family_dict[family_id] = vcf_id
			
	elif affected_flag == "0":
		total_unaffected_vcf_id_list.append(vcf_id)
		
		if family_id in unaffected_carriers_per_family_dict:
			unaffected_carriers_per_family_dict[family_id] = unaffected_carriers_per_family_dict[family_id]+'/'+vcf_id
		else:
			unaffected_carriers_per_family_dict[family_id] = vcf_id
		
	if family_id in all_carriers_per_family_dict:
		all_carriers_per_family_dict[family_id] = all_carriers_per_family_dict[family_id]+'/'+vcf_id
	else:
		all_carriers_per_family_dict[family_id] = vcf_id

for family in total_family_id_list:
	if family in affected_carriers_per_family_dict:
		affected_list = affected_carriers_per_family_dict[family].split('/')
	else:
		affected_list = ''
	if family in unaffected_carriers_per_family_dict:
		unaffected_list = unaffected_carriers_per_family_dict[family].split('/')
	else:
		unaffected_list = ''
	print('processing family '+str(family)+' with '+str(len(affected_list))+' affected and '+str(len(unaffected_list))+' unaffected individuals')

# print output header
opened_output_file.write('chr\tpos\tref\talt\tcanonical_flag\tconsequence_canon\tIMPACT_canon\tGene_canon\tSYMBOL_canon\tHGVSc_canon\tHGVSp_canon\tCADD_PHRED\t')
opened_output_file.write('gnomAD_AF_NFE\tAF_1kG\tEUR_AF_1kG\tMAX_AF_1kG\tClinVar')
opened_output_file.write('\tExisting_variation\tshared_families\tcarrier_pass\tcarrier_filtered\tnoncarrier_pass\tnoncarrier_filtered\tno_call\ttotal_affected_count')
for family in total_family_id_list:
	opened_output_file.write('\t'+family+'_count')
opened_output_file.write('\n')
			
print('now reading in variants from input vcf...')

chromosome_list = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrMT']

consequence_severity_dict = { "transcript_ablation":1, "splice_acceptor_variant":2, "splice_donor_variant":3, "stop_gained":4, "frameshift_variant":5,
"stop_lost":6, "start_lost":7, "transcript_amplification":8, "inframe_insertion":9, "inframe_deletion":10, "missense_variant":11,
"protein_altering_variant":12, "splice_region_variant":13, "incomplete_terminal_codon_variant":14, "start_retained_variant":15,
"stop_retained_variant":16, "synonymous_variant":17, "coding_sequence_variant":18, "mature_miRNA_variant":19, "5_prime_UTR_variant":20,
"3_prime_UTR_variant":21, "non_coding_transcript_exon_variant":22, "intron_variant":23, "NMD_transcript_variant":24,
"non_coding_transcript_variant":25, "upstream_gene_variant":26, "downstream_gene_variant":27, "TFBS_ablation":28, "TFBS_amplification":29,
"TF_binding_site_variant":30, "regulatory_region_ablation":31, "regulatory_region_amplification":32, "feature_elongation":33,
"regulatory_region_variant":34, "feature_truncation":35, "intergenic_variant":36
}

sample_dict = {}
vep_lookup_hash = {}

for line in opened_input_vcf:
	fields = line.split('\t')
	
	# get INFO fields for VEP and read in order
	if line.startswith('##INFO=<ID=CSQ'):
		temp1 = line.split(':')
		temp2 = temp1[1].strip().split('"')
		vep_fields = temp2[0].strip().split('|')
		
		vep_ticker = 0
		
		for anno in vep_fields:
			vep_lookup_hash[anno] = vep_ticker
			vep_ticker += 1
			print(anno)
	
	# parse vcf header to obtain sample list
	
	if line.startswith('#CHROM'):
		vcf_sample_header = line.split('\t')
		vcf_samples_list = vcf_sample_header[9:]
		
		header_ticker = 0
		
		for sample in vcf_samples_list:
			if "variant" in sample:
				sample_split = sample.split(".")
				sample_out = sample_split[0]
			else:
				sample_out = sample
			
			print(sample_out)
			sample_dict[sample_out.strip()] = header_ticker
			header_ticker += 1
		
	if fields[0] in chromosome_list or 'chr'+fields[0] in chromosome_list:
		chr = fields[0].strip()
		pos = fields[1].strip()
		id = fields[2].strip()
		ref = fields[3].strip()
		alt = fields[4].strip()
		qual = fields[5].strip()
		filter = fields[6].strip()
		info = fields[7].strip()
		format = fields[8].strip()
		
		# skip variants that are not "PASS" or "."
		if filter == "PASS" or filter == ".":
			pass
		else:
			continue
		
		#print(line)
		
		var = chr+':'+pos+'_'+ref+'/'+alt
		alt_split = alt.split(',')
		alt_list = []
		
		for alt in alt_split:
			if alt != "*":
				alt_list.append(alt)
		
		# parse variant FORMAT fields
		format_fields = format.split(':')
		format_ticker = 0
		for query in format_fields:
			if query == "GT":
				genotype = format_ticker
			if query == "AD":
				allelic_depth = format_ticker
			if query == "DP":
				total_depth = format_ticker
			if query == "GQ":
				genotype_quality = format_ticker
			format_ticker += 1
		#GT:AD:DP:GQ:MIN_DP:PGT:PID:PL:PS:RGQ:SB		
		#GT:AD:DP:GQ:PL
		
		# parse VEP fields
		# note: as --allele_number flag specified in VEP command we can assume for muiltiallelic SNPs alleles are ordered by their order in alt column (https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html#opt_allele_number)
		# note: this avoids having to convert alleles to ensembl format to check against each VEP string (https://www.ensembl.org/info/docs/tools/vep/vep_formats.html)
		info_split = info.split(';')
		#print(line)
		for query in info_split:
			if query.startswith('CSQ='):
				vep_split = query.split(',')
		
		if alt == "*":
			continue
		
		# create list of alleles to iterate over in vep_split (if we know order of alleles and length of vep_split)
		vep_allele_list = []
		for i in range(0, int(float(len(vep_split))/float(len(alt_list)))):
			for alt in alt_list:
				#print(alt)
				vep_allele_list.append(alt)
		#print(vep_allele_list)
		
		canonical_flag_dict = {}
		consequence_canon_dict = {}
		IMPACT_canon_dict = {}
		Gene_canon_dict = {}
		SYMBOL_canon_dict = {}
		HGVSc_canon_dict = {}
		HGVSp_canon_dict = {}
		CADD_PHRED_dict = {}
		gnomADg_non_cancer_AF_nfe_dict = {}
		ClinVar_dict = {}
		Existing_variation_dict = {}
		gnomAD_AF_dict = {}
		MAX_AF_dict = {}
		EUR_AF_dict = {}
		AF_dict = {}
		
		vep_allele_ticker = 0
		
		for allele in vep_allele_list:
			vep_string = vep_split[vep_allele_ticker]
			
			vep_string_split = vep_string.split('|')
			
			canonical_flag = vep_string_split[vep_lookup_hash['CANONICAL']]
			consequence = vep_string_split[vep_lookup_hash['Consequence']]
			IMPACT = vep_string_split[vep_lookup_hash['IMPACT']]
			
			if 'CADD_PHRED' in vep_lookup_hash:
				CADD_PHRED = vep_string_split[vep_lookup_hash['CADD_PHRED']]
			else:
				CADD_PHRED = 'NA'
			#gnomADg_non_cancer_AF_nfe = vep_string_split[vep_lookup_hash['gnomADg_non_cancer_AF_nfe']]
			
			if 'ClinVar' in vep_lookup_hash:
				ClinVar = vep_string_split[vep_lookup_hash['ClinVar']]
			elif 'CLIN_SIG' in vep_lookup_hash:
				ClinVar = vep_string_split[vep_lookup_hash['CLIN_SIG']]
			else:
				ClinVar = 'NA'
				
			Existing_variation = vep_string_split[vep_lookup_hash['Existing_variation']]
			Gene = vep_string_split[vep_lookup_hash['Gene']]
			SYMBOL = vep_string_split[vep_lookup_hash['SYMBOL']]
			HGVSc = vep_string_split[vep_lookup_hash['HGVSc']]
			HGVSp = vep_string_split[vep_lookup_hash['HGVSp']]
			
			if 'gnomADg_AF_nfe' in vep_lookup_hash:
				gnomAD_AF = vep_string_split[vep_lookup_hash['gnomADg_AF_nfe']]
			elif 'gnomAD_NFE_AF' in vep_lookup_hash:
				gnomAD_AF = vep_string_split[vep_lookup_hash['gnomAD_NFE_AF']]
			else:
				gnomAD_AF = 'NA'
				
			MAX_AF = vep_string_split[vep_lookup_hash['MAX_AF']]
			EUR_AF = vep_string_split[vep_lookup_hash['EUR_AF']]
			AF_1kG = vep_string_split[vep_lookup_hash['AF']]
			#print(gnomAD_AF)
			
			#gnomADg_non_cancer_AF_nfe_split = gnomADg_non_cancer_AF_nfe.split("&")
			#
			#gnomADg_non_cancer_AF_nfe = gnomADg_non_cancer_AF_nfe_split[0]
			#
			## to detect instances where gnomad is returning two rsids/AFs for a given variant (e.g. if liftover from two different build37 positions to 1 build38 position)
			#if len(gnomADg_non_cancer_AF_nfe_split) >= 2:
			#	print('multiple gnomAD matches for '+str(var))
			#
			## in case of multiple frequencies (due to liftover issue...), return highest
			#for gnomad in gnomADg_non_cancer_AF_nfe_split:
			#	#print(gnomad)
			#	if gnomad != "" and float(gnomad) > float(gnomADg_non_cancer_AF_nfe):
			#		gnomADg_non_cancer_AF_nfe = gnomad
			
			MAX_AF_split = MAX_AF.split("&")
			
			MAX_AF = MAX_AF_split[0]
			
			# to detect instances where AF_1kG is returning two rsids/AFs for a given variant
			if len(MAX_AF_split) >= 2:
				print('multiple MAX_AF matches for '+str(var))
			
			# in case of multiple frequencies, return highest
			for AF in MAX_AF_split:
				#print(gnomad)
				if AF != "" and float(AF) > float(MAX_AF):
					MAX_AF = AF
					
			EUR_AF_split = EUR_AF.split("&")
			
			EUR_AF = EUR_AF_split[0]
			
			# to detect instances where AF_1kG is returning two rsids/AFs for a given variant
			if len(EUR_AF_split) >= 2:
				print('multiple EUR_AF matches for '+str(var))
			
			# in case of multiple frequencies, return highest
			for AF in EUR_AF_split:
				#print(gnomad)
				if AF != "" and float(AF) > float(EUR_AF):
					EUR_AF = AF
					
			AF_split = AF_1kG.split("&")
			
			AF_1kG = AF_split[0]
			
			# to detect instances where AF_1kG is returning two rsids/AFs for a given variant
			if len(AF_split) >= 2:
				print('multiple AF matches for '+str(var))
			
			# in case of multiple frequencies, return highest
			for AF in AF_split:
				#print(gnomad)
				if AF != "" and float(AF) > float(AF_1kG):
					AF_1kG = AF
				
			if allele not in consequence_canon_dict:
				#canonical_flag_dict[allele] = 'NA'
				consequence_canon_dict[allele] = 'NA'
				IMPACT_canon_dict[allele] = 'NA'
				Gene_canon_dict[allele] = 'NA'
				SYMBOL_canon_dict[allele] = 'NA'
				HGVSc_canon_dict[allele] = 'NA'
				HGVSp_canon_dict[allele] = 'NA'
			
			consequence_split = consequence.split('&')
			
			# checking consequence severity for canonical transcript, retaining annotations from most severe canonical consequence
			
			for consequence in consequence_split:
				if canonical_flag == "YES":
					if allele not in consequence_canon_dict or consequence_canon_dict[allele] == "NA":
						consequence_canon_dict[allele] = consequence
						IMPACT_canon_dict[allele] = IMPACT
						Gene_canon_dict[allele] = Gene
						SYMBOL_canon_dict[allele] = SYMBOL
						HGVSc_canon_dict[allele] = HGVSc
						HGVSp_canon_dict[allele] = HGVSp
						
					elif consequence_canon_dict[allele] != "NA" and int(consequence_severity_dict[consequence]) < int(consequence_severity_dict[consequence_canon_dict[allele]]):
						consequence_canon_dict[allele] = consequence
						IMPACT_canon_dict[allele] = IMPACT
						Gene_canon_dict[allele] = Gene
						SYMBOL_canon_dict[allele] = SYMBOL
						HGVSc_canon_dict[allele] = HGVSc
						HGVSp_canon_dict[allele] = HGVSp	
			
			CADD_PHRED_dict[allele] = CADD_PHRED
			#gnomADg_non_cancer_AF_nfe_dict[allele] = gnomADg_non_cancer_AF_nfe
			gnomAD_AF_dict[allele] = gnomAD_AF
			ClinVar_dict[allele] = ClinVar
			Existing_variation_dict[allele] = Existing_variation
			MAX_AF_dict[allele] = MAX_AF
			EUR_AF_dict[allele] = EUR_AF
			AF_dict[allele] = AF_1kG
			
			vep_allele_ticker += 1
		
		alt_counter = 0
		
		for allele in alt_list:
			alt_counter += 1
			# parse genotypes
			geno_samples_list = fields[9:]
			#print(var)
			
			carrier_pass = ''
			carrier_filtered = ''
			noncarrier_pass = ''
			noncarrier_filtered = ''
			no_call = ''
			family_count_dict = {}
			total_affected_count = 0
			
			for family in total_family_id_list:
				#print(family)
				family_count_dict[family] = 0
				#print(family_count_dict[family])
				
			for sample in total_vcf_id_list:
				sample.strip()
				
				sample_geno = geno_samples_list[sample_dict[sample.strip()]]
				
				sample_geno_split = sample_geno.split(':')
	
				if len(sample_geno_split) <= genotype_quality + 1: # sanity check to ignore genotypes that don't have format fields up to GQ
					#print(sample_geno)
					continue
					
				sample_genotype = sample_geno_split[genotype]
				sample_allelic_depth = sample_geno_split[allelic_depth]
				sample_total_depth = sample_geno_split[total_depth]
				sample_genotype_quality = sample_geno_split[genotype_quality]
				
				geno_filter = ''
				
							
				if sample_genotype_quality == "." or sample_total_depth == "." or int(sample_genotype_quality) < int(gq_threshold) or int(sample_total_depth) < int(depth_threshold):
					geno_filter = 'YES'
				
				if sample_genotype == "./." or sample_genotype == ".|.":
					if no_call == '':
						no_call = sample+':'+sample_genotype+':'+sample_allelic_depth+':'+sample_total_depth+':'+sample_genotype_quality
					else:
						no_call = no_call+';'+sample+':'+sample_genotype+':'+sample_allelic_depth+':'+sample_total_depth+':'+sample_genotype_quality
				elif sample_genotype == "0/0" or sample_genotype == "0|0":
					if geno_filter == "YES":
						if noncarrier_filtered == '':
							noncarrier_filtered = sample+':'+sample_genotype+':'+sample_allelic_depth+':'+sample_total_depth+':'+sample_genotype_quality
						else:
							noncarrier_filtered = noncarrier_filtered+';'+sample+':'+sample_genotype+':'+sample_allelic_depth+':'+sample_total_depth+':'+sample_genotype_quality
					else:
						if noncarrier_pass == '':
							noncarrier_pass = sample+':'+sample_genotype+':'+sample_allelic_depth+':'+sample_total_depth+':'+sample_genotype_quality
						else:
							noncarrier_pass = noncarrier_pass+';'+sample+':'+sample_genotype+':'+sample_allelic_depth+':'+sample_total_depth+':'+sample_genotype_quality
				else:
					# only consider non-ref genotypes with at least one allele matching alt under consideration
					if "|" in sample_genotype:
						sample_genotype_split = sample_genotype.split('|')
						a1 = sample_genotype_split[0]
						a2 = sample_genotype_split[1]
					elif "/" in sample_genotype:
						sample_genotype_split = sample_genotype.split('/')
						a1 = sample_genotype_split[0]
						a2 = sample_genotype_split[1]
					else:
						print('unable to determine delimiter of geno '+sample_genotype)
						geno_filter = "YES"
						
					#print('considering alt alleles '+str(a1)+' and '+str(a2)+' for geno '+str(sample_genotype)+' and counter '+str(alt_counter))
					
					if int(a1) != int(alt_counter) and int(a2) != int(alt_counter):
						print('ignoring alleles '+str(a1)+' and '+str(a2)+' for geno '+str(sample_genotype)+' and counter '+str(alt_counter))
						continue
				
					if geno_filter == "YES":
						if carrier_filtered == '':
							carrier_filtered = sample+':'+sample_genotype+':'+sample_allelic_depth+':'+sample_total_depth+':'+sample_genotype_quality
						else:
							carrier_filtered = carrier_filtered+';'+sample+':'+sample_genotype+':'+sample_allelic_depth+':'+sample_total_depth+':'+sample_genotype_quality
					else:
						if carrier_pass == '':
							carrier_pass = sample+':'+sample_genotype+':'+sample_allelic_depth+':'+sample_total_depth+':'+sample_genotype_quality
						else:
							carrier_pass = carrier_pass+';'+sample+':'+sample_genotype+':'+sample_allelic_depth+':'+sample_total_depth+':'+sample_genotype_quality
							
						if sample in total_affected_vcf_id_list:
							total_affected_count += 1
							family_count_dict[family_id_dict[sample]] += 1
							#print(family_count_dict[family_id_dict[sample]])
							#print(total_affected_count)
				
			#print('no_calls: '+no_call)
			#print('noncarrier_filtered: '+noncarrier_filtered)
			#print('noncarrier_pass: '+noncarrier_pass)
			#print('carrier_filtered: '+carrier_filtered)
			#print('carrier_pass: '+carrier_pass)
			
			segregating_family_str = ''
			
			for family_id in total_family_id_list:
				affected_list = affected_carriers_per_family_dict[family_id].split('/')
			
				#print(family_count_dict[family_id])
				#print(len(affected_list))
				
				if family_id in family_count_dict and family_count_dict[family_id] == len(affected_list) and len(affected_list) > 1:
					if segregating_family_str == '':
						segregating_family_str = family_id
					else:
						segregating_family_str = segregating_family_str+'/'+family_id
			
			# write to output
		
			#print(gnomADg_non_cancer_AF_nfe_dict[allele])
			#print(GnoMAD_NFE_threshold)
			
			#print(str(chr)+'\t'+str(pos)+'\t'+str(ref)+'\t'+str(alt))
			
			if allele not in canonical_flag_dict:
				canonical_flag_dict[allele] = 'NA'
			
			if family_sharing_flag == "ALL_AFFECTED" and len(segregating_family_str) < 1:
				continue
			
			if gnomAD_AF_dict[allele] != "" and float(gnomAD_AF_dict[allele]) > float(GnoMAD_NFE_threshold):
				continue
				
			if MAX_AF_dict[allele] != "" and float(MAX_AF_dict[allele]) > float(GnoMAD_NFE_threshold):
				continue
			
			if EUR_AF_dict[allele] != "" and float(EUR_AF_dict[allele]) > float(GnoMAD_NFE_threshold):
				continue
				
			if AF_dict[allele] != "" and float(AF_dict[allele]) > float(GnoMAD_NFE_threshold):
				continue
				
			if int(total_affected_count) < 2:
				continue
				
			if IMPACT_canon_dict[allele] != "HIGH" and IMPACT_canon_dict[allele] != "MODERATE":
				continue
			
			#print(line)
			opened_output_file.write(str(chr)+'\t'+str(pos)+'\t'+str(ref)+'\t'+str(allele)+'\t'+str(canonical_flag_dict[allele])+'\t'+str(consequence_canon_dict[allele])+'\t'+str(IMPACT_canon_dict[allele]))
			opened_output_file.write('\t'+str(Gene_canon_dict[allele])+'\t'+str(SYMBOL_canon_dict[allele])+'\t'+str(HGVSc_canon_dict[allele])+'\t'+str(HGVSp_canon_dict[allele]))
			opened_output_file.write('\t'+str(CADD_PHRED_dict[allele])+'\t'+str(gnomAD_AF_dict[allele])+'\t'+str(AF_dict[allele])+'\t'+str(EUR_AF_dict[allele])+'\t'+str(MAX_AF_dict[allele]))
			opened_output_file.write('\t'+str(ClinVar_dict[allele])+'\t'+str(Existing_variation_dict[allele])+'\t'+str(segregating_family_str))
			opened_output_file.write('\t'+str(carrier_pass)+'\t'+str(carrier_filtered)+'\t'+str(noncarrier_pass)+'\t'+str(noncarrier_filtered)+'\t'+str(no_call))
			opened_output_file.write('\t'+str(total_affected_count))
			
			for family in total_family_id_list:
				opened_output_file.write('\t'+str(family_count_dict[family]))
			opened_output_file.write('\n')

opened_output_file.close()

print('now reviewing for potential shared recessive variants and compound heterozygotes in the same gene')

opened_output_recessive = open(output_recessive, 'w')

print('writing to output file '+output_recessive)

opened_output_as_input = open(output_file)

header = opened_output_as_input.readline()

opened_output_recessive.write(header)

recessive_gene_list = []
recessive_sample_list = []
gene_carrier_count_dict = {}
gene_carrier_variant_dict = {}
variant_dict = {}

for line in opened_output_as_input:
	line = line.strip()
	fields = line.split('\t')
	
	ensembl_gene = fields[7]
	gene_symbol = fields[8]
	carrier_pass = fields[19]
	
	#carrier_pass = carrier_pass+';'+sample+'_'+sample_genotype+'_'+sample_allelic_depth+'_'+sample_total_depth+'_'+sample_genotype_quality
	carrier_split = carrier_pass.split(';')
	
	if ensembl_gene not in gene_carrier_count_dict:
		gene_carrier_count_dict[ensembl_gene] = {}
	
	if ensembl_gene not in gene_carrier_variant_dict:
		gene_carrier_variant_dict[ensembl_gene] = {}
	
	for carrier in carrier_split:
		fields = carrier.split(':')
		sample = fields[0]
		sample_genotype = fields[1]
		
		# we are considering "recessive" carriers as those with 2 alleles in a given gene, either for the same variant ("classifical" bi-allelic recessive inheriance) or two alleles from different variants in the same gene (e.g. compound heterozygosity)
		
		if sample not in gene_carrier_count_dict[ensembl_gene]:
			gene_carrier_count_dict[ensembl_gene][sample] = 0
		
		if "|" in sample_genotype:
			sample_genotype_split = sample_genotype.split('|')
			a1 = sample_genotype_split[0]
			a2 = sample_genotype_split[1]
		elif "/" in sample_genotype:
			sample_genotype_split = sample_genotype.split('/')
			a1 = sample_genotype_split[0]
			a2 = sample_genotype_split[1]
		
		if int(a1) == int(a2):
			# 2 copies of allele
			gene_carrier_count_dict[ensembl_gene][sample] += 2
		else:
			# 1 copy of allele
			gene_carrier_count_dict[ensembl_gene][sample] += 1
		
		if sample not in gene_carrier_variant_dict[ensembl_gene]:
			gene_carrier_variant_dict[ensembl_gene][sample] = []
		
		gene_carrier_variant_dict[ensembl_gene][sample].append(line)
		#print(gene_carrier_variant_dict[ensembl_gene][sample])
	
		if gene_carrier_count_dict[ensembl_gene][sample] >= 2 and ensembl_gene not in recessive_gene_list:
			recessive_gene_list.append(ensembl_gene)
			print('adding potential recessive gene '+gene_symbol)
			
		if gene_carrier_count_dict[ensembl_gene][sample] >= 2 and sample not in recessive_sample_list:
			recessive_sample_list.append(sample)
			print('adding potential recessive carrier '+sample)

for gene in recessive_gene_list:
	for sample in recessive_sample_list:
		if sample in gene_carrier_count_dict[gene] and gene_carrier_count_dict[gene][sample] >= 2:
			for line in gene_carrier_variant_dict[gene][sample]:
				if line not in variant_dict:
					opened_output_recessive.write(line+'\n')
					variant_dict[line] = 1
	



	