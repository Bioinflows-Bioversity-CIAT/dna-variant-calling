workflow GenotypeGVCFs_WF {
	String genome_ref

	String listChromosomeFile
	String listChromosomeLengthFile
	Array[String] listChromosome = read_lines(listChromosomeFile)
	String? listScaffoldFile
	String? listScaffoldLengthFile
	Boolean withScaffolds
	String basename

	String workDir
	String gDB_base
	Array[String] gDB_chr
	String? gDB_Sca
	String tmp
	String log

	String threads
	String gatk
	String tabix

	call getIntervals as intervals_chr {
		input:
			contigs=listChromosomeFile,
			lengths=listChromosomeLengthFile,
			basename=basename,
			log=log
	}

	scatter (interval in intervals_chr.out) {
		call GenotypeGVCFs as GTChromosome {
			input: 
				ref = genome_ref,
				outDir = workDir,
				gDB_base = gDB_base,
				gDb_suffix=interval[0],
				basename=basename,
				log = log,
				tmp = tmp,
				interval_contig = interval[0],
				interval_start = interval[1],
				interval_end = interval[2],
				threads = 2,
				gatk = gatk,
				tabix = tabix
		}
	}

	if(withScaffolds){
		call getIntervals as intervals_sca {
			input:
				contigs=listScaffoldFile,
				lengths=listScaffoldLengthFile,
				basename=basename,
				log=log
		}

		scatter (interval in intervals_sca.out) {
			call GenotypeGVCFs as GTScaffold {
				input: 
					ref = genome_ref,
					outDir = workDir,
					gDB_base=gDB_base,
					gDb_suffix = "Scaffolds",
					basename=basename,
					log = log,
					tmp = tmp,
					interval_contig = interval[0],
					interval_start = interval[1],
					interval_end = interval[2],
					threads = 1,
					gatk = gatk,
					tabix = tabix
			}
		}
	}

	call GatherVCFs as GatherVCFs_chr{
		input:
			ref = genome_ref,
			outDir = workDir,
			basename=basename,
			log = log,
			tmp = tmp,
			vcfs = GTChromosome.vcf,
			vcfsIndex = GTChromosome.vcfIndex,
			suffix = "chromosomes",
			threads = 2,
			gatk = gatk,
			tabix = tabix
	}

	if(withScaffolds){
		call GatherVCFs as GatherVCFs_sca{
			input:
				ref = genome_ref,
				outDir = workDir,
				basename=basename,
				log = log,
				tmp = tmp,
				vcfs = GTScaffold.vcf,
				vcfsIndex = GTScaffold.vcfIndex,
				suffix = "scaffolds",
				threads = 2,
				gatk = gatk,
				tabix = tabix
		}
		call GatherVCFs_final as GatherVCFs_final_sca{
			input:
				ref = genome_ref,
				outDir = workDir,
				basename=basename,
				log = log,
				tmp = tmp,
				chr = GatherVCFs_chr.vcf,
				sca = GatherVCFs_sca.vcf,
				threads = 2,
				gatk = gatk,
				tabix = tabix
		}
	}
}

#################################################################
task getIntervals {
	String contigs
	String lengths
	String basename
	String log

	command <<<
		paste ${contigs} ${lengths} | while IFS=$'\t' read -r contig length ; do
			start=1 ;
			end=1000000;
			while (( end<length )); do
				echo -e "$contig\t$start\t$end" ;
				(( start=start+1000000 )) ;
				(( end=end+1000000 )) ;
			done;
			echo -e "$contig\t$start\t$length" ;
		done ;
		exit 0 ;
	>>>
	runtime { 
		cpu: 2
		memory: "512 MB"
		my_job_name: "VD_GetIntervals_${basename}"
		my_log: "${log}/VD_GetIntervals_${basename}"
	}
	
	output {
		Array[Array[String]] out = read_tsv(stdout())
	}
}

task GenotypeGVCFs {
	String ref
	String outDir
	String gDB_base
	String gDb_suffix
	String basename
	String log
	String tmp
	String interval_contig
	String interval_start
	String interval_end
	String gatk
	String tabix

	Int threads
	Int mem=10000
	Int mems=8000

	command <<<
		mkdir -p ${outDir} ${tmp} ;
		VCF="${outDir}/${basename}_${interval_contig}_${interval_start}-${interval_end}.raw.vcf.gz" ;
		VCFINDEX="${outDir}/${basename}_${interval_contig}_${interval_start}-${interval_end}.raw.vcf.gz.tbi" ;
		if [[ ! -f $VCF || ! -f $VCFINDEX ]]; then
			${gatk} --java-options "-Xms${mems}M -Xmx${mem}M" GenotypeGVCFs \
				-R ${ref} \
				-L ${interval_contig}:${interval_start}-${interval_end} \
				-V gendb://${gDB_base}/${gDb_suffix} \
				-O ${outDir}/${basename}_${interval_contig}_${interval_start}-${interval_end}.raw.vcf.gz \
				--seconds-between-progress-updates 60 \
				--max-alternate-alleles 4 \
				--only-output-calls-starting-in-intervals \
				--tmp-dir ${tmp} ;
		else
			echo -e "${outDir}/${basename}_${interval_contig}_${interval_start}-${interval_end}.raw.vcf.gz exists" ;
		fi ;
		exit 0 ;
	>>>
	runtime { 
		cpu: threads
		memory: 2.0*mem/threads + " MB"
		my_job_name: "VD_GenotypeGVCFs_${basename}_${interval_contig}_${interval_start}-${interval_end}"
		my_log: "${log}/VD_GenotypeGVCFs_${basename}_${interval_contig}_${interval_start}-${interval_end}"
		maxRetries: 1
	}

	output {
		File vcf = "${outDir}/${basename}_${interval_contig}_${interval_start}-${interval_end}.raw.vcf.gz"
		File vcfIndex = "${outDir}/${basename}_${interval_contig}_${interval_start}-${interval_end}.raw.vcf.gz.tbi"
	}
}

task GatherVCFs {
	String ref
	String outDir
	String basename
	String log
	String tmp
	Array[String] vcfs
	Array[String] vcfsIndex
	String suffix
	String gatk
	String tabix

	Int threads
	Int mem=64000
	Int mems=60000

	command <<<
		mkdir -p ${outDir} ${tmp} ;
		VCF="${outDir}/${basename}.${suffix}.raw.vcf.gz" ;
		VCFINDEX="${outDir}/${basename}.${suffix}.raw.vcf.gz.tbi" ;
		if [[ ! -f $VCF || ! -f $VCFINDEX ]]; then
			${gatk} --java-options "-Xms${mems}M -Xmx${mem}M" GatherVcfs \
				-R ${ref} \
				--TMP_DIR ${tmp} \
				--VALIDATION_STRINGENCY LENIENT \
				--INPUT ${sep=' --INPUT ' vcfs} \
				--OUTPUT ${outDir}/${basename}.${suffix}.raw.vcf.gz ;
		fi ;
		${tabix} -p vcf ${outDir}/${basename}.${suffix}.raw.vcf.gz ;
		exit 0 ;
	>>>
	runtime { 
		cpu: threads
		memory: 2.0*mem/threads + " MB"
		my_job_name: "VD_GatherVCFs_${basename}_${suffix}"
		my_log: "${log}/VD_GatherVCFs_${basename}_${suffix}"
		maxRetries: 1
	}

	output {
		File vcf = "${outDir}/${basename}.${suffix}.raw.vcf.gz"
		File vcfIndex = "${outDir}/${basename}.${suffix}.raw.vcf.gz.tbi"
	}
}

task GatherVCFs_final {
	String ref
	String outDir
	String basename
	String log
	String tmp
	String chr
	String sca
	String gatk
	String tabix

	Int threads
	Int mem=32000
	Int mems=30000

	command <<<
		mkdir -p ${outDir} ${tmp} ;
		VCF="${outDir}/${basename}.raw.vcf.gz" ;
		VCFINDEX="${outDir}/${basename}.raw.vcf.gz.tbi" ;
		if [[ ! -f $VCF || ! -f $VCFINDEX ]]; then
			${gatk} --java-options "-Xms${mems}M -Xmx${mem}M" GatherVcfs \
				-R ${ref} \
				--TMP_DIR ${tmp} \
				--VALIDATION_STRINGENCY LENIENT \
				--INPUT ${chr} --INPUT ${sca} \
				--OUTPUT ${outDir}/${basename}.raw.vcf.gz ;
		fi ;
		${tabix} -p vcf ${outDir}/${basename}.raw.vcf.gz ;
		exit 0 ;
	>>>
	runtime { 
		cpu: threads
		memory: 2.0*mem/threads + " MB"
		my_job_name: "VD_GatherVCFs_${basename}"
		my_log: "${log}/VD_GatherVCFs_${basename}"
		maxRetries: 1
	}

	output {
		File vcf = "${outDir}/${basename}.raw.vcf.gz"
		File vcfIndex = "${outDir}/${basename}.raw.vcf.gz.tbi"
	}
}
