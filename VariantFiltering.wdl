workflow VariantFiltering_WF {
	String genome_ref
	String genome_ref_index
	String genome_ref_dict
	String? repeats

	String inputVcfFile
	String inputVcfFileIndex
	String basename

	Boolean? onlyBiallelic = true
	Boolean filterMissingMAF = false
	Float? filterMissingValue = 0.90
	Float? filterMAFValue = 0.05

	String? minDPValue = 2
	String? maxDPValue = 0
	String? GQValue = 0

	String? QD = 2.0	# <
	String? MQ = 40.0 	# <
	String? FS = 60.0	# >
	String? SOR = 3.0	# >
	String? MQR = -12.5	# <
	String? RPR = -8.0	# <

	String filterDensityScript
	String calculateStatsScript
	String calculateIStatsScript
	String calculateHistogramsScript

	String workDir
	String log = "${workDir}/log"
	String tmp = "${workDir}/tmp"
	
	String threads
	String java
	String picard
	String gatk
	String tabix
	String vcftools
	String bcftools
	String groovy

	String vcf_out = "${workDir}/vcf"

	call getSNPs {
			input:
				ref = genome_ref,
				ref_index = genome_ref_index,
				ref_dict = genome_ref_dict,
				inputVcfFile = inputVcfFile,
				inputVcfFileIndex = inputVcfFileIndex,
				basename = basename,
				onlyBiallelic = onlyBiallelic,
				outDir = vcf_out,
				log = log,
				tmp = tmp,
				threads = 2,
				bcftools = bcftools,
				vcftools = vcftools,
				tabix = tabix,
				calculateStatsScript = calculateStatsScript,
				calculateIStatsScript = calculateIStatsScript,
				calculateHistogramsScript = calculateHistogramsScript
	}

	if(defined(repeats)) {
		call getRepeats {
			input:
				ref = genome_ref,
				ref_index = genome_ref_index,
				ref_dict = genome_ref_dict,
				repeats = repeats,
				inputVcfFile = inputVcfFile,
				inputVcfFileIndex = inputVcfFileIndex,
				basename = basename,
				outDir = vcf_out,
				log = log,
				tmp = tmp,
				threads = 2,
				gatk = gatk,
				vcftools = vcftools
		}

		call filterINFO {
			input:
				ref = genome_ref,
				ref_index = genome_ref_index,
				ref_dict = genome_ref_dict,
				inputVcfFile = getSNPs.vcf,
				inputVcfFileIndex = getSNPs.index,
				repMaskVcfFile = getRepeats.vcf,
				repMaskVcfFileIndex = getRepeats.index,
				basename = basename(getSNPs.vcf,".vcf.gz"),
				outDir = vcf_out,
				log = log,
				tmp = tmp,
				threads = 2,
				gatk = gatk,
				vcftools = vcftools,
				calculateStatsScript = calculateStatsScript,
				calculateIStatsScript = calculateIStatsScript,
				calculateHistogramsScript = calculateHistogramsScript,
				QD = QD,
				MQ = MQ,
				FS = FS,
				SOR = SOR,
				MQR = MQR,
				RPR = RPR
		}

		call filterGT {
			input:
				inputVcfFile = filterINFO.vcf,
				inputVcfFileIndex = filterINFO.index,
				basename = basename(filterINFO.vcf,".vcf.gz"),
				minDPValue = minDPValue,
				maxDPValue = maxDPValue,
				GQValue = GQValue,
				outDir = vcf_out,
				log = log,
				tmp = tmp,
				threads = 2,
				vcftools = vcftools,
				tabix = tabix,
				calculateStatsScript = calculateStatsScript,
				calculateIStatsScript = calculateIStatsScript,
				calculateHistogramsScript = calculateHistogramsScript
		}

		if(filterMissingMAF){
			call filterMissing {
				input:
					inputVcfFile = filterGT.vcf,
					inputVcfFileIndex = filterGT.index,
					basename = basename(filterGT.vcf,".vcf.gz"),
					filterMissingValue = filterMissingValue,
					filterDensityScript = filterDensityScript,
					outDir = vcf_out,
					log = log,
					tmp = tmp,
					threads = 2,
					vcftools = vcftools,
					groovy = groovy,
					tabix = tabix,
					calculateStatsScript = calculateStatsScript,
					calculateIStatsScript = calculateIStatsScript,
					calculateHistogramsScript = calculateHistogramsScript
			}

			call filterMAF {
				input:
					inputVcfFile = filterMissing.vcf,
					inputVcfFileIndex = filterMissing.index,
					basename = basename(filterMissing.vcf,".vcf.gz"),
					filterMAFValue = filterMAFValue,
					outDir = vcf_out,
					log = log,
					tmp = tmp,
					threads = 2,
					vcftools = vcftools,
					tabix = tabix,
					calculateStatsScript = calculateStatsScript,
					calculateIStatsScript = calculateIStatsScript,
					calculateHistogramsScript = calculateHistogramsScript
			}
		}
	}

	if(!defined(repeats)) {
		call filterINFO_noRM {
			input:
				ref = genome_ref,
				ref_index = genome_ref_index,
				ref_dict = genome_ref_dict,
				inputVcfFile = getSNPs.vcf,
				inputVcfFileIndex = getSNPs.index,
				basename = basename(getSNPs.vcf,".vcf.gz"),
				outDir = vcf_out,
				log = log,
				tmp = tmp,
				threads = 2,
				gatk = gatk,
				vcftools = vcftools,
				calculateStatsScript = calculateStatsScript,
				calculateIStatsScript = calculateIStatsScript,
				calculateHistogramsScript = calculateHistogramsScript,
				QD = QD,
				MQ = MQ,
				FS = FS,
				SOR = SOR,
				MQR = MQR,
				RPR = RPR
		}

		call filterGT as filterGT_noRM {
			input:
				inputVcfFile = filterINFO_noRM.vcf,
				inputVcfFileIndex = filterINFO_noRM.index,
				basename = basename(filterINFO_noRM.vcf,".vcf.gz"),
				minDPValue = minDPValue,
				maxDPValue = maxDPValue,
				GQValue = GQValue,
				outDir = vcf_out,
				log = log,
				tmp = tmp,
				threads = 2,
				vcftools = vcftools,
				tabix = tabix,
				calculateStatsScript = calculateStatsScript,
				calculateIStatsScript = calculateIStatsScript,
				calculateHistogramsScript = calculateHistogramsScript
		}

		if(filterMissingMAF){
			call filterMissing as filterMissing_noRM {
				input:
					inputVcfFile = filterGT_noRM.vcf,
					inputVcfFileIndex = filterGT_noRM.index,
					basename = basename(filterGT_noRM.vcf,".vcf.gz"),
					filterMissingValue = filterMissingValue,
					filterDensityScript = filterDensityScript,
					outDir = vcf_out,
					log = log,
					tmp = tmp,
					threads = 2,
					vcftools = vcftools,
					groovy = groovy,
					tabix = tabix,
					calculateStatsScript = calculateStatsScript,
					calculateIStatsScript = calculateIStatsScript,
					calculateHistogramsScript = calculateHistogramsScript
			}

			call filterMAF as filterMAF_noRM {
				input:
					inputVcfFile = filterMissing_noRM.vcf,
					inputVcfFileIndex = filterMissing_noRM.index,
					basename = basename(filterMissing_noRM.vcf,".vcf.gz"),
					filterMAFValue = filterMAFValue,
					outDir = vcf_out,
					log = log,
					tmp = tmp,
					threads = 2,
					vcftools = vcftools,
					tabix = tabix,
					calculateStatsScript = calculateStatsScript,
					calculateIStatsScript = calculateIStatsScript,
					calculateHistogramsScript = calculateHistogramsScript
			}
		}
	}

}

############################################################################
task getRepeats {
	File ref
	File ref_index
	File ref_dict
	File repeats
	File inputVcfFile
	File inputVcfFileIndex
	String basename
	String outDir
	String log
	String tmp
	String gatk
	String vcftools

	Int? threads=2
	Int mem=30000
	Int mems=25000

	command <<<
		mkdir -p ${log} ${tmp} ${outDir} ;
		FILE="${outDir}/${basename}.reps.vcf.gz" ;
		INDEX="${outDir}/${basename}.reps.vcf.gz.tbi" ;
		if [[ ! -f $FILE || ! -f $INDEX ]]; then
			${gatk} --java-options "-Xms${mems}M -Xmx${mem}M" SelectVariants \
				-R ${ref} \
				-L ${repeats} \
 				-V ${inputVcfFile} \
 				-O ${outDir}/${basename}.reps.vcf.gz ;
 			${vcftools} --gzvcf ${outDir}/${basename}.reps.vcf.gz > ${outDir}/${basename}.reps.vcf.log 2>&1 ;
		else
			echo -e "${outDir}/${basename}.reps.vcf.gz exists" ;
		fi ;
		exit 0 ;
	>>>
	runtime { 
		cpu: threads
		memory: 2.0*mem/threads + " MB"
		my_job_name: "VF_GetRepeats_${basename}"
		my_log: "${log}/VF_GetRepeats_${basename}"
	}

	output {
		File vcf = "${outDir}/${basename}.reps.vcf.gz"
		File index = "${outDir}/${basename}.reps.vcf.gz.tbi"
	}
}

task getSNPs {
	File ref
	File ref_index
	File ref_dict
	File inputVcfFile
	File inputVcfFileIndex
	String basename
	String outDir
	String log
	String tmp
	Boolean? onlyBiallelic=true
	String bcftools
	String vcftools
	String tabix
	String calculateStatsScript
	String calculateIStatsScript
	String calculateHistogramsScript

	Int? threads=4
	Int mem=30000
	Int mems=25000

	command <<<
		mkdir -p ${log} ${tmp} ${outDir} ;
		FILE="${outDir}/${basename}.snps.vcf.gz" ;
		INDEX="${outDir}/${basename}.snps.vcf.gz.tbi" ;
		if [[ ! -f $FILE || ! -f $INDEX ]]; then
			if ${onlyBiallelic} ; then
				ALLELES="-m2 -M2 --exclude-types indels,mnps,other,bnd" ;
			else
				ALLELES="--exclude-types indels,mnps,other,bnd" ;
			fi ;
			${bcftools} view \
				$ALLELES \
				-Ov ${inputVcfFile} \
				| sed 's/=nan\([;\t]\)/=NaN\1/g' \
				| bgzip -c > ${outDir}/${basename}.snps.vcf.gz ;
			${tabix} -p vcf ${outDir}/${basename}.snps.vcf.gz ;
			${vcftools} --gzvcf ${outDir}/${basename}.snps.vcf.gz > ${outDir}/${basename}.snps.vcf.log 2>&1 ;
		else
			echo -e "${outDir}/${basename}.snps.vcf.gz exists" ;
		fi ;

		#if [[ -f ${outDir}/${basename}.snps.vcf.gz && ! -f ${outDir}/${basename}.snps.vcf.stats ]]; then
		#	bash ${calculateStatsScript} ${outDir}/${basename}.snps.vcf.gz ${outDir}/${basename}.snps.vcf.stats ;
		#fi ;
		#if [[ -f ${outDir}/${basename}.snps.vcf.stats && ! -f ${outDir}/${basename}.snps.vcf.Ho_He_PIC_MAF_Mis.hist.pdf ]]; then
		#	bash ${calculateHistogramsScript} ${outDir}/${basename}.snps.vcf.stats ${outDir}/${basename}.snps.vcf ;
		#fi ;
		#if [[ -f ${outDir}/${basename}.snps.vcf.gz && ! -f ${outDir}/${basename}.snps.vcf.istats ]]; then
		#	bash ${calculateIStatsScript} ${outDir}/${basename}.snps.vcf.gz ${outDir}/${basename}.snps.vcf.istats 32 ;
		#fi ;

		exit 0 ;
	>>>
	runtime { 
		cpu: threads
		memory: 2.0*mem/threads + " MB"
		my_job_name: "VF_GetSNPs_${basename}"
		my_log: "${log}/VF_GetSNPs_${basename}"
	}

	output {
		File vcf = "${outDir}/${basename}.snps.vcf.gz"
		File index = "${outDir}/${basename}.snps.vcf.gz.tbi"
	}
}

task filterINFO {
	File ref
	File ref_index
	File ref_dict
	File inputVcfFile
	File inputVcfFileIndex
	File repMaskVcfFile
	File repMaskVcfFileIndex
	String basename
	String outDir
	String log
	String tmp
	String gatk
	String vcftools
	String calculateStatsScript
	String calculateIStatsScript
	String calculateHistogramsScript

	String? QD=2.0
	String? MQ=40.0
	String? FS=60.0
	String? SOR=3.0
	String? MQR=-12.5
	String? RPR=-8.0

	Int? threads=2
	Int mem=30000
	Int mems=25000

	command <<<
		mkdir -p ${log} ${tmp} ${outDir} ;
		FILE="${outDir}/${basename}.filter_info.vcf.gz" ;
		INDEX="${outDir}/${basename}.filter_info.vcf.gz.tbi" ;
		if [[ ! -f $FILE || ! -f $INDEX ]]; then
			${gatk} --java-options "-Xms${mems}M -Xmx${mem}M" VariantFiltration \
				--verbosity ERROR \
				-R ${ref} \
 				-V ${inputVcfFile} \
 				--filter-name "LowQD" --filter-expression "QD<${QD}" \
 				--filter-name "LowMQ" --filter-expression "MQ<${MQ}" \
 				--filter-name "HiFS" --filter-expression "FS>${FS}" \
 				--filter-name "HiSOR" --filter-expression "SOR>${SOR}" \
 				--filter-name "LowMQRankSum" --filter-expression "MQRankSum<${MQR}" \
 				--filter-name "LowReadPosRankSum" --filter-expression "ReadPosRankSum<${RPR}" \
 				--mask-name "RepMask" --mask ${repMaskVcfFile} \
 				-O ${outDir}/${basename}.filter_info_marked.vcf.gz ;
 			${gatk} --java-options "-Xms${mems}M -Xmx${mem}M" SelectVariants \
				-R ${ref} \
				--exclude-filtered \
 				-V ${outDir}/${basename}.filter_info_marked.vcf.gz \
 				-O ${outDir}/${basename}.filter_info.vcf.gz ;
 			${vcftools} --gzvcf ${outDir}/${basename}.filter_info.vcf.gz > ${outDir}/${basename}.filter_info.vcf.log 2>&1 ;
 			rm ${outDir}/${basename}.filter_info_marked.vcf.gz ;
 			rm ${outDir}/${basename}.filter_info_marked.vcf.gz.tbi ;
		else
			echo -e "${outDir}/${basename}.filter_info.vcf.gz exists" ;
		fi ;

		#if [[ -f ${outDir}/${basename}.filter_info.vcf.gz && ! -f ${outDir}/${basename}.filter_info.vcf.stats ]]; then
		#	bash ${calculateStatsScript} ${outDir}/${basename}.filter_info.vcf.gz ${outDir}/${basename}.filter_info.vcf.stats ;
		#fi ;
		#if [[ -f ${outDir}/${basename}.filter_info.vcf.stats && ! -f ${outDir}/${basename}.filter_info.vcf.Ho_He_PIC_MAF_Mis.hist.pdf ]]; then
		#	bash ${calculateHistogramsScript} ${outDir}/${basename}.filter_info.vcf.stats ${outDir}/${basename}.filter_info.vcf ;
		#fi ;
		#if [[ -f ${outDir}/${basename}.filter_info.vcf.gz && ! -f ${outDir}/${basename}.filter_info.vcf.istats ]]; then
		#	bash ${calculateIStatsScript} ${outDir}/${basename}.filter_info.vcf.gz ${outDir}/${basename}.filter_info.vcf.istats 32 ;
		#fi ;

		exit 0 ;
	>>>
	runtime { 
		cpu: threads
		memory: 2.0*mem/threads + " MB"
		my_job_name: "VF_FilterINFO_${basename}"
		my_log: "${log}/VF_FilterINFO_${basename}"
	}

	output {
		File vcf = "${outDir}/${basename}.filter_info.vcf.gz"
		File index = "${outDir}/${basename}.filter_info.vcf.gz.tbi"
	}
}

task filterINFO_noRM {
	File ref
	File ref_index
	File ref_dict
	File inputVcfFile
	File inputVcfFileIndex
	String basename
	String outDir
	String log
	String tmp
	String gatk
	String vcftools
	String calculateStatsScript
	String calculateIStatsScript
	String calculateHistogramsScript

	String? QD=2.0
	String? MQ=40.0
	String? FS=60.0
	String? SOR=3.0
	String? MQR=-12.5
	String? RPR=-8.0

	Int? threads=2
	Int mem=30000
	Int mems=25000

	command <<<
		mkdir -p ${log} ${tmp} ${outDir} ;
		FILE="${outDir}/${basename}.filter_info.vcf.gz" ;
		INDEX="${outDir}/${basename}.filter_info.vcf.gz.tbi" ;
		if [[ ! -f $FILE || ! -f $INDEX ]]; then
			${gatk} --java-options "-Xms${mems}M -Xmx${mem}M" VariantFiltration \
				--verbosity ERROR \
				-R ${ref} \
 				-V ${inputVcfFile} \
 				--filter-name "LowQD" --filter-expression "QD<${QD}" \
 				--filter-name "LowMQ" --filter-expression "MQ<${MQ}" \
 				--filter-name "HiFS" --filter-expression "FS>${FS}" \
 				--filter-name "HiSOR" --filter-expression "SOR>${SOR}" \
 				--filter-name "LowMQRankSum" --filter-expression "MQRankSum<${MQR}" \
 				--filter-name "LowReadPosRankSum" --filter-expression "ReadPosRankSum<${RPR}" \
 				-O ${outDir}/${basename}.filter_info_marked.vcf.gz ;
 			${gatk} --java-options "-Xms${mems}M -Xmx${mem}M" SelectVariants \
				-R ${ref} \
				--exclude-filtered \
 				-V ${outDir}/${basename}.filter_info_marked.vcf.gz \
 				-O ${outDir}/${basename}.filter_info.vcf.gz ;
 			${vcftools} --gzvcf ${outDir}/${basename}.filter_info.vcf.gz > ${outDir}/${basename}.filter_info.vcf.log 2>&1 ;
 			rm ${outDir}/${basename}.filter_info_marked.vcf.gz ;
 			rm ${outDir}/${basename}.filter_info_marked.vcf.gz.tbi ;
		else
			echo -e "${outDir}/${basename}.filter_info.vcf.gz exists" ;
		fi ;

		#if [[ -f ${outDir}/${basename}.filter_info.vcf.gz && ! -f ${outDir}/${basename}.filter_info.vcf.stats ]]; then
		#	bash ${calculateStatsScript} ${outDir}/${basename}.filter_info.vcf.gz ${outDir}/${basename}.filter_info.vcf.stats ;
		#fi ;
		#if [[ -f ${outDir}/${basename}.filter_info.vcf.stats && ! -f ${outDir}/${basename}.filter_info.vcf.Ho_He_PIC_MAF_Mis.hist.pdf ]]; then
		#	bash ${calculateHistogramsScript} ${outDir}/${basename}.filter_info.vcf.stats ${outDir}/${basename}.filter_info.vcf ;
		#fi ;
		#if [[ -f ${outDir}/${basename}.filter_info.vcf.gz && ! -f ${outDir}/${basename}.filter_info.vcf.istats ]]; then
		#	bash ${calculateIStatsScript} ${outDir}/${basename}.filter_info.vcf.gz ${outDir}/${basename}.filter_info.vcf.istats 32 ;
		#fi ;

		exit 0 ;
	>>>
	runtime { 
		cpu: threads
		memory: 2.0*mem/threads + " MB"
		my_job_name: "VF_FilterINFO_${basename}"
		my_log: "${log}/VF_FilterINFO_${basename}"
	}

	output {
		File vcf = "${outDir}/${basename}.filter_info.vcf.gz"
		File index = "${outDir}/${basename}.filter_info.vcf.gz.tbi"
	}
}

task filterGT{
	File inputVcfFile
	File inputVcfFileIndex
	String basename
	String minDPValue
	String maxDPValue
	String GQValue
	String outDir
	String log
	String tmp
	String vcftools
	String tabix
	String calculateStatsScript
	String calculateIStatsScript
	String calculateHistogramsScript

	Int? threads=2
	Int mem=30000
	Int mems=25000

	command <<<
		mkdir -p ${log} ${tmp} ${outDir} ;
		if (( ${minDPValue} > 0 )) && (( ${maxDPValue} > 0 )) && (( ${GQValue} > 0 )) ; then
			FILE="${outDir}/${basename}_minDP${minDPValue}_maxDP${maxDPValue}_GQ${GQValue}.vcf" ;
			INDEX="$FILE.gz.tbi" ;
			if [[ ! -f $FILE.gz || ! -f $INDEX ]]; then
				${vcftools} --gzvcf ${inputVcfFile} --minDP ${minDPValue} --maxDP ${maxDPValue} --minGQ ${GQValue} --recode --recode-INFO-all --stdout | bgzip -c > $FILE.gz ;
				${vcftools} --gzvcf $FILE.gz > $FILE.log 2>&1 ;
 				${tabix} -p vcf $FILE.gz ;
 				echo -e "$FILE" ;
			else
				echo -e "$FILE.gz exists" ;
			fi ;
		elif (( ${minDPValue} > 0 )) && (( ${maxDPValue} > 0 )) ; then
			FILE="${outDir}/${basename}_minDP${minDPValue}_maxDP${maxDPValue}.vcf" ;
			INDEX="$FILE.gz.tbi" ;
			if [[ ! -f $FILE.gz || ! -f $INDEX ]]; then
				${vcftools} --gzvcf ${inputVcfFile} --minDP ${minDPValue} --maxDP ${maxDPValue} --recode --recode-INFO-all --stdout | bgzip -c > $FILE.gz ;
				${vcftools} --gzvcf $FILE.gz > $FILE.log 2>&1 ;
 				${tabix} -p vcf $FILE.gz ;
 				echo -e "$FILE" ;
			else
				echo -e "$FILE.gz exists" ;
			fi ;
		elif (( ${minDPValue} > 0 )) ; then
			FILE="${outDir}/${basename}_minDP${minDPValue}.vcf" ;
			INDEX="$FILE.gz.tbi" ;
			if [[ ! -f $FILE.gz || ! -f $INDEX ]]; then
				${vcftools} --gzvcf ${inputVcfFile} --minDP ${minDPValue} --recode --recode-INFO-all --stdout | bgzip -c > $FILE.gz ;
				${vcftools} --gzvcf $FILE.gz > $FILE.log 2>&1 ;
 				${tabix} -p vcf $FILE.gz ;
 				echo -e "$FILE" ;
			else
				echo -e "$FILE.gz exists" ;
			fi ;
		else
			FILE="${outDir}/${basename}.vcf" ;
			INDEX="$FILE.gz.tbi" ;
			if [[ ! -f $FILE.gz || ! -f $INDEX ]]; then
				${vcftools} --gzvcf ${inputVcfFile} --recode --recode-INFO-all --stdout | bgzip -c > $FILE.gz ;
				${vcftools} --gzvcf $FILE.gz > $FILE.log 2>&1 ;
 				${tabix} -p vcf $FILE.gz ;
 				echo -e "$FILE" ;
			else
				echo -e "$FILE.gz exists" ;
			fi ;
		fi ;

		exit 0 ;
	>>>
	runtime { 
		cpu: threads
		memory: 2.0*mem/threads + " MB"
		my_job_name: "VF_FilterGT_${basename}"
		my_log: "${log}/VF_FilterGT_${basename}"
	}

	output {
		File vcf = read_string(stdout()) + ".gz"
		File index = read_string(stdout()) + ".gz.tbi"
	}
}

task filterMissing {
	File inputVcfFile
	File inputVcfFileIndex
	String basename
	Float filterMissingValue
	String filterDensityScript
	String outDir
	String log
	String tmp
	String vcftools
	String groovy
	String tabix
	String calculateStatsScript
	String calculateIStatsScript
	String calculateHistogramsScript

	Int? threads=2
	Int mem=30000
	Int mems=25000

	command <<<
		mkdir -p ${log} ${tmp} ${outDir} ;
		FILE="${outDir}/${basename}_den${filterMissingValue}.vcf.gz" ;
		INDEX="${outDir}/${basename}_den${filterMissingValue}.vcf.gz.tbi" ;
		if [[ ! -f $FILE || ! -f $INDEX ]]; then
			export JAVA_OPTS="$JAVA_OPTS -Xms${mems}M -Xmx${mem}M" ;
			${groovy} ${filterDensityScript} ${inputVcfFile} ${outDir}/${basename}_den${filterMissingValue}.vcf ${filterMissingValue} ;
			bgzip ${outDir}/${basename}_den${filterMissingValue}.vcf ;

 			${vcftools} --gzvcf ${outDir}/${basename}_den${filterMissingValue}.vcf.gz > ${outDir}/${basename}_den${filterMissingValue}.vcf.log 2>&1 ;
 			${tabix} -p vcf ${outDir}/${basename}_den${filterMissingValue}.vcf.gz ;
		else
			echo -e "${outDir}/${basename}_den${filterMissingValue}.vcf.gz exists" ;
		fi ;

		#if [[ -f ${outDir}/${basename}_den${filterMissingValue}.vcf.gz && ! -f ${outDir}/${basename}_den${filterMissingValue}.vcf.stats ]]; then
		#	bash ${calculateStatsScript} ${outDir}/${basename}_den${filterMissingValue}.vcf.gz ${outDir}/${basename}_den${filterMissingValue}.vcf.stats ;
		#fi ;
		#if [[ -f ${outDir}/${basename}_den${filterMissingValue}.vcf.stats && ! -f ${outDir}/${basename}_den${filterMissingValue}.vcf.Ho_He_PIC_MAF_Mis.hist.pdf ]]; then
		#	bash ${calculateHistogramsScript} ${outDir}/${basename}_den${filterMissingValue}.vcf.stats ${outDir}/${basename}_den${filterMissingValue}.vcf ;
		#fi ;
		#if [[ -f ${outDir}/${basename}_den${filterMissingValue}.vcf.gz && ! -f ${outDir}/${basename}_den${filterMissingValue}.vcf.istats ]]; then
		#	bash ${calculateIStatsScript} ${outDir}/${basename}_den${filterMissingValue}.vcf.gz ${outDir}/${basename}_den${filterMissingValue}.vcf.istats 32 ;
		#fi ;

		exit 0 ;
	>>>
	runtime { 
		cpu: threads
		memory: 2.0*mem/threads + " MB"
		my_job_name: "VF_FilterMissing_${basename}"
		my_log: "${log}/VF_FilterMissing_${basename}"
	}

	output {
		File vcf = "${outDir}/${basename}_den${filterMissingValue}.vcf.gz"
		File index = "${outDir}/${basename}_den${filterMissingValue}.vcf.gz.tbi"
	}
}

task filterMAF {
	File inputVcfFile
	File inputVcfFileIndex
	String basename
	Float filterMAFValue
	String outDir
	String log
	String tmp
	String vcftools
	String tabix
	String calculateStatsScript
	String calculateIStatsScript
	String calculateHistogramsScript

	Int? threads=2
	Int mem=30000
	Int mems=25000

	command <<<
		mkdir -p ${log} ${tmp} ${outDir} ;
		FILE="${outDir}/${basename}_MAF${filterMAFValue}.vcf.gz" ;
		INDEX="${outDir}/${basename}_MAF${filterMAFValue}.vcf.gz.tbi" ;
		if [[ ! -f $FILE || ! -f $INDEX ]]; then
			${vcftools} --gzvcf ${inputVcfFile} --maf ${filterMAFValue} --recode --recode-INFO-all --stdout \
				| bgzip -c > ${outDir}/${basename}_MAF${filterMAFValue}.vcf.gz ;

 			${vcftools} --gzvcf ${outDir}/${basename}_MAF${filterMAFValue}.vcf.gz > ${outDir}/${basename}_MAF${filterMAFValue}.vcf.log 2>&1 ;
 			${tabix} -p vcf ${outDir}/${basename}_MAF${filterMAFValue}.vcf.gz ;
		else
			echo -e "${outDir}/${basename}_MAF${filterMAFValue}.vcf.gz exists" ;
		fi ;

		#if [[ -f ${outDir}/${basename}_MAF${filterMAFValue}.vcf.gz && ! -f ${outDir}/${basename}_MAF${filterMAFValue}.vcf.stats ]]; then
		#	bash ${calculateStatsScript} ${outDir}/${basename}_MAF${filterMAFValue}.vcf.gz ${outDir}/${basename}_MAF${filterMAFValue}.vcf.stats ;
		#fi ;
		#if [[ -f ${outDir}/${basename}_MAF${filterMAFValue}.vcf.stats && ! -f ${outDir}/${basename}_MAF${filterMAFValue}.vcf.Ho_He_PIC_MAF_Mis.hist.pdf ]]; then
		#	bash ${calculateHistogramsScript} ${outDir}/${basename}_MAF${filterMAFValue}.vcf.stats ${outDir}/${basename}_MAF${filterMAFValue}.vcf ;
		#fi ;
		#if [[ -f ${outDir}/${basename}_MAF${filterMAFValue}.vcf.gz && ! -f ${outDir}/${basename}_MAF${filterMAFValue}.vcf.istats ]]; then
		#	bash ${calculateIStatsScript} ${outDir}/${basename}_MAF${filterMAFValue}.vcf.gz ${outDir}/${basename}_MAF${filterMAFValue}.vcf.istats 32 ;
		#fi ;

		exit 0 ;
	>>>
	runtime { 
		cpu: threads
		memory: 2.0*mem/threads + " MB"
		my_job_name: "VF_FilterMAF_${basename}"
		my_log: "${log}/VF_FilterMAF_${basename}"
	}

	output {
		File vcf = "${outDir}/${basename}_MAF${filterMAFValue}.vcf.gz"
		File index = "${outDir}/${basename}_MAF${filterMAFValue}.vcf.gz.tbi"
	}
}
