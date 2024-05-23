"""
# Make logs folder 
# Make cluster folder inside logs folder
#Updated on 0712203 - PON is used in this case
#before running snakemake, do in tmux terminal:
ml snakemake/5.19.2-foss-2019b-Python-3.7.4
ml Java/11.0.2
ml picard/2.25.1-Java-11
ml GATK/4.1.8.1-GCCcore-8.3.0-Java-11
#ml GATK/4.2.6.1-GCCcore-11.2.0
ml tabix/0.2.6-GCCcore-8.3.0
ml annovar

#command to run snakemake (remove -np at end when done validating):
snakemake -s haplotypeCaller.snakefile --latency-wait 60 --restart-times 0 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output}" -j 100 -np
"""

configfile: "config/samples.yaml"
configfile: "config/config.yaml"


rule all:
    input:
        expand("results/{sample}.germline.vcf", sample=config["samples"]), 
        expand("results/{sample}.germline.hg38_multianno.vcf", sample=config["samples"])

rule runHaplotypeCaller:
    input:
        # pull in bam, bedfile (fields are two separate values under each sample in samples.yaml)
        bam = lambda wildcards: config["samples"][wildcards.sample]['normalBamPath'],
        reference = config["reference_genome"]
    output:
        vcf = "results/{sample}.germline.vcf"
    params:
        gatk = config["gatk"]
    log:
        "logs/cluster/haplotypeCaller_{sample}.log"
    shell:
        """
        {params.gatk} HaplotypeCaller -R {input.reference} -I {input.bam} -O {output.vcf} --tmp-dir /fh/fast/ha_g/projects/ProstateTAN/analysis_mutational/haplotypeCaller/tmp/ 2> {log}
        """

rule runAnnovar_germline:
    input:
        input_vcfs = "results/{sample}.germline.vcf"
    output:
        output_annovar_vcf = "results/{sample}.germline.hg38_multianno.vcf"
    params:
        annovar_python_script = config["annovar_python_script"],
        interpreter = config["interpreter"]
    log:
        "logs/runAnnovar_mutect2/{sample}_runAnnovar_snvs_mutect2.txt"
    shell:
        "({params.interpreter} {params.annovar_python_script} --input_vcf_file_path {input.input_vcfs} ) 2> {log}"
