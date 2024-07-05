rule subset_window:
    input:
        vcf = config["snps"],
    output:
        "results/{window}.vcf",
    params:
        extra = get_coords,
    wrapper:
        "v3.13.3/bio/vcftools/filter"

rule vcf2fasta:
    input:
        vcf = "results/{window}.vcf",
    output:
        fasta = "results/{window}.min4.fasta",
    log:
        "logs/{window}.vcf2fasta.log",
    resources:
        mem_mb=10000,
    shell:
        """
        workflow/scripts/vcf2phylip.py -i {input.vcf} -f
        """


rule make_phylogeny:
    input:
        fasta = "results/{window}.min4.fasta",
    output:
        phylo = "results/{window}.raxml.bestTree",
    log:
        "logs/{window}.raxmlng.log",
    threads: config["threads"],
    resources:
        mem_mb=20000,
    conda:
        "../envs/raxml_ng.yaml"
    shell:
        """
        raxml-ng -msa {input.fasta} \
        --model GTGTR4+G \
        --threads {threads} \
        --seed 2 \
        --search \
        --prefix "results/{wildcards.window}"
        """
