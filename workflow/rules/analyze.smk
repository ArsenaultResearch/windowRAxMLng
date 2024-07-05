rule compute_metrics:
    input:
        expand("results/{window}.raxml.bestTree",window=WINDOWS),
        window_file = config["windows"],
        sample_file = config["samples"],
    output:
        out = "results/metrics.tsv",
    log:
        "logs/R_metrics.log",
    conda:
        "../envs/Ranalyze.yaml"
    resources:
        mem_mb=20000,
    shell:
        """
        Rscript workflow/scripts/GenerateMetrics.R {input.window_file} {input.sample_file}
        """
