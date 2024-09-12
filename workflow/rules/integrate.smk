
# integrate normalized data using reComBat
rule integrate_recombat:
    input:
        normalized_data = os.path.join(result_path,'{split}','norm{norm_method}.csv'),
        annotation = os.path.join(result_path,'{split}','annotation.csv'),
    output:
        integrated_data = os.path.join(result_path,'{split}','norm{norm_method}_reComBat.csv'),
    threads: config.get("threads", 1)
    resources:
        mem_mb=config.get("mem", "16000"),
    conda:
        "../envs/reComBat.yaml",
    log:
        "logs/rules/integrate_recombat_{split}_norm{norm_method}.log"
    script:
        "../scripts/integrate_reComBat.py"