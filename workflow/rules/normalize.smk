
# normalize by configured methods using edgeR::CalcNormFactors
rule norm_edgeR:
    input:
        filtered_counts = os.path.join(result_path,'{split}','filtered.csv'),
    output:
        normalized_counts = expand(os.path.join(result_path,'{{split}}','{method}.csv'), method=norm_edgeR_methods),
    params:
        partition = config.get("partition"),
        result_path = result_path,
        split = lambda w: "{}".format(w.split),
    threads: config.get("threads", 1)
    resources:
        mem_mb=config.get("mem", "16000"),
    conda:
        "../envs/edgeR.yaml",
    log:
        "logs/rules/norm_edgeR_{split}.log"
    script:
        "../scripts/norm_edgeR.R"
        
        
# normalize using CQN
rule norm_cqn:
    input:
        filtered_counts = os.path.join(result_path,'{split}','filtered.csv'),
        annotation = os.path.join(result_path,'{split}','annotation.csv'),
    output:
        normalized_counts = os.path.join(result_path,'{split}','normCQN.csv'),
        cqn_plot = report(os.path.join(result_path,'{split}','plots','normCQN_QRfit.png'),
                          caption="../report/CQN_plot.rst", 
                          category="{}_{}".format(config["project_name"], module_name),
                          subcategory="{split}",
                          labels={
                              "name": "normCQN",
                              "type": "QR fit plot"
                          }),
    params:
        partition = config.get("partition"),
        split = lambda w: "{}".format(w.split),
    threads: config.get("threads", 1)
    resources:
        mem_mb=config.get("mem", "16000"),
    conda:
        "../envs/cqn.yaml",
    log:
        "logs/rules/norm_cqn_{split}.log"
    script:
        "../scripts/norm_cqn.R"

        
# normalize using limma::voom
rule norm_voom:
    input:
        filtered_counts = os.path.join(result_path,'{split}','filtered.csv'),
    output:
        normalized_counts = os.path.join(result_path,'{split}','normVOOM.csv'),
        voom_plot = os.path.join(result_path,'{split}','normVOOM_mean_variance_trend.png'),
    params:
        partition = config.get("partition"),
        split = lambda w: "{}".format(w.split),
    threads: config.get("threads", 1)
    resources:
        mem_mb=config.get("mem", "16000"),
    conda:
        "../envs/edgeR.yaml",
    log:
        "logs/rules/norm_voom_{split}.log"
    script:
        "../scripts/norm_voom.R"

