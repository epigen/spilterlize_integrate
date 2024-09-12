

# split data according to config and annotation
rule split:
    input:
        data = config["data"],
        annotation = config["annotation"],
    output:
        counts = expand(os.path.join(result_path,'{split}','counts.csv'), split=splits),
        annotations = expand(os.path.join(result_path,'{split}','annotation.csv'), split=splits),
    params:
        split_by = ["all"] + config["split_by"],
        result_path = result_path,
    threads: config.get("threads", 1)
    resources:
        mem_mb=config.get("mem", "16000"),
    conda:
        "../envs/python_basics.yaml",
    log:
        "logs/rules/split_data.log"
    script:
        "../scripts/split_data.py"
        
        
        
# filter features by "expression"
rule filter_features:
    input:
        data = os.path.join(result_path,'{split}','counts.csv'),
        annotation = os.path.join(result_path,'{split}','annotation.csv'),
    output:
        filtered_counts = os.path.join(result_path,'{split}','filtered.csv'),
    params:
        filter_parameters = config["filter_parameters"]
    threads: config.get("threads", 1)
    resources:
        mem_mb=config.get("mem", "16000"),
    conda:
        "../envs/edgeR.yaml",
    log:
        "logs/rules/filter_{split}.log"
    script:
        "../scripts/filter_features.R"
        
        
# select Highly Variable Features (HVF)
rule select_hvf:
    input:
        data = os.path.join(result_path,'{split}','{transformed_data}.csv'),
    output:
        hvf_data = os.path.join(result_path,'{split}','{transformed_data}_HVF.csv'),
        hvf_plot = report(os.path.join(result_path,'{split}','plots','{transformed_data}_HVF_selection.png'),
                          caption="../report/HVF_plot.rst",
                          category="{}_{}".format(config["project_name"], module_name),
                          subcategory="{split}",
                          labels={
                              "name": "{transformed_data}",
                              "type": "HVF selection"
                          }),
    threads: config.get("threads", 1)
    resources:
        mem_mb=config.get("mem", "16000"),
    conda:
        "../envs/python_basics.yaml",
    log:
        "logs/rules/hvf_{split}_{transformed_data}.log"
    script:
        "../scripts/select_hvf.py"