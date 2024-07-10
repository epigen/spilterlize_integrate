        
# plot diagnostics of each processing step
rule plot_diagnostics:
    input:
        data = os.path.join(result_path,'{split}','{label}.csv'),
        annotation = os.path.join(result_path,'{split}','annotation.csv'),
    output:
        diag_plot = report(os.path.join(result_path,'{split}','plots','{label}.png'),
                      caption="../report/diagnostic_plot.rst", 
                      category="{}_{}".format(config["project_name"], module_name),
                      subcategory="{split}",
                      labels={
                          "name": "{label}",
                          "type": "diagnostic plot"
                      }),
        cfa_plot = report(os.path.join(result_path,'{split}','plots','{label}_CFA.png'),
                      caption="../report/cfa_plot.rst", 
                      category="{}_{}".format(config["project_name"], module_name),
                      subcategory="{split}",
                      labels={
                          "name": "{label}",
                          "type": "confounding factor analysis plot"
                      }),
    params:
        partition = config.get("partition"),
    threads: config.get("threads", 1)
    resources:
        mem_mb=config.get("mem", "16000"),
    conda:
        "../envs/ggplot.yaml",
    log:
        "logs/rules/plot_diagnostics_{split}_{label}.log"
    script:
        "../scripts/plot_diagnostics.R"

# plot correlation heatmap of data after each processing step
rule plot_heatmap:
    input:
        data = os.path.join(result_path,'{split}','{label}.csv'),
        metadata = os.path.join(result_path,'{split}','annotation.csv'),
    output:
        heatmap_clustered = report(os.path.join(result_path,'{split}','plots','{label}_heatmap_clustered.png'),
                      caption="../report/correlation_heatmap.rst", 
                      category="{}_{}".format(config["project_name"], module_name),
                      subcategory="{split}",
                      labels={
                          "name": "{label}",
                          "type": "correlation heatmap clustered"
                      }),
        heatmap_sorted = report(os.path.join(result_path,'{split}','plots','{label}_heatmap_sorted.png'),
                      caption="../report/correlation_heatmap.rst", 
                      category="{}_{}".format(config["project_name"], module_name),
                      subcategory="{split}",
                      labels={
                          "name": "{label}",
                          "type": "correlation heatmap sorted"
                      }),
    params:
        partition = config.get("partition"),
    threads: config.get("threads", 1)
    resources:
        mem_mb=config.get("mem", "16000"),
    conda:
        "../envs/ComplexHeatmap.yaml",
    log:
        "logs/rules/plot_heatmap_{split}_{label}.log"
    script:
        "../scripts/plot_heatmap.R"
