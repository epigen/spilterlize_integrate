        
# plot diagnostics of each processing step
rule plot_diagnostics:
    input:
        data = os.path.join(result_path,'{split}','{label}.csv'),
        annotation = os.path.join(result_path,'{split}','annotation.csv'),
    output:
        plot = report(os.path.join(result_path,'{split}','plots','{label}.png'),
                      caption="../report/diagnostic_plot.rst", 
                      category="{}_{}".format(config["project_name"], module_name),
                      subcategory="{split}",
                      labels={
                          "name": "{label}",
                          "type": "diagnostic plot"
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
