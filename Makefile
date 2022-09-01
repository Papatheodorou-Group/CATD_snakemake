prep : bin/getDag.sh bin/runPip.sh bin/basicSetup.sh bin/Snakefile bin/config.yaml
	@mkdir -p Input Input/Psuedobulks Input/Cell_splits && \
		cp bin/getDag.sh bin/runPip.sh bin/basicSetup.sh bin/Snakefile bin/config.yaml . && \
		bash basicSetup.sh && \
		rm basicSetup.sh

clean :
	@rm -rf Input Output Metrics Plots Benchmarks getDag.sh runPip.sh Snakefile config.yaml && \
		find . -maxdepth 1 -name '*scratch*' -exec rm -rf {} \; && \
		find . -maxdepth 1 -name '*benchmarks*' -exec rm -rf {} \; && \
		find . -maxdepth 1 -name '*pipOut*' -exec rm -rf {} \;
		find . -maxdepth 1 -name '*runSum*' -exec rm -rf {} \;
