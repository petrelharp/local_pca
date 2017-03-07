

%/run-summary.html : %/mds_coords.csv %/config.json
	Rscript -e 'templater::render_template("summarize_run.Rmd",change.rootdir=TRUE,output="$*/run-summary.html")'
