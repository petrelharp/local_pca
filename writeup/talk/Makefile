.PHONY : clean

PDFLATEX=pdflatex -interaction nonstopmode

clean : 
	-rm *.{aux,log,bbl,blg,out}

%.html : %.Rmd
	Rscript -e "templater::render_template(\"$<\", output=\"$@\", change.rootdir=TRUE)"

%.pdf : %.tex
	while ( $(PDFLATEX) $<;  grep -q "Rerun to get" $*.log ) do true ; done

evolution-2017-talk.pdf : files/what_is_covariance.pdf ../the-method-diagram.pdf files/the-method-diagram-modified.png

%.svg : %.pdf
	inkscape $< --export-plain-svg=$@

%.png : %.pdf
	convert -density 300 $< -flatten $@

%.pdf : %.ink.svg
	inkscape $< --export-pdf=$@

