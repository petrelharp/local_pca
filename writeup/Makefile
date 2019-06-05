.PHONY : figs all

SHELL := /bin/bash

DISPLAYDIR = .

PDFLATEX := pdflatex

# save inkscape svg files as .ink.svg and this'll do the right thing
$(DISPLAYDIR)/%.svg : %.ink.svg
	inkscape $< --export-plain-svg=$@

$(DISPLAYDIR)/%.pdf : %.ink.svg
	inkscape $< --export-pdf=$@

$(DISPLAYDIR)/%.svg : %.pdf
	inkscape $< --export-plain-svg=$@

$(DISPLAYDIR)/%.png : %.pdf
	convert -density 300 $< -flatten $@

$(DISPLAYDIR)/%.eps : %.pdf
	inkscape $< --export-eps=$@

all : with-line-numbers.pdf cover-letter.pdf review-responses.pdf

with-line-numbers.pdf : local_pca_paper.pdf
	pdfjam --outfile $@ $< 1-54

cover-letter.pdf : local_pca_paper.pdf
	pdfjam --outfile $@ $< 55

review-responses.pdf : local_pca_paper.pdf
	pdfjam --outfile $@ $< 56-

# these make the pdf render too slow
local_pca_paper.pdf : review-responses.tex FigS_Together_MDS_plot_allchr.png FigS_pca_plot_allchr_human.png Fig1_allchr_Together_MDS_plot_compact_with_ChrX_inv.png MDS_allchr_Together_plot_samples_noinv.png Medicago_MDS_plot_allchr_win103_PCAk2.png

figs : supp_figs
	for x in $$(tail -n +$$( grep -n '\\begin{document' local_pca_paper.tex|cut -f 1 -d:) local_pca_paper.tex | grep includegr | grep -v png | tr '{}' ',,' | cut -f 2 -d ','); do make $$x.eps; done
	readarray FIGS < <(head -n $$(grep -n '\\appendix' local_pca_paper.tex|cut -f 1 -d:) local_pca_paper.tex | tail -n +$$( grep -n '\\begin{document' local_pca_paper.tex|cut -f 1 -d:) | grep includegr | tr '{}' ',,' | cut -f 2 -d ','| sed -e 's/ //g' | grep -v png); for index in "$${!FIGS[@]}"; do x=$$(echo $${FIGS[index]} | tr -d ' '); ln -s $$x.eps Fig_$$((1+index)).eps; done
	# readarray FIGS < <(tail -n +$$(grep -n '\\appendix' local_pca_paper.tex|cut -f 1 -d:) local_pca_paper.tex | grep includegr | tr '{}' ',,' | cut -f 2 -d ','| sed -e 's/ //g' | grep -v png); for index in "$${!FIGS[@]}"; do x=$$(echo $${FIGS[index]} | tr -d ' '); ln -s $$x.eps Fig_S$$((1+index)).eps; done
	mv *.eps publication_stuff
	mv publication_stuff/Fig_6.eps publication_stuff/Fig_6A.eps
	mv publication_stuff/Fig_7.eps publication_stuff/Fig_6B.eps
	mv publication_stuff/Fig_8.eps publication_stuff/Fig_7.eps

supp_figs : local_pca_paper.pdf
	pdfjam --outfile publication_stuff/Fig_S1.pdf local_pca_paper.pdf 35
	pdfjam --outfile publication_stuff/Fig_S2.pdf local_pca_paper.pdf 36
	pdfjam --outfile publication_stuff/Table_S1.pdf local_pca_paper.pdf 37
	for x in $$(seq 17); do pdfjam --outfile publication_stuff/Fig_S$$((2 + $$x)).pdf local_pca_paper.pdf $$((37 + $$x)); done

%.pdf : %.tex %.bbl
	while ( $(PDFLATEX) $<;  grep -q "Rerun to get" $*.log ) do true ; done

%.bbl : %.tex
	-$(PDFLATEX) $<
	bibtex $*.aux

clean :
	-rm *.aux *.bbl *.blg *.fff *.lof *.log *.out

