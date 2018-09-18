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
	pdfjam --outfile $@ $< 1-53

cover-letter.pdf : local_pca_paper.pdf
	pdfjam --outfile $@ $< 54

review-responses.pdf : local_pca_paper.pdf
	pdfjam --outfile $@ $< 55-

# these make the pdf render too slow
local_pca_paper.pdf : review-responses.tex FigS_Together_MDS_plot_allchr.png FigS_pca_plot_allchr_human.png Fig1_allchr_Together_MDS_plot_compact_with_ChrX_inv.png MDS_allchr_Together_plot_samples_noinv.png Medicago_MDS_plot_allchr_win103_PCAk2.png

figs :
	echo "This doesn't make Figure 6 correctly yet (FIXME)."
	for x in $$(cat local_pca_paper.tex | grep includegr | tr '{}' ',,' | cut -f 2 -d ','); do make $$x.eps; done
	rm -f eps_figs/*
	readarray FIGS < <(head -n $$(grep -n '\\appendix' local_pca_paper.tex|cut -f 1 -d:) local_pca_paper.tex | grep includegr | tr '{}' ',,' | cut -f 2 -d ','| sed -e 's/ //g' ); for index in "$${!FIGS[@]}"; do x=$$(echo $${FIGS[index]} | tr -d ' '); ln -s $$x.eps Fig_$$((1+index)).eps; done
	readarray FIGS < <(tail -n +$$(grep -n '\\appendix' local_pca_paper.tex|cut -f 1 -d:) local_pca_paper.tex | grep includegr | tr '{}' ',,' | cut -f 2 -d ','| sed -e 's/ //g' ); for index in "$${!FIGS[@]}"; do x=$$(echo $${FIGS[index]} | tr -d ' '); ln -s $$x.eps Fig_S$$((1+index)).eps; done
	mv *.eps eps_figs

%.pdf : %.tex %.bbl
	while ( $(PDFLATEX) $<;  grep -q "Rerun to get" $*.log ) do true ; done

%.bbl : %.tex
	-$(PDFLATEX) $<
	bibtex $*.aux

clean :
	-rm *.aux *.bbl *.blg *.fff *.lof *.log *.out

