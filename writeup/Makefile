.PHONY : figs

SHELL := /bin/bash

DISPLAYDIR = .

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

figs :
	echo "This doesn't make Figure 6 correctly yet (FIXME)."
	for x in $$(cat local_pca_paper.tex | grep includegr | tr '{}' ',,' | cut -f 2 -d ','); do make $$x.eps; done
	rm -f eps_figs/*
	readarray FIGS < <(head -n $$(grep -n '\\appendix' local_pca_paper.tex|cut -f 1 -d:) local_pca_paper.tex | grep includegr | tr '{}' ',,' | cut -f 2 -d ','| sed -e 's/ //g' ); for index in "$${!FIGS[@]}"; do x=$$(echo $${FIGS[index]} | tr -d ' '); ln -s $$x.eps Fig_$$((1+index)).eps; done
	readarray FIGS < <(tail -n +$$(grep -n '\\appendix' local_pca_paper.tex|cut -f 1 -d:) local_pca_paper.tex | grep includegr | tr '{}' ',,' | cut -f 2 -d ','| sed -e 's/ //g' ); for index in "$${!FIGS[@]}"; do x=$$(echo $${FIGS[index]} | tr -d ' '); ln -s $$x.eps Fig_S$$((1+index)).eps; done
	mv *.eps eps_figs


