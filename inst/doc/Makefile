
## There is an old bug in texidvi that makes it not swallow the ~
## marker used to denote whitespace. This is actually due to fixing
## another bug they whereby you could not run texidvi on directory
## names containing a tilde (as we happen to do for Debian builds of R
## alpha/beta/rc releases). The 'tilde' bug will go away as it
## reportedly has been squashed upstream but I am still bitten by it
## on Ubuntu so for Dirk will insist on pdflatex and this variable helps
whoami=$(shell whoami)

all: clean genefu.pdf

pdfclean:
	rm -f genefu.pdf

clean:
	if test -d cache; then rm -fr cache; fi
	rm -f *.tex *.bbl *.blg *.aux *.out *.log *.tikz *.toc *.map *.sh *.dvi
	rm -f genefu-*.pdf
	cp -f vignette-src/genefu-fake.Rnw genefu.Rnw

genefu.pdf: vignette-src/genefu.Rnw
	cp -f vignette-src/genefu.Rnw .
	Rscript --vanilla -e "library(pgfSweave); pgfSweave('genefu.Rnw');"
	pdflatex genefu
	pdflatex genefu
	if test -d cache; then rm -fr cache; fi
	rm -f *.tex *.bbl *.blg *.aux *.out *.log *.tikz *.toc *.map *.sh *.dvi
	rm -f genefu-*.pdf
	cp -f vignette-src/genefu-fake.Rnw genefu.Rnw
