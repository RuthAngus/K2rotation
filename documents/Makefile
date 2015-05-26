LATEX       = pdflatex
BASH        = bash -c
ECHO        = echo
RM          = rm -rf
TMP_SUFFS   = pdf aux bbl blg log dvi ps eps out
CHECK_RERUN = grep Rerun $*.log

NAME = K2rotation

all: ${NAME}.pdf

%.pdf: %.tex
	${LATEX} $<
	bibtex ${NAME}
	${LATEX} $<
	${LATEX} $<
	open ${NAME}.pdf

clean:
	${RM} $(foreach suff, ${TMP_SUFFS}, ${NAME}.${suff})

publish:
	tar -czf publish.tar.gz K2rotation.tex K2rotation.bib \
		rawvbg_201183188.pdf 201545182.pdf 201211472.pdf \
		201433687.pdf 201444854.pdf 201607835.pdf \
		rotation201133037.pdf K2_rotation_201133037.pdf \
		rotation201142023.pdf K2_rotation_201142023.pdf sip_hist.pdf \
		vbg_hist.pdf RR_201339783.pdf EB_201473612.pdf \
		planet_201637175.pdf K2rotation.bbl
