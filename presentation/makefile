# remove temporary files

files := presentation.pdf
VIEWPDF := xdg-open
LATEX := pdflatex

show: ${files}
	${foreach var, $?, ${VIEWPDF} ${var};}

all: ${files}

%.pdf : %.tex
	${LATEX} $?
	mpost ${basename $?}_mp.mp
	${LATEX} $?

clean:
	rm *_mp* *.mp -f
	rm *.aux *.log *.nav *.snm *.out *.toc *.pdf -f
