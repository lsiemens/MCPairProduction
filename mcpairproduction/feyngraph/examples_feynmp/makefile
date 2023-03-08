# remove temporary files

files := minimal.pdf multiple.pdf
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
	rm *_mp* *.aux *.log *.pdf *.mp -f
