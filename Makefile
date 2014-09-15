main: angelino.pdf clean
angelino.pdf: angelino.tex refs.bib files/*.tex latex_macros.tex # thesis.sty
	pdflatex angelino
	pdflatex angelino
	bibtex angelino
	pdflatex angelino
	pdflatex angelino
clean:
	rm -f angelino.aux angelino.bbl angelino.log angelino.blg angelino.out
flush: clean
	rm -f angelino.pdf


