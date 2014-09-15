main: protocol.pdf clean
protocol.pdf: protocol.tex
	latex protocol.tex
	dvips protocol.dvi
	ps2pdf protocol.ps
	open protocol.pdf
clean:
	rm -f protocol.aux protocol.dvi protocol.log protocol.ps
flush: clean
	rm -f protocol.pdf


