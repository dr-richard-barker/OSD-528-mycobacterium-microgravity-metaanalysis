# Makefile for Manuscript Compilation

all: pdf docx

pdf:
	pdflatex -interaction=nonstopmode main.tex
	bibtex main
	pdflatex -interaction=nonstopmode main.tex
	pdflatex -interaction=nonstopmode main.tex
	cp main.pdf OSD528_Microbial_Microgravity_Manuscript.pdf

docx:
	python3 build_manuscript.py

clean:
	rm -f *.aux *.bbl *.blg *.log *.out *.toc
