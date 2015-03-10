# Makefile for Rmarkdown files

all: Rmarkdown_example.html Exercises.html Manipulation.html nihexporter.html

R_OPTS=--no-save --no-restore --no-init-file --no-site-file # vanilla, but with --environ

nihexporter.html:
	R ${R_OPTS} -e "rmarkdown::render('nihexporter.Rmd')"

Rmarkdown_example.html:
	R ${R_OPTS} -e "rmarkdown::render('Rmarkdown_example.Rmd')"

Exercises.html:
	R ${R_OPTS} -e "rmarkdown::render('Exercises.Rmd')"

Manipulation.html:
	R ${R_OPTS} -e "rmarkdown::render('Manipulation.Rmd')"

clean:
	rm -rf Rmarkdown_example.html Exercises.html Manipulation.html nihexporter.html

