# these examples rely on the R/qtl package, www.rqtl.org
# install R/qtl package, within R, by typing install.packages("qtl")

all: Rmarkdown_example.html Exercises.html 

R_OPTS=--no-save --no-restore --no-init-file --no-site-file # vanilla, but with --environ

Rmarkdown_example.html:
	R ${R_OPTS} -e "rmarkdown::render('Rmarkdown_example.Rmd')"

Exercises.html:
	R ${R_OPTS} -e "rmarkdown::render('Exercises.Rmd')"

clean:
	rm -rf Rmarkdown_example.html Exercises.html 

