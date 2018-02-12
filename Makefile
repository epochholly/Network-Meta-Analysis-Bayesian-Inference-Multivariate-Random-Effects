objects := $(wildcard R/*.R) DESCRIPTION
version := $(shell grep "Version" DESCRIPTION | sed "s/Version: //")
pkg := $(shell grep "Package" DESCRIPTION | sed "s/Package: //")
tar := $(pkg)_$(version).tar.gz
checkLog := $(pkg).Rcheck/00check.log
yr := $(shell date +"%Y")
dt := $(shell date +"%Y-%m-%d")


.PHONY: check
check: $(checkLog)

.PHONY: build
build: $(tar)

.PHONY: install
install: $(tar)
	R CMD INSTALL $(tar)

$(tar): $(objects)
	@if [ "$$(uname)" == "Darwin" ];\
	then echo "remeber to update date and version number";\
	else make -s updateMeta;\
	fi;\
	Rscript -e "library(methods); roxygen2::roxygenize();";
	R CMD build .

$(checkLog): $(tar)
# R CMD check --as-cran $(tar)
	R CMD check $(tar)

## update copyright year in HEADER, R script and date in DESCRIPTION
.PHONY: updateMeta
updateMeta:
	@sed -i "s/Date: [0-9]\{4\}-[0-9]\{1,2\}-[0-9]\{1,2\}/Date: $(dt)/" DESCRIPTION
# @sed -i "s/version [0-9]\.[0-9]\.[0-9]\(\.[0-9][0-9]*\)*/version $(version)/" $(citation)
# @sed -i "s/20[0-9]\{2\}/$(yr)/" $(citation)


.PHONY: clean
clean:
	rm -rf *~ */*~ *.Rhistroy *.tar.gz *.Rcheck/ .\#*
