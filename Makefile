.PHONY: rd readme test build check clean

rd:
	Rscript -e 'roxygen2::roxygenize(load_code = "source")'

readme:
	Rscript -e 'lib <- tempfile("xcnv-readme-lib-"); dir.create(lib); on.exit(unlink(lib, recursive = TRUE), add = TRUE); utils::install.packages(".", repos = NULL, type = "source", lib = lib, quiet = TRUE); .libPaths(c(lib, .libPaths())); rmarkdown::render("README.Rmd", output_file = "README.md", quiet = TRUE)'

test:
	Rscript -e 'tinytest::build_install_test(".")'

build: rd readme
	mkdir -p ../xcnv-build
	(cd ../xcnv-build && R CMD build --no-manual $(CURDIR))

check: build
	R CMD check --no-manual ../xcnv-build/*.tar.gz

clean:
	rm -rf ../xcnv-build
