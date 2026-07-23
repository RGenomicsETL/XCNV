.PHONY: rd readme test build check clean

rd:
	Rscript -e 'roxygen2::roxygenize(load_code = "source")'

readme:
	Rscript -e 'rmarkdown::render("README.Rmd", output_file = "README.md", quiet = TRUE)'

test:
	Rscript -e 'tinytest::build_install_test(".")'

build: rd readme
	mkdir -p ../xcnv-build
	(cd ../xcnv-build && R CMD build --no-manual $(CURDIR))

check: build
	R CMD check --no-manual ../xcnv-build/*.tar.gz

clean:
	rm -rf ../xcnv-build
