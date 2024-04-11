.PHONY:  test

repo=$(shell basename $(CURDIR))


test: 
	python -m coverage run -m pytest -s -k test_low_level 
	python -m coverage html
	-open htmlcov/index.html
	


# set up jupyter dev kernel
jupyter:
	-yes | jupyter kernelspec uninstall $(repo)
	python -m ipykernel install --user --name $(repo)