.PHONY:  test

repo=$(shell basename $(CURDIR))

test:
	python -m pytest -s -k test_low_level


# set up jupyter dev kernel
jupyter:
	-yes | jupyter kernelspec uninstall $(repo)
	python -m ipykernel install --user --name $(repo)