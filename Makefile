.PHONY:  test

repo=$(shell basename $(CURDIR))
client="mock"

test: 
	python -m coverage run -m pytest -s -x -k test_managed_data --client $(client)

	python -m coverage html
	-open htmlcov/index.html
	


# set up jupyter dev kernel
jupyter:
	-yes | jupyter kernelspec uninstall $(repo)
	python -m ipykernel install --user --name $(repo)