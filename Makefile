.PHONY:  test

repo=$(shell basename $(CURDIR))


test: 
	#python -m coverage run -m pytest -s -x -k test_cli
	python -m coverage run -m pytest -s -x -k "test_config or test_feature_flags or test_managed_data"
	#python -m coverage run -m pytest -s -x -k test_context
	#python -m coverage run -m pytest -s -x -k test_variables

	python -m coverage html
	-open htmlcov/index.html
	


# set up jupyter dev kernel
jupyter:
	-yes | jupyter kernelspec uninstall $(repo)
	python -m ipykernel install --user --name $(repo)