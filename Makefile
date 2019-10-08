default:
	@echo "\"make lint\"?"

lint:
	awesome-lint README.md
	awesome_bot --white-list https://awesome.re,https://portal.nersc.gov/project/sparse/superlu/ README.md
