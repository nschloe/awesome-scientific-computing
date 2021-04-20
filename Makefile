default:
	@echo "\"make lint\"?"

lint:
	awesome-lint README.md
	awesome_bot --white-list https://awesome.re README.md
