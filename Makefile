default:
	@echo "\"make lint\"?"

lint:
	awesome-lint README.md
	# Allow redirect for awesome.re
	awesome_bot --allow-redirect README.md
