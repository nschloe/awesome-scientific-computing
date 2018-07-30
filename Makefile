default:
	@echo "\"make lint\"?"

lint:
	awesome-lint README.md
	# Allow SSL errors for Trilinos, allow redirect for awesome.re
	awesome_bot --allow-ssl --allow-redirect README.md
